#include "simple_bench.hpp"
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

// PETSc/MPI lifecycle
#include "test_petsc_guard.hpp"
static PetscTestGuard _petsc_guard; // ensures PETSc (and MPI if available) is (de)initialized

#include "PressurePoisson.hpp"   // fluids::make_poisson
#include "Program.hpp"           // KV, BcTable consumed inside make_poisson
#include "master/FieldCatalog.hpp"
#include "master/RunContext.hpp"
#include "master/Views.hpp"

// Byte-stride helper
static inline std::array<std::ptrdiff_t,3> strides_bytes(int nx, int ny) {
    const std::ptrdiff_t s0 = (std::ptrdiff_t)sizeof(double);
    const std::ptrdiff_t s1 = (std::ptrdiff_t)nx * s0;
    const std::ptrdiff_t s2 = (std::ptrdiff_t)nx * (std::ptrdiff_t)ny * s0;
    return {s0,s1,s2};
}

// Build MAC faces so (1/dt) div(u*) = div(∇p*) with β=1
static void build_faces_from_pstar(const std::vector<double>& pstar,
                                   std::vector<double>& u, std::vector<double>& v, std::vector<double>& w,
                                   int nx, int ny, int nz, int ng, double dx, double dy, double dz, double dt)
{
    const int nxc_tot = nx + 2*ng, nyc_tot = ny + 2*ng, nzc_tot = nz + 2*ng;
    const int nxu_tot = nxc_tot + 1, nyu_tot = nyc_tot,     nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot,     nyv_tot = nyc_tot + 1, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot,     nyw_tot = nyc_tot,     nzw_tot = nzc_tot + 1;
    auto cidx = [&](int I,int J,int K){ return (size_t)I + (size_t)nxc_tot * ((size_t)J + (size_t)nyc_tot * (size_t)K); };
    auto aidx = [&](int I,int J,int K,int nxT,int nyT){ return (size_t)I + (size_t)nxT * ((size_t)J + (size_t)nyT * (size_t)K); };
    // u(i+1/2,j,k)
    for(int K=0;K<nzu_tot;++K) for(int J=0;J<nyu_tot;++J) for(int I=0;I<nxu_tot;++I){
        if(I==0||I==nxu_tot-1){ u[aidx(I,J,K,nxu_tot,nyu_tot)] = 0.0; continue; }
        const int ic=I-1,jc=J,kc=K;
        const double pc = pstar[cidx(ic,   jc,kc)];
        const double pe = pstar[cidx(ic+1, jc,kc)];
        u[aidx(I,J,K,nxu_tot,nyu_tot)] = -(dt/dx)*(pe-pc);
    }
    // v(i,j+1/2,k)
    for(int K=0;K<nzv_tot;++K) for(int J=0;J<nyv_tot;++J) for(int I=0;I<nxv_tot;++I){
        if(J==0||J==nyv_tot-1){ v[aidx(I,J,K,nxv_tot,nyv_tot)] = 0.0; continue; }
        const int ic=I,jc=J-1,kc=K;
        const double ps = pstar[cidx(ic,jc,   kc)];
        const double pn = pstar[cidx(ic,jc+1, kc)];
        v[aidx(I,J,K,nxv_tot,nyv_tot)] = -(dt/dy)*(pn-ps);
    }
    // w(i,j,k+1/2)
    for(int K=0;K<nzw_tot;++K) for(int J=0;J<nyw_tot;++J) for(int I=0;I<nxw_tot;++I){
        if(K==0||K==nzw_tot-1){ w[aidx(I,J,K,nxw_tot,nyw_tot)] = 0.0; continue; }
        const int ic=I,jc=J,kc=K-1;
        const double pb = pstar[cidx(ic,jc,kc   )];
        const double pt = pstar[cidx(ic,jc,kc+1 )];
        w[aidx(I,J,K,nxw_tot,nyw_tot)] = -(dt/dz)*(pt-pb);
    }
}

int main() {
    // Problem size & geometry
    const int nx=64, ny=64, nz=64, ng=1;
    const int nxc_tot = nx + 2*ng, nyc_tot = ny + 2*ng, nzc_tot = nz + 2*ng;
    const int nxu_tot = nxc_tot + 1, nyu_tot = nyc_tot,     nzu_tot = nzc_tot;
    const int nxv_tot = nxc_tot,     nyv_tot = nyc_tot + 1, nzv_tot = nzc_tot;
    const int nxw_tot = nxc_tot,     nyw_tot = nyc_tot,     nzw_tot = nzc_tot + 1;
    const double dx=1.0, dy=1.0, dz=1.0, dt=1.0;

    // Allocate fields like the tests
    std::vector<double> p   ((size_t)nxc_tot*nyc_tot*nzc_tot, 0.0);
    std::vector<double> rho ((size_t)nxc_tot*nyc_tot*nzc_tot, 1.0);
    std::vector<double> u   ((size_t)nxu_tot*nyu_tot*nzu_tot, 0.0);
    std::vector<double> v   ((size_t)nxv_tot*nyv_tot*nzv_tot, 0.0);
    std::vector<double> w   ((size_t)nxw_tot*nyw_tot*nzw_tot, 0.0);
    std::vector<double> pstar = p; // manufactured p*

    // Manufactured p*
    for (int K=0; K<nzc_tot; ++K)
      for (int J=0; J<nyc_tot; ++J)
        for (int I=0; I<nxc_tot; ++I)
          pstar[(size_t)I + (size_t)nxc_tot*((size_t)J + (size_t)nyc_tot*(size_t)K)]
            = std::sin(2*M_PI*(I-0.5)/nx) + 0.3*std::cos(2*M_PI*(J-0.5)/ny) + 0.2*std::cos(2*M_PI*(K-0.5)/nz);

    // Build MAC faces so (1/dt)div(u*) = div(∇p*) with β=1
    build_faces_from_pstar(pstar, u, v, w, nx, ny, nz, ng, dx, dy, dz, dt);

    // One-tile view like the tests
    core::master::MeshTileView tile{};
    tile.box.lo[0]=0; tile.box.lo[1]=0; tile.box.lo[2]=0;
    tile.box.hi[0]=nx; tile.box.hi[1]=ny; tile.box.hi[2]=nz;

    core::master::FieldCatalog fields;
    fields.register_scalar("p",   p.data(),   sizeof(double), {nxc_tot,nyc_tot,nzc_tot}, strides_bytes(nxc_tot,nyc_tot));
    fields.register_scalar("u",   u.data(),   sizeof(double), {nxu_tot,nyu_tot,nzu_tot}, strides_bytes(nxu_tot,nyu_tot));
    fields.register_scalar("v",   v.data(),   sizeof(double), {nxv_tot,nyv_tot,nzv_tot}, strides_bytes(nxv_tot,nyv_tot));
    fields.register_scalar("w",   w.data(),   sizeof(double), {nxw_tot,nyw_tot,nzw_tot}, strides_bytes(nxw_tot,nyw_tot));
    fields.register_scalar("rho", rho.data(), sizeof(double), {nxc_tot,nyc_tot,nzc_tot}, strides_bytes(nxc_tot,nyc_tot));

    // Build the MG Poisson actio
    core::master::plugin::KV kv{
        {"dx", std::to_string(dx)}, {"dy", std::to_string(dy)}, {"dz", std::to_string(dz)},
        {"rho","1.0"}, {"iters","50"}, {"div_tol","1e-10"},
        // Dirichlet on all sides to pin the gauge
        {"p.west","dirichlet:0"},{"p.east","dirichlet:0"},
        {"p.south","dirichlet:0"},{"p.north","dirichlet:0"},
        {"p.bottom","dirichlet:0"},{"p.top","dirichlet:0"}
    };
    core::master::RunContext rc{};
    if (PetscTestGuard::petsc_uses_mpi()) {
        rc.mpi_comm = const_cast<void*>(_petsc_guard.mpi_comm_ptr()); // else leave nullptr
    }
    auto poisson = fluids::make_poisson(kv, rc);

    // Time a full action execute (this runs MG inside PressurePoisson)
    auto [mean_ms, std_ms] = bench::run([&]{ poisson->execute(tile, fields, dt); });
    bench::report("fluids_poisson_mg_fullpath_64^3", mean_ms, std_ms, /*bytes=*/0.0);
    return 0;
}