#include "MacOps.hpp"
#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>
using namespace numerics::kernels;

using Catch::Approx;

static inline int center_tot(int n, int ng)
{
    return n + 2 * ng;
}
static inline int face_tot(int n, int ng)
{
    return n + 1 + 2 * ng;
}
static inline std::size_t idx(int I, int J, int K, int nx, int ny)
{
    return std::size_t(I) + std::size_t(nx) * (std::size_t(J) + std::size_t(ny) * std::size_t(K));
}

TEST_CASE("Diffusion (FE, MAC) preserves constant velocity field on interior faces",
          "[fluids][diffuse][mac]")
{
    const int nx = 12, ny = 7, nz = 5, ng = 2;
    const double dx = 1.0, dy = 1.0, dz = 1.0, dt = 0.3;

    // MAC totals
    const int nxc_tot = center_tot(nx, ng), nyc_tot = center_tot(ny, ng),
              nzc_tot = center_tot(nz, ng);
    const int nxu_tot = face_tot(nx, ng), nyu_tot = center_tot(ny, ng),
              nzu_tot = center_tot(nz, ng);
    const int nxv_tot = center_tot(nx, ng), nyv_tot = face_tot(ny, ng),
              nzv_tot = center_tot(nz, ng);
    const int nxw_tot = center_tot(nx, ng), nyw_tot = center_tot(ny, ng),
              nzw_tot = face_tot(nz, ng);

    // Allocate faces and centers
    std::vector<double> u((size_t) nxu_tot * nyu_tot * nzu_tot, 3.14);
    std::vector<double> v((size_t) nxv_tot * nyv_tot * nzv_tot, 3.14);
    std::vector<double> w((size_t) nxw_tot * nyw_tot * nzw_tot, 3.14);
    std::vector<double> nu_eff((size_t) nxc_tot * nyc_tot * nzc_tot, 1.0e-3);
    std::vector<double> us(u.size(), 0.0), vs(v.size(), 0.0), ws(w.size(), 0.0);

    diffuse_fe(
        /*u in*/ u.data(), nxu_tot, nyu_tot, nzu_tot,
        /*v in*/ v.data(), nxv_tot, nyv_tot, nzv_tot,
        /*w in*/ w.data(), nxw_tot, nyw_tot, nzw_tot,
        /*nu  */ nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
        /*geo */ ng, dx, dy, dz, dt,
        /*out */ us.data(), vs.data(), ws.data());

    // Check interior faces only (ghosts belong to BC/exchange)
    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + 1 + ng; ++I) // u is x-face ⇒ +1 in x
                REQUIRE(us[idx(I, J, K, nxu_tot, nyu_tot)] == Approx(3.14));

    for (int K = ng; K < nz + ng; ++K)
        for (int J = ng; J < ny + 1 + ng; ++J) // v is y-face ⇒ +1 in y
            for (int I = ng; I < nx + ng; ++I)
                REQUIRE(vs[idx(I, J, K, nxv_tot, nyv_tot)] == Approx(3.14));

    for (int K = ng; K < nz + 1 + ng; ++K) // w is z-face ⇒ +1 in z
        for (int J = ng; J < ny + ng; ++J)
            for (int I = ng; I < nx + ng; ++I)
                REQUIRE(ws[idx(I, J, K, nxw_tot, nyw_tot)] == Approx(3.14));
}
