#include <cgnslib.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#ifndef CGNS_ENUMT
#define CGNS_ENUMT(t) t
#endif
#ifndef CGNS_ENUMV
#define CGNS_ENUMV(v) v
#endif

namespace fs = std::filesystem;

static void die(const std::string& msg, int code = 2)
{
    std::cerr << "ERROR: " << msg << "\n";
    std::exit(code);
}

static std::vector<double> read_field(int f, int B, int Z, int S, const char* name, int nx, int ny,
                                      int nz)
{
    std::vector<double> buf((size_t) nx * ny * nz);
    cgsize_t rmin[3] = {1, 1, 1};
    cgsize_t rmax[3] = {(cgsize_t) nx, (cgsize_t) ny, (cgsize_t) nz};
    if (cg_field_read(f, B, Z, S, name, CGNS_ENUMV(RealDouble), rmin, rmax, buf.data()))
        die(std::string("cg_field_read failed for ") + name, 12);
    return buf;
}

int main(int argc, char** argv)
{
    if (argc < 10)
    {
        std::cerr << "usage: " << argv[0]
                  << " <out_dir> <nu> <Lx> <Ly> <Lz> <dt> <t_end> <U0> <rel_tol>\n";
        return 2;
    }
    const fs::path out_dir = argv[1];
    const double nu = std::atof(argv[2]);
    const double Lx = std::atof(argv[3]);
    const double Ly = std::atof(argv[4]);
    const double Lz = std::atof(argv[5]);
    const double dt = std::atof(argv[6]);
    const double tEnd = std::atof(argv[7]);
    const double U0 = std::atof(argv[8]);
    const double tol = std::atof(argv[9]);

    // find the single .cgns (like existing checker)
    fs::path file;
    for (auto& e : fs::directory_iterator(out_dir))
    {
        if (e.is_regular_file() && e.path().extension() == ".cgns")
        {
            file = e.path();
            break;
        }
    }
    if (file.empty())
        die("No .cgns found in " + out_dir.string(), 3);

    int f = -1;
    if (cg_open(file.string().c_str(), CG_MODE_READ, &f))
        die("cg_open failed", 4);

    int nbases = 0;
    if (cg_nbases(f, &nbases) || nbases < 1)
        die("no bases", 5);
    int cell_dim = 0, phys_dim = 0;
    char bname[33] = {0};
    if (cg_base_read(f, 1, bname, &cell_dim, &phys_dim))
        die("cg_base_read", 6);
    if (cell_dim != 3 || phys_dim != 3)
        die("unexpected base dims", 7);

    int nzones = 0;
    if (cg_nzones(f, 1, &nzones) || nzones < 1)
        die("no zones", 8);
    char zname[33] = {0};
    cgsize_t size[9] = {0};
    if (cg_zone_read(f, 1, 1, zname, size))
        die("cg_zone_read", 9);

    // Structured; cells are size[3],size[4],size[5] in file order (I,J,K)
    const int nx = (int) size[3], ny = (int) size[4], nz = (int) size[5];
    int nsols = 0;
    if (cg_nsols(f, 1, 1, &nsols) || nsols < 1)
        die("no FlowSolution", 10);
    const int S = nsols; // take the last write as "final"

    // read u and v at final step
    auto u = read_field(f, 1, 1, S, "u", nx, ny, nz);
    auto v = read_field(f, 1, 1, S, "v", nx, ny, nz);
    cg_close(f);

    // mean(u^2 + v^2)
    const size_t N = (size_t) nx * ny * nz;
    long double accum = 0.0L;
    for (size_t i = 0; i < N; ++i)
    {
        const long double uu = u[i];
        const long double vv = v[i];
        accum += uu * uu + vv * vv;
    }
    const long double mean_u2v2 = accum / (long double) N;
    const double U_hat = 2.0 * std::sqrt((double) mean_u2v2);

    // theory (continuous)
    const double kx = 2.0 * M_PI / Lx, ky = 2.0 * M_PI / Ly, kz = 2.0 * M_PI / Lz;
    const double lam_cont = (kx * kx + ky * ky + kz * kz);
    const double U_cont = U0 * std::exp(-nu * lam_cont * tEnd);

    const int steps = (int) std::llround(tEnd / dt);

    // pass if we're within tol of the continuous (theoretical) TG decay
    const auto rel = [](double a, double b)
    {
        const double denom = std::max(1e-16, std::abs(b));
        return std::abs(a - b) / denom;
    };
    const double err_cont = rel(U_hat, U_cont);

    std::cout << "OK: " << file << "\n"
              << "  grid: " << nx << "x" << ny << "x" << nz << ", steps=" << steps << "\n"
              << "  U_hat(final)   = " << U_hat << "\n"
              << "  U_cont(theory) = " << U_cont << "  (rel err=" << err_cont << ")\n"
              << "  using rel tol  = " << tol << "\n";

    if (err_cont > tol)
    {
        std::cerr << "FAIL: amplitude outside tolerance (rel err=" << err_cont << " > " << tol
                  << ")\n";
        return 13;
    }
    return 0;
}
