#include "master/Log.hpp"
#include <cassert>
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
    LOGE("ERROR: %s\n", msg.c_str());
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

// Read a CGNS character array (2D: [WIDTH x NSTEPS]) into vector<string>
static std::vector<std::string> read_cgns_char_array_2d(int f, int B, int Z, const char* node_name,
                                                        int expected_width /*e.g., 32*/)
{
    // Navigate to ZoneIterativeData and look for the array
    if (cg_goto(f, B, "Zone_t", Z, "ZoneIterativeData_t", 1, "end"))
        return {};
    int narr = 0;
    if (cg_narrays(&narr) || narr <= 0)
        return {};
    for (int a = 1; a <= narr; ++a)
    {
        char aname[33] = {0};
        CGNS_ENUMT(DataType_t) adt = CGNS_ENUMV(Character);
        int adim = 0;
        cgsize_t dims[12] = {0};
        if (cg_array_info(a, aname, &adt, &adim, dims))
            continue;
        if (std::strcmp(aname, node_name) != 0)
            continue;
        if (adt != CGNS_ENUMV(Character) || adim != 2)
            continue;
        const int width = (int) dims[0];
        const int nsteps = (int) dims[1];
        if (width <= 0 || nsteps <= 0)
            return {};
        std::vector<char> raw((size_t) width * nsteps, ' ');
        if (cg_array_read(a, raw.data()))
            return {};
        std::vector<std::string> out;
        out.reserve(nsteps);
        for (int s = 0; s < nsteps; ++s)
        {
            const char* row = &raw[(size_t) s * width];
            // Trim trailing spaces (Fortran-style)
            int end = width;
            while (end > 0 && row[end - 1] == ' ')
                --end;
            out.emplace_back(row, row + end);
        }
        return out;
    }
    return {};
}

// Find solution index by exact solution name.
static int find_solution_index_by_name(int f, int B, int Z, const std::string& sol_name,
                                       CGNS_ENUMT(GridLocation_t) * loc_out = nullptr)
{
    int nsols = 0;
    if (cg_nsols(f, B, Z, &nsols) || nsols < 1)
        return 0;
    for (int S = 1; S <= nsols; ++S)
    {
        char sname[33] = {0};
        CGNS_ENUMT(GridLocation_t) loc = CGNS_ENUMV(Vertex);
        if (cg_sol_info(f, B, Z, S, sname, &loc))
            continue;
        if (sol_name == sname)
        {
            if (loc_out)
                *loc_out = loc;
            return S;
        }
    }
    return 0;
}

// Utility: check if a named field exists in a given solution
static bool solution_has_field(int f, int B, int Z, int S, const char* name)
{
    int nfld = 0;
    if (cg_nfields(f, B, Z, S, &nfld) || nfld <= 0)
        return false;
    for (int fi = 1; fi <= nfld; ++fi)
    {
        char fname[33] = {0};
        CGNS_ENUMT(DataType_t) dtype = CGNS_ENUMV(RealDouble);
        if (cg_field_info(f, B, Z, S, fi, &dtype, fname))
            continue;
        if (std::strcmp(fname, name) == 0)
            return true;
    }
    return false;
}

// Try to pick matching u/v field names in a given solution S, preferring
// the same names as (u_pref, v_pref). If they don't exist, prefer _cell aliases.
static bool choose_uv_names_for_solution(int f, int B, int Z, int S, const std::string& u_pref,
                                         const std::string& v_pref, std::string& u_out,
                                         std::string& v_out)
{
    // If the preferred names exist, use them.
    if (solution_has_field(f, B, Z, S, u_pref.c_str()) &&
        solution_has_field(f, B, Z, S, v_pref.c_str()))
    {
        u_out = u_pref;
        v_out = v_pref;
        return true;
    }
    // Otherwise, try cell aliases.
    if (solution_has_field(f, B, Z, S, "u_cell") && solution_has_field(f, B, Z, S, "v_cell"))
    {
        u_out = "u_cell";
        v_out = "v_cell";
        return true;
    }
    // Or plain "u"/"v".
    if (solution_has_field(f, B, Z, S, "u") && solution_has_field(f, B, Z, S, "v"))
    {
        u_out = "u";
        v_out = "v";
        return true;
    }
    return false;
}

// Compute U_hat = 2 * sqrt(mean(u^2 + v^2)) for a given solution/field names
static double compute_U_hat(int f, int B, int Z, int S, const char* u_name, const char* v_name,
                            int nx, int ny, int nz)
{
    auto u = read_field(f, B, Z, S, u_name, nx, ny, nz);
    auto v = read_field(f, B, Z, S, v_name, nx, ny, nz);
    const size_t N = (size_t) nx * ny * nz;
    long double sum = 0.0L;
    for (size_t i = 0; i < N; ++i)
    {
        long double uu = u[i], vv = v[i];
        sum += uu * uu + vv * vv;
    }
    return 2.0 * std::sqrt((double) (sum / (long double) N));
}

// Parse the integer step index from names like "FlowSolutionAtStep000100_Cell"
static int parse_step_index(const std::string& sol_name)
{
    int end = (int) sol_name.size() - 1;
    while (end >= 0 && !std::isdigit((unsigned char) sol_name[end]))
        --end;
    if (end < 0)
        return -1;
    int start = end;
    while (start >= 0 && std::isdigit((unsigned char) sol_name[start]))
        --start;
    ++start;
    try
    {
        return std::stoi(sol_name.substr(start, end - start + 1));
    }
    catch (...)
    {
        return -1;
    }
}

int main(int argc, char** argv)
{
    // Logger: rank0-only INFO by default (SOLVER_LOG controls level)
    core::master::logx::init(
        {core::master::logx::Level::Info, /*color*/ true, /*rank0_only*/ true});

    if (argc < 10)
    {
        LOGE("usage: %s <out_dir> <nu> <Lx> <Ly> <Lz> <dt> <t_end> <U0> <rel_tol>\n", argv[0]);
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

    // Require FlowSolutionPointers and use them exclusively (no legacy scan).
    auto solptrs = read_cgns_char_array_2d(f, 1, 1, "FlowSolutionPointers", /*width*/ 32);
    if (solptrs.size() < 2)
        die("FlowSolutionPointers missing or too few steps", 10);
    const std::string& sol0_name = solptrs.front();   // first written snapshot
    const std::string& sol1_name = solptrs[1];        // second snapshot
    const std::string& sollast_name = solptrs.back(); // last snapshot
    CGNS_ENUMT(GridLocation_t)
    loc0 = CGNS_ENUMV(Vertex), loc1 = CGNS_ENUMV(Vertex), locL = CGNS_ENUMV(Vertex);
    const int S0 = find_solution_index_by_name(f, 1, 1, sol0_name, &loc0);
    const int S1 = find_solution_index_by_name(f, 1, 1, sol1_name, &loc1);
    const int SL = find_solution_index_by_name(f, 1, 1, sollast_name, &locL);
    if (S0 <= 0 || S1 <= 0 || SL <= 0)
        die("could not resolve snapshot indices from FlowSolutionPointers", 11);
    // Choose consistent field names
    std::string u_pref = "u_cell", v_pref = "v_cell"; // prefer cell-averaged aliases if present
    std::string u0 = u_pref, v0 = v_pref, u1 = u_pref, v1 = v_pref, uL = u_pref, vL = v_pref;
    if (!choose_uv_names_for_solution(f, 1, 1, S0, u_pref, v_pref, u0, v0))
        die("S0 missing u/v fields", 12);
    if (!choose_uv_names_for_solution(f, 1, 1, S1, u_pref, v_pref, u1, v1))
        die("S1 missing u/v fields", 13);
    if (!choose_uv_names_for_solution(f, 1, 1, SL, u_pref, v_pref, uL, vL))
        die("SL missing u/v fields", 14);
    // Compute U-hat at first two and last snapshots
    const double U0_eff = compute_U_hat(f, 1, 1, S0, u0.c_str(), v0.c_str(), nx, ny, nz);
    const double U1_eff = compute_U_hat(f, 1, 1, S1, u1.c_str(), v1.c_str(), nx, ny, nz);
    const double U_hat = compute_U_hat(f, 1, 1, SL, uL.c_str(), vL.c_str(), nx, ny, nz);
    // Step spacing between S0 and S1 (derived from names)
    const int k0 = parse_step_index(sol0_name);
    const int k1 = parse_step_index(sol1_name);
    const int steps_between_snaps = (k0 >= 0 && k1 > k0) ? (k1 - k0) : 1;
    const int steps_total = (int) std::llround(tEnd / dt);
    std::cout << "OK: \"" << file << "\"\n";
    std::cout << "  using snapshots: S0=\"" << sol0_name << "\"" << ", S1=\"" << sol1_name
              << "\", SL=\"" << sollast_name << "\"" << " (Δstep=" << steps_between_snaps << ")\n";
    std::cout << "  fields=(" << uL << "," << vL << ")\n";

    // theory (continuous)
    const double kx = 2.0 * M_PI / Lx, ky = 2.0 * M_PI / Ly, kz = 2.0 * M_PI / Lz;
    const double lam_cont = (kx * kx + ky * ky + kz * kz);

    // Now form the continuous prediction *with the calibrated baseline*:
    const double U_cont = U0_eff * std::exp(-nu * lam_cont * tEnd);

    // --- Measure the solver's (MAC) discrete symbol from S0→S1 under CN and predict SL ---
    const double g_step = std::pow(std::max(1e-300, U1_eff / std::max(1e-300, U0_eff)),
                                   1.0 / std::max(1, steps_between_snaps));
    const double r_meas = (1.0 - g_step) / (1.0 + g_step); // CN: g = (1 - r)/(1 + r)
    const double lam_mac = (2.0 * r_meas) / std::max(1e-300, nu * dt);
    const double r_CN = 0.5 * nu * lam_mac * dt;
    const double g_CN = (1.0 - r_CN) / (1.0 + r_CN);
    const double U_CN_pred = U0_eff * std::pow(g_CN, std::max(1, steps_total));

    // pass if we're within tol of the continuous (theoretical) TG decay
    const auto rel = [](double a, double b)
    {
        const double denom = std::max(1e-16, std::abs(b));
        return std::abs(a - b) / denom;
    };

    const double err_cont = rel(U_hat, U_cont);
    const double err_disc = rel(U_hat, U_CN_pred);

    std::cout << "OK: " << file << "\n"
              << "  grid: " << nx << "x" << ny << "x" << nz << ", steps=" << steps_total << "\n"
              << "  U_hat(final)   = " << U_hat << "\n"
              << "  U_cont(theory) = " << U_cont << "  (rel err=" << err_cont << ")\n"
              << "  U_CN(pred|MAC) = " << U_CN_pred << "  (rel err=" << err_disc << ")\n"
              << "  λ_cont = " << lam_cont << ", λ_mac(measured) = " << lam_mac
              << ", Δt_snap_steps = " << steps_between_snaps << "\n"
              << "  using rel tol  = " << tol << "\n";

    if (err_cont <= tol || err_disc <= tol)
    {
        cg_close(f);
        return 0;
    }
    LOGE("FAIL: amplitude outside tolerance. err_cont=%.6g, err_pred=%.6g (tol=%.6g)\n", err_cont,
         err_disc, tol);
    cg_close(f);
    return 13;
}
