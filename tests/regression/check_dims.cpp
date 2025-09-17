#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <cgnslib.h>

// ---- CGNS compatibility shims ----
#ifndef CGNS_ENUMT
#define CGNS_ENUMT(t) t
#endif
#ifndef CGNS_ENUMV
#define CGNS_ENUMV(v) v
#endif

namespace fs = std::filesystem;

struct Dims
{
    int a{}, b{}, c{};
};

static bool read_file(const std::string& path, std::string& out)
{
    std::ifstream in(path);
    if (!in)
        return false;
    out.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    return true;
}

static std::string lower(std::string s)
{
    for (auto& c : s)
        c = (char) std::tolower((unsigned char) c);
    return s;
}

static bool match_either_order(const Dims& got, int x, int y, int z)
{
    const bool kji = (got.a == z && got.b == y && got.c == x);
    const bool ijk = (got.a == x && got.b == y && got.c == z);
    return kji || ijk;
}

/* ---------------- XDMF helpers (existing behavior) ---------------- */
static bool parse_topology_kji(const std::string& xml, Dims& topo, std::string& topo_type)
{
    std::regex re(
        R"(<Topology[^>]*TopologyType\s*=\s*\"([^\"]+)\"[^>]*Dimensions\s*=\s*\"(\d+)\s+(\d+)\s+(\d+)\")",
        std::regex::icase);
    std::smatch m;
    if (!std::regex_search(xml, m, re))
        return false;
    topo_type = m[1].str();
    topo.a = std::stoi(m[2]);
    topo.b = std::stoi(m[3]);
    topo.c = std::stoi(m[4]);
    return true;
}
struct Attr
{
    std::string name, center;
    Dims data;
};
static std::vector<Attr> parse_attributes(const std::string& xml)
{
    std::regex re_attr(
        R"(<Attribute[^>]*Name\s*=\s*\"([^\"]+)\"[^>]*Center\s*=\s*\"(Cell|Node)\"[^>]*>([\s\S]*?)</Attribute>)",
        std::regex::icase);
    std::regex re_data(R"(<DataItem[^>]*Dimensions\s*=\s*\"(\d+)\s+(\d+)\s+(\d+)\"[^>]*>)",
                       std::regex::icase);
    std::vector<Attr> out;
    for (auto it = std::sregex_iterator(xml.begin(), xml.end(), re_attr);
         it != std::sregex_iterator(); ++it)
    {
        Attr a;
        a.name = (*it)[1].str();
        a.center = (*it)[2].str();
        std::smatch m;
        const std::string inner = (*it)[3].str();
        if (std::regex_search(inner, m, re_data))
        {
            a.data.a = std::stoi(m[1]);
            a.data.b = std::stoi(m[2]);
            a.data.c = std::stoi(m[3]);
            out.push_back(a);
        }
    }
    return out;
}

/* ---------------- CGNS checker ---------------- */
static int check_cgns_file(const fs::path& file, int nx, int ny, int nz)
{
    int f = -1;
    if (cg_open(file.string().c_str(), CG_MODE_READ, &f))
    {
        std::cerr << "CGNS: could not open " << file << "\n";
        return 20;
    }

    int nbases = 0;
    if (cg_nbases(f, &nbases) || nbases < 1)
    {
        std::cerr << "CGNS: no bases\n";
        cg_close(f);
        return 21;
    }

    int cell_dim = 0, phys_dim = 0;
    char bname[33] = {0};
    if (cg_base_read(f, 1, bname, &cell_dim, &phys_dim) || cell_dim != 3 || phys_dim != 3)
    {
        std::cerr << "CGNS: unexpected base dims (cell_dim=" << cell_dim
                  << ", phys_dim=" << phys_dim << ")\n";
        cg_close(f);
        return 22;
    }

    int nzones = 0;
    if (cg_nzones(f, 1, &nzones) || nzones < 1)
    {
        std::cerr << "CGNS: no zones in base\n";
        cg_close(f);
        return 23;
    }

    char zname[33] = {0};
    cgsize_t size[9] = {0};
    if (cg_zone_read(f, 1, 1, zname, size))
    {
        std::cerr << "CGNS: cg_zone_read failed\n";
        cg_close(f);
        return 24;
    }
    CGNS_ENUMT(ZoneType_t) zt{};
    if (cg_zone_type(f, 1, 1, &zt) || zt != CGNS_ENUMV(Structured))
    {
        std::cerr << "CGNS: zone is not Structured\n";
        cg_close(f);
        return 25;
    }

    const Dims vertices{(int) size[0], (int) size[1], (int) size[2]};
    const Dims cells{(int) size[3], (int) size[4], (int) size[5]};

    const int ex_nx_nodes = nx + 1, ex_ny_nodes = ny + 1, ex_nz_nodes = nz + 1;
    if (!match_either_order(vertices, ex_nx_nodes, ex_ny_nodes, ex_nz_nodes))
    {
        std::cerr << "CGNS: vertex sizes mismatch. Got (" << vertices.a << "," << vertices.b << ","
                  << vertices.c << "), expected (" << ex_nx_nodes << "," << ex_ny_nodes << ","
                  << ex_nz_nodes << ") in any order.\n";
        cg_close(f);
        return 26;
    }
    if (!match_either_order(cells, nx, ny, nz))
    {
        std::cerr << "CGNS: cell sizes mismatch. Got (" << cells.a << "," << cells.b << ","
                  << cells.c << "), expected (" << nx << "," << ny << "," << nz
                  << ") in any order.\n";
        cg_close(f);
        return 27;
    }

    int ncoords = 0;
    if (cg_ncoords(f, 1, 1, &ncoords) || ncoords < 3)
    {
        std::cerr << "CGNS: expected at least 3 coordinates (X,Y,Z), got " << ncoords << "\n";
        cg_close(f);
        return 28;
    }

    int nsols = 0;
    if (cg_nsols(f, 1, 1, &nsols) || nsols == 0)
    {
        std::cerr << "CGNS: no FlowSolution nodes found\n";
        cg_close(f);
        return 29;
    }
    int bad_loc = 0, empty_sol = 0;
    for (int is = 1; is <= nsols; ++is)
    {
        char sname[33] = {0};
        CGNS_ENUMT(GridLocation_t) loc{};
        if (cg_sol_info(f, 1, 1, is, sname, &loc))
            continue;
        if (loc != CGNS_ENUMV(CellCenter))
        {
            std::cerr << "CGNS: FlowSolution \"" << sname << "\" not at CellCenter\n";
            ++bad_loc;
        }
        int nfld = 0;
        if (cg_nfields(f, 1, 1, is, &nfld) || nfld <= 0)
        {
            std::cerr << "CGNS: FlowSolution \"" << sname << "\" has no fields\n";
            ++empty_sol;
        }
    }
    cg_close(f);
    if (bad_loc)
        return 30;
    if (empty_sol)
        return 31;

    std::cout << "OK: " << file << "\n"
              << "  Zone (cells): " << cells.a << "," << cells.b << "," << cells.c << "\n"
              << "  Verified " << ncoords << " coordinate array(s), " << nsols
              << " FlowSolution(s) at CellCenter with fields.\n";
    return 0;
}

/* ---------------- main: pick .xmf or .cgns ---------------- */
int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "usage: " << argv[0] << " <out_dir> <nx> <ny> <nz>\n";
        return 2;
    }
    const fs::path out_dir = argv[1];
    const int nx = std::stoi(argv[2]);
    const int ny = std::stoi(argv[3]);
    const int nz = std::stoi(argv[4]);

    if (!fs::exists(out_dir))
    {
        std::cerr << "Missing output directory: " << out_dir << "\n";
        return 3;
    }

    std::vector<fs::path> outs;
    for (auto& e : fs::directory_iterator(out_dir))
    {
        if (!e.is_regular_file())
            continue;
        const auto ext = lower(e.path().extension().string());
        if (ext == ".xmf" || ext == ".cgns")
            outs.push_back(e.path());
    }

    if (outs.empty())
    {
        std::cerr << "Missing .xmf or .cgns in: " << out_dir << "\n";
        return 4;
    }
    if (outs.size() > 1)
    {
        std::cerr << "Found multiple output files in " << out_dir << ":\n";
        for (auto& p : outs)
            std::cerr << "  - " << p << "\n";
        return 5;
    }

    const fs::path file = outs.front();
    const std::string ext = lower(file.extension().string());

    if (ext == ".cgns")
    {
        return check_cgns_file(file, nx, ny, nz);
    }

    // ----- XDMF path -----
    std::string xml;
    if (!read_file(file.string(), xml))
    {
        std::cerr << "Could not read " << file << "\n";
        return 6;
    }

    Dims topo{};
    std::string topo_type;
    if (!parse_topology_kji(xml, topo, topo_type))
    {
        std::cerr << "Could not find <Topology ... Dimensions=\"K J I\"> in " << file << "\n";
        return 7;
    }
    if (!std::regex_search(topo_type, std::regex("3DCoRectMesh", std::regex::icase)))
    {
        std::cerr << "Unexpected TopologyType=\"" << topo_type << "\" (expected 3DCoRectMesh)\n";
        return 8;
    }

    const int ex_nx_nodes = nx + 1, ex_ny_nodes = ny + 1, ex_nz_nodes = nz + 1;
    if (!match_either_order(topo, ex_nx_nodes, ex_ny_nodes, ex_nz_nodes))
    {
        std::cerr << "Topology mismatch. Parsed Dimensions=(" << topo.a << "," << topo.b << ","
                  << topo.c << ")\n"
                  << "Expected nodes (+1): (" << ex_nz_nodes << "," << ex_ny_nodes << ","
                  << ex_nx_nodes << ") or (" << ex_nx_nodes << "," << ex_ny_nodes << ","
                  << ex_nz_nodes << ")\n";
        return 9;
    }

    auto attrs = parse_attributes(xml);
    if (attrs.empty())
    {
        std::cerr << "No <Attribute> blocks found in " << file << "\n";
        return 10;
    }
    int errors = 0;
    for (const auto& a : attrs)
    {
        if (std::regex_search(a.center, std::regex("Cell", std::regex::icase)))
        {
            if (!match_either_order(a.data, nx, ny, nz))
            {
                std::cerr << "Attribute \"" << a.name << "\" Center=Cell has Dimensions=("
                          << a.data.a << "," << a.data.b << "," << a.data.c << ") "
                          << "but expected (" << nz << "," << ny << "," << nx << ") or (" << nx
                          << "," << ny << "," << nz << ")\n";
                ++errors;
            }
        }
        else if (std::regex_search(a.center, std::regex("Node", std::regex::icase)))
        {
            const int ex_nx_nodes2 = nx + 1, ex_ny_nodes2 = ny + 1, ex_nz_nodes2 = nz + 1;
            if (!match_either_order(a.data, ex_nx_nodes2, ex_ny_nodes2, ex_nz_nodes2))
            {
                std::cerr << "Attribute \"" << a.name << "\" Center=Node has Dimensions=("
                          << a.data.a << "," << a.data.b << "," << a.data.c << ") "
                          << "but expected nodes (+1): (" << ex_nz_nodes2 << "," << ex_ny_nodes2
                          << "," << ex_nx_nodes2 << ") or (" << ex_nx_nodes2 << "," << ex_ny_nodes2
                          << "," << ex_nz_nodes2 << ")\n";
                ++errors;
            }
        }
        else
        {
            std::cerr << "Attribute \"" << a.name << "\" has unknown Center=\"" << a.center
                      << "\"\n";
            ++errors;
        }
    }
    if (errors)
        return 11;

    std::cout << "OK: " << file << "\n"
              << "  Topology (nodes): " << topo.a << "," << topo.b << "," << topo.c << "\n"
              << "  Checked " << attrs.size() << " attribute(s) with Option A rules.\n";
    return 0;
}
