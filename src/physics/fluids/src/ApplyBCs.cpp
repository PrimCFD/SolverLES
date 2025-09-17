#include "ApplyBCs.hpp"
#include "master/FieldCatalog.hpp"
#include "mesh/Boundary.hpp" // Array3DView + apply_* kernels
#include "mesh/Field.hpp"    // core typed view
#include <memory>
#include <string>

using namespace core::master;
using namespace core::master::plugin;
using core::mesh::apply_mirror_vector;
using core::mesh::apply_scalar_bc;
using core::mesh::Array3DView;
using core::mesh::Axis;
using core::mesh::BCOp;
using core::mesh::Field;
using core::mesh::MirrorMask;

namespace fluids
{

// ------------ helpers -------------
static inline const BcSpec* find_bc(const BcTable& b, const std::string& key)
{
    auto it = b.find(key);
    return (it == b.end() ? nullptr : &it->second);
}
static inline BCOp to_op(BcSpec::Type t)
{
    switch (t)
    {
    case BcSpec::Type::dirichlet:
        return BCOp::Dirichlet;
    case BcSpec::Type::neumann:
        return BCOp::NeumannZero;
    case BcSpec::Type::extrap:
        return BCOp::Extrapolate1;
    case BcSpec::Type::mirror:
        return BCOp::Mirror;
    }
    return BCOp::Dirichlet;
}
static inline MirrorMask mask_for(Axis ax)
{
    switch (ax)
    {
    case Axis::I:
        return MirrorMask{{-1, +1, +1}}; // flip u_n
    case Axis::J:
        return MirrorMask{{+1, -1, +1}}; // flip v_n
    case Axis::K:
        return MirrorMask{{+1, +1, -1}}; // flip w_n
    }
    return MirrorMask{{+1, +1, +1}};
}

// ------------ action wiring -------------
ApplyBCs::ApplyBCs(const Params& p) : bcs_(p.bcs)
{
    info_.name = "apply_bcs";
    info_.phases = plugin::Phase::PreExchange | plugin::Phase::Interior | plugin::Phase::PostBC;
}

std::shared_ptr<IAction> make_apply_bcs(const KV& kv)
{
    return std::make_shared<ApplyBCs>(parse_params(kv));
}

void ApplyBCs::execute(const MeshTileView& tile, FieldCatalog& fields, double)
{

    // Pull one view just to compute totals & ng (uniform ghosts)
    auto pick_any = [&]()
    {
        if (fields.contains("u"))
            return fields.view("u");
        if (fields.contains("v"))
            return fields.view("v");
        if (fields.contains("w"))
            return fields.view("w");
        return fields.view("p");
    };
    auto any = pick_any();

    const int nx = tile.box.hi[0] - tile.box.lo[0];
    const int ny = tile.box.hi[1] - tile.box.lo[1];
    const int nz = tile.box.hi[2] - tile.box.lo[2];
    const int nx_tot = any.extents[0];
    const int ny_tot = any.extents[1];
    const int nz_tot = any.extents[2];
    const int ng = (nx_tot - nx) / 2; // uniform halo width

    // Build core Field<T> wrappers (typed, non-owning) + Array3DView<T>
    std::unique_ptr<Field<double>> fu, fv, fw, fp;
    Array3DView<double> U, V, W, P;

    if (fields.contains("u"))
    {
        auto vu = fields.view("u");
        fu = std::make_unique<Field<double>>(
            static_cast<double*>(vu.host_ptr),
            std::array<int, 3>{vu.extents[0], vu.extents[1], vu.extents[2]}, ng);
        U = Array3DView<double>(*fu);
    }
    if (fields.contains("v"))
    {
        auto vv = fields.view("v");
        fv = std::make_unique<Field<double>>(
            static_cast<double*>(vv.host_ptr),
            std::array<int, 3>{vv.extents[0], vv.extents[1], vv.extents[2]}, ng);
        V = Array3DView<double>(*fv);
    }
    if (fields.contains("w"))
    {
        auto vw = fields.view("w");
        fw = std::make_unique<Field<double>>(
            static_cast<double*>(vw.host_ptr),
            std::array<int, 3>{vw.extents[0], vw.extents[1], vw.extents[2]}, ng);
        W = Array3DView<double>(*fw);
    }
    if (fields.contains("p"))
    {
        auto vp = fields.view("p");
        fp = std::make_unique<Field<double>>(
            static_cast<double*>(vp.host_ptr),
            std::array<int, 3>{vp.extents[0], vp.extents[1], vp.extents[2]}, ng);
        P = Array3DView<double>(*fp);
    }

    struct Face
    {
        const char* key;
        Axis ax;
        int sgn;
    };
    static constexpr Face faces[6] = {{"west", Axis::I, -1},   {"east", Axis::I, +1},
                                      {"south", Axis::J, -1},  {"north", Axis::J, +1},
                                      {"bottom", Axis::K, -1}, {"top", Axis::K, +1}};

    auto have = [&](const char* n) { return fields.contains(n); };
    auto is_mirror = [](const BcSpec* s) { return s && s->type == BcSpec::Type::mirror; };

    for (const auto& f : faces)
    {
        const BcSpec* su = find_bc(bcs_, std::string("u.") + f.key);
        const BcSpec* sv = find_bc(bcs_, std::string("v.") + f.key);
        const BcSpec* sw = find_bc(bcs_, std::string("w.") + f.key);
        const BcSpec* sp = find_bc(bcs_, std::string("p.") + f.key);

        // Vector mirror if requested on all 3 velocity components and present
        if (have("u") && have("v") && have("w") && is_mirror(su) && is_mirror(sv) && is_mirror(sw))
        {
            apply_mirror_vector(U, V, W, f.ax, f.sgn, mask_for(f.ax));
        }
        else
        {
            // Component-wise scalar ops
            if (have("u") && su)
                apply_scalar_bc(U, f.ax, f.sgn, to_op(su->type), su->value);
            if (have("v") && sv)
                apply_scalar_bc(V, f.ax, f.sgn, to_op(sv->type), sv->value);
            if (have("w") && sw)
                apply_scalar_bc(W, f.ax, f.sgn, to_op(sw->type), sw->value);
        }

        // Pressure: treat mirror as NeumannZero (scalar mirror is ignored)
        if (have("p") && sp)
        {
            BCOp op = to_op(sp->type);
            if (op == BCOp::Mirror)
                op = BCOp::NeumannZero;
            apply_scalar_bc(P, f.ax, f.sgn, op, sp->value);
        }
    }
}

} // namespace fluids
