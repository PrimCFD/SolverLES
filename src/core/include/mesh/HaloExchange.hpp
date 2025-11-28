#pragma once
#include "Field.hpp"
#include "Mesh.hpp"
#include "master/Views.hpp"

#include <mpi.h>
#include <array>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <assert.h>

namespace core::mesh
{

// 26-neighbor halo exchange with MAC-staggering support.
// - T in {float,double}
// - comm must be a 3D Cartesian communicator whose periodicity matches Mesh::periodic
// - tag_base must reserve a window >= 26 tags (HaloOps reserves 64 per field)
//
// Tag contract:
//   For each direction d with index id = dir_index(ox,oy,oz),
//   and opposite direction d' with opp_id = dir_index(-ox,-oy,-oz):
//     recv tag = tag_base + id      (expecting data coming *from* direction d)
//     send tag = tag_base + opp_id  (so our neighbor’s recv tag matches)
//
// That way rank A sending east uses the same tag as rank B receiving west, etc.
template <class T>
bool exchange_ghosts(Field<T>& f,
                     const Mesh& m,
                     MPI_Comm comm,
                     int tag_base = 200,
                     core::master::Stagger stagger = core::master::Stagger::Cell)
{
    static_assert(std::is_same_v<T,float> || std::is_same_v<T,double>,
                  "exchange_ghosts<T>: T must be float or double");

    // Use the mesh ghost width (or f.ng(); they should match)
    const int ng = m.ng;
    const auto e  = f.extents();

    // Interior extents (cell or face); for faces this is N_cells_local+1 on the normal axis.
    const int nxI = e[0] - 2*ng;
    const int nyI = e[1] - 2*ng;
    const int nzI = e[2] - 2*ng;

    if (ng <= 0 || nxI <= 0 || nyI <= 0 || nzI <= 0) return true;

    // MAC staggering flags: used only to pick which interior slab to send
    // along the normal direction; they do NOT change the interior size.
    const bool faceI = (stagger == core::master::Stagger::IFace);
    const bool faceJ = (stagger == core::master::Stagger::JFace);
    const bool faceK = (stagger == core::master::Stagger::KFace);

    int init = 0;
    MPI_Initialized(&init);
    if (!init || comm == MPI_COMM_NULL) return false;
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);

    int topo = MPI_UNDEFINED;
    MPI_Topo_test(comm, &topo);
    if (topo != MPI_CART) return false;

    int cdims[3]{}, cper[3]{}, ccoords[3]{};
    MPI_Cart_get(comm, 3, cdims, cper, ccoords);
    if ((cper[0] != (m.periodic[0]?1:0)) ||
        (cper[1] != (m.periodic[1]?1:0)) ||
        (cper[2] != (m.periodic[2]?1:0))) return false;

    // --- Tag safety: ensure tag_base..tag_base+25 fit under MPI_TAG_UB ---
    {
        int flag = 0;
        int *p_ub = nullptr;
        MPI_Comm_get_attr(comm, MPI_TAG_UB, &p_ub, &flag);
        if (flag && p_ub) {
            const int tag_ub = *p_ub;
            const int max_tag_needed = tag_base + 25;
            if (max_tag_needed > tag_ub) {
                return false;
            }
        }
        // If !flag, assume tag_base was chosen sanely by the caller.
    }

    int myrank = -1;
    MPI_Comm_rank(comm, &myrank);
    MPI_Datatype Tmpi = std::is_same_v<T,double> ? MPI_DOUBLE : MPI_FLOAT;

    auto rank_of = [&](int ox,int oy,int oz)->int {
        int c[3] = { ccoords[0]+ox, ccoords[1]+oy, ccoords[2]+oz };
        for (int a=0; a<3; ++a){
            if (c[a] < 0 || c[a] >= cdims[a]){
                if (cper[a]) {
                    if (c[a] < 0)                c[a] += cdims[a];
                    else if (c[a] >= cdims[a])   c[a] -= cdims[a];
                } else {
                    return MPI_PROC_NULL;
                }
            }
        }
        int r; MPI_Cart_rank(comm, c, &r); return r;
    };

    // Map (ox,oy,oz) ∈ {-1,0,1}^3 \ {(0,0,0)} → [0..25],
    // paired with opposite(ox,oy,oz)
    auto dir_index = [](int ox,int oy,int oz)->int {
        int ix = ox+1, iy = oy+1, iz = oz+1;   // each in {0,1,2}
        int id = ix*9 + iy*3 + iz;            // 0..26
        return (id>13 ? id-1 : id);           // skip center (id==13)
    };
    auto opposite = [](int ox,int oy,int oz){ return std::array<int,3>{-ox,-oy,-oz}; };

    struct Box { int i0,j0,k0, sx,sy,sz; }; // start + size (in interior indexing)

    // n  = interior extent along this axis (cells or faces)
    // g  = number of ghost layers
    // o  = neighbor offset along this axis (-1,0,+1)
    // isFaceNormal = true if this axis is the MAC-normal one for the field
    auto src_start_axis = [&](int n, int g, int o, bool isFaceNormal)->int
    {
        if (o == 0) return 0;

        // Cell-centred (or tangential-to-axis) case: original behaviour.
        if (!isFaceNormal)
            return (o < 0 ? 0 : (n - g));

        // Face-centred in the normal direction:
        // here n = N_cells_local + 1.  Let nc be the number of cells.
        const int nc = n - 1;

        // To the minus neighbour we must send the first g faces *after*
        // the shared interface face (local index 0); i = 1..g.
        if (o < 0)
            return 1;

        // To the plus neighbour we must send the last g faces *before*
        // the shared interface face at local index nc; i = nc-g..nc-1.
        // (Assuming g <= nc; which holds for your runs.)
        return nc - g;
    };

    // Destination: always write into canonical ghost locations; no MAC shift on
    // the receiving side. The indexing convention is:
    //   interior: [0 .. n-1]
    //   ghosts- : [-g .. -1]
    //   ghosts+ : [ n .. n+g-1]
    auto dst_start_axis = [&](int n, int g, int o)->int {
        if (o == 0) return 0;
        return (o < 0 ? -g : n);
    };

    auto src_box = [&](int ox,int oy,int oz)->Box {
        const int sx = (ox==0 ? nxI : ng);
        const int sy = (oy==0 ? nyI : ng);
        const int sz = (oz==0 ? nzI : ng);
        const int i0 = (ox==0 ? 0 : src_start_axis(nxI, ng, ox, faceI));
        const int j0 = (oy==0 ? 0 : src_start_axis(nyI, ng, oy, faceJ));
        const int k0 = (oz==0 ? 0 : src_start_axis(nzI, ng, oz, faceK));
    #ifndef NDEBUGEX
        auto ok = [&](int n, int s0, int s)->bool { return s0 >= 0 && s0 + s <= n; };
        assert(ok(nxI, i0, sx) && "src i out of range");
        assert(ok(nyI, j0, sy) && "src j out of range");
        assert(ok(nzI, k0, sz) && "src k out of range");
    #endif
        return { i0, j0, k0, sx, sy, sz };
    };

    auto dst_box = [&](int ox,int oy,int oz)->Box {
        const int sx = (ox==0 ? nxI : ng);
        const int sy = (oy==0 ? nyI : ng);
        const int sz = (oz==0 ? nzI : ng);
        const int i0d = (ox==0 ? 0 : dst_start_axis(nxI, ng, ox));
        const int j0d = (oy==0 ? 0 : dst_start_axis(nyI, ng, oy));
        const int k0d = (oz==0 ? 0 : dst_start_axis(nzI, ng, oz));
    #ifndef NDEBUGEX
        auto ghost_ok = [&](int n, int g, int o, int s0, int s)->bool {
            if (o == 0) return (s0 >= 0 && s0 + s <= n);
            if (o < 0)  return (s0 == -g && s == g);
            else        return (s0 ==  n && s == g);
        };
        assert(ghost_ok(nxI, ng, ox, i0d, sx) && "dst i must lie in ghosts");
        assert(ghost_ok(nyI, ng, oy, j0d, sy) && "dst j must lie in ghosts");
        assert(ghost_ok(nzI, ng, oz, k0d, sz) && "dst k must lie in ghosts");
    #endif
        return { i0d, j0d, k0d, sx, sy, sz };
    };

    struct Dir {
        int ox,oy,oz, peer, tag_recv, tag_send;
        Box s, d;
        std::size_t count;
    };
    std::array<Dir,26> dirs{};

    {   // build all 26 directions
        int p=0;
        for (int oz=-1; oz<=1; ++oz)
        for (int oy=-1; oy<=1; ++oy)
        for (int ox=-1; ox<=1; ++ox){
            if (ox==0 && oy==0 && oz==0) continue;
            const int peer = rank_of(ox,oy,oz);
            const int id   = dir_index(ox,oy,oz);
            auto opp       = opposite(ox,oy,oz);
            const int opp_id = dir_index(opp[0],opp[1],opp[2]);

            Dir d;
            d.ox=ox; d.oy=oy; d.oz=oz; d.peer=peer;
            d.s = src_box(ox,oy,oz);
            d.d = dst_box(ox,oy,oz);
            d.count = std::size_t(d.s.sx)*std::size_t(d.s.sy)*std::size_t(d.s.sz);
            // Irecv expects base+id; Isend uses peer’s expected tag = base+opp_id
            d.tag_recv = tag_base + id;
            d.tag_send = tag_base + opp_id;
            dirs[p++] = d;
        }
    }

    // helpers
    auto pack = [&](const Dir& d, std::vector<T>& buf){
        buf.resize(d.count);
        std::size_t q=0;
        for (int kk=0; kk<d.s.sz; ++kk){
            const int k = d.s.k0 + kk;
            for (int jj=0; jj<d.s.sy; ++jj){
                const int j = d.s.j0 + jj;
                for (int ii=0; ii<d.s.sx; ++ii){
                    const int i = d.s.i0 + ii;
                    buf[q++] = f(i,j,k);
                }
            }
        }
    };
    auto unpack = [&](const Dir& d, const std::vector<T>& buf){
        std::size_t q=0;
        for (int kk=0; kk<d.d.sz; ++kk){
            const int k = d.d.k0 + kk;
            for (int jj=0; jj<d.d.sy; ++jj){
                const int j = d.d.j0 + jj;
                for (int ii=0; ii<d.d.sx; ++ii){
                    const int i = d.d.i0 + ii;
                    f(i,j,k) = buf[q++];
                }
            }
        }
    };

    // PROC_NULL: clamp from this rank’s interior (zero-gradient)
    auto clamp_fill = [&](const Dir& d){
        std::vector<T> tmp; tmp.resize(d.count);
        std::size_t q=0;
        for (int kk=0; kk<d.d.sz; ++kk){
            int k = std::clamp(d.d.k0 + kk, 0, std::max(1,nzI)-1);
            for (int jj=0; jj<d.d.sy; ++jj){
                int j = std::clamp(d.d.j0 + jj, 0, std::max(1,nyI)-1);
                for (int ii=0; ii<d.d.sx; ++ii){
                    int i = std::clamp(d.d.i0 + ii, 0, std::max(1,nxI)-1);
                    tmp[q++] = f(i,j,k);
                }
            }
        }
        unpack(d, tmp);
    };

    // post all receives
    std::array<std::vector<T>,26> rbuf{}, sbuf{};
    std::vector<MPI_Request> reqs; reqs.reserve(52);
    for (int idx=0; idx<26; ++idx){
        const auto& d = dirs[idx];
        if (d.peer!=MPI_PROC_NULL && d.peer!=myrank){
            rbuf[idx].resize(d.count);
            MPI_Request r{};
            MPI_Irecv(rbuf[idx].data(),
                      static_cast<int>(d.count), Tmpi,
                      d.peer, d.tag_recv, comm, &r);
            reqs.push_back(r);
        }
    }

    // pack & send (or local paths)
    for (int idx=0; idx<26; ++idx){
        const auto& d = dirs[idx];
        if (d.peer == MPI_PROC_NULL){
            clamp_fill(d);
        } else if (d.peer == myrank){
            std::vector<T> tmp; pack(d,tmp); unpack(d,tmp);
        } else {
            pack(d, sbuf[idx]);
            MPI_Request s{};
            MPI_Isend(sbuf[idx].data(),
                      static_cast<int>(d.count), Tmpi,
                      d.peer, d.tag_send, comm, &s);
            reqs.push_back(s);
        }
    }

    if (!reqs.empty())
        MPI_Waitall(static_cast<int>(reqs.size()),
                    reqs.data(), MPI_STATUSES_IGNORE);

    // unpack remote receives
    for (int idx=0; idx<26; ++idx){
        const auto& d = dirs[idx];
        auto &rb = rbuf[idx];
        if (!rb.empty()) unpack(d, rb);
    }
    return true;
}

} // namespace core::mesh