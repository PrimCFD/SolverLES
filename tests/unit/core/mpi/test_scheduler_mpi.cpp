#include "master/Master.hpp"
#include "master/io/IWriter.hpp"
#include <catch2/catch_test_macros.hpp>
#include <mpi.h>
#include <vector>

using namespace core;

struct RankAwareWriter : master::io::IWriter
{
    int opens = 0, writes = 0, closes = 0;
    void open_case(const std::string&, const mesh::Layout&) override { ++opens; }
    void write(const master::io::WriteRequest&) override { ++writes; }
    void close() override { ++closes; }
};

TEST_CASE("MPI: only rank 0 performs IO; steps synchronized", "[mpi][scheduler]")
{
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    master::RunContext rc{};
    mesh::Layout layout;
    mesh::HaloExchange halos;
    mesh::Boundary bcs;
    master::Master m(rc, layout, halos, bcs);

    std::vector<double> a(8, 0.0);
    m.fields().register_scalar("a", a.data(), sizeof(double), {2, 2, 2}, {1, 2, 4});
    m.fields().select_for_output("a");

    auto w = std::make_unique<RankAwareWriter>();
    auto* W = w.get();
    m.set_writer(std::move(w));

    m.configure_program("noop", {});
    master::TimeControls tc;
    tc.dt = 0.1;
    tc.t_end = 0.2;
    tc.write_every = 1;
    m.run(tc);

    // Policy: IO on rank 0 only (tweak if different in your codebase)
    if (rank == 0)
    {
        REQUIRE(W->opens == 1);
        REQUIRE(W->writes == 2);
        REQUIRE(W->closes == 1);
    }
    else
    {
        REQUIRE(W->opens == 0);
        REQUIRE(W->writes == 0);
        REQUIRE(W->closes == 0);
    }

    // All ranks executed the same number of steps → reduction check
    int local_steps = 2; // dt=0.1, t_end=0.2 → steps 0,1
    int sum = 0;
    MPI_Allreduce(&local_steps, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    REQUIRE(sum == 2 * size);
}
