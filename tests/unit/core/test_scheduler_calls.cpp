#include "master/Master.hpp"
#include "master/io/IWriter.hpp"
#include "mesh/Mesh.hpp"
#include <catch2/catch_test_macros.hpp>
#include <vector>

using namespace core;

struct CountingWriter : master::io::IWriter
{
    int opens = 0, writes = 0, closes = 0;
    void open_case(const std::string&) override { ++opens; }
    void write(const master::io::WriteRequest&) override { ++writes; }
    void close() override { ++closes; }
};

static void setup_one_field(master::Master& m)
{
    std::vector<double> a(8, 0.0);
    m.fields().register_scalar("a", a.data(), sizeof(double), {2, 2, 2}, {1, 2, 4});
    m.fields().select_for_output("a");
}

static void run_case(int steps, int write_every, int expected_writes)
{
    master::RunContext rc{};

    mesh::Mesh mesh;
    mesh.local = {2, 2, 2};
    mesh.ng = 0;

    master::Master m(rc, mesh);
    setup_one_field(m);

    auto writer = std::make_unique<CountingWriter>();
    auto* W = writer.get();
    m.set_writer(std::move(writer));

    m.configure_program("noop", {}); // builtin program

    master::TimeControls tc;
    tc.dt = 1.0;
    tc.t_end = steps; // steps = round(t_end/dt)
    tc.write_every = write_every;

    m.run(tc);

    REQUIRE(W->opens == 1);
    REQUIRE(W->writes == expected_writes);
    REQUIRE(W->closes == 1);
}

TEST_CASE("Scheduler I/O cadence edge cases", "[scheduler][io]")
{
    SECTION("write_every=2 over 3 steps -> writes at steps 0 and 2")
    {
        run_case(/*steps=*/3, /*write_every=*/2, /*expected_writes=*/2);
    }
    SECTION("write_every=1 -> write every step (3)")
    {
        run_case(/*steps=*/3, /*write_every=*/1, /*expected_writes=*/3);
    }
    SECTION("write_every > nSteps -> only step 0 write")
    {
        run_case(/*steps=*/3, /*write_every=*/100, /*expected_writes=*/1);
    }
}
