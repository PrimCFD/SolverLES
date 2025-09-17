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

TEST_CASE("Scheduler runs builtin noop program and triggers I/O cadence", "[scheduler][noop]")
{
    master::RunContext rc{};

    // New API: Master requires a Mesh
    mesh::Mesh mesh;
    mesh.local = {2, 2, 2};
    mesh.ng = 0; // no ghosts in this simple test

    master::Master m(rc, mesh);

    // one tiny field selected for output (totals == interior when ng==0)
    std::vector<double> a(8, 0.0);
    m.fields().register_scalar("a", a.data(), sizeof(double), {2, 2, 2}, {1, 2, 4});
    m.fields().select_for_output("a");

    auto writer = std::make_unique<CountingWriter>();
    auto* wptr = writer.get();
    m.set_writer(std::move(writer));

    // builtin program from PluginHost ctor
    m.configure_program("noop", {});

    master::TimeControls tc;
    tc.dt = 0.1;
    tc.t_end = 0.3;     // 3 steps
    tc.write_every = 2; // write on cadence multiples after stepping: step 2 → 1 write

    m.run(tc);

    REQUIRE(wptr->opens == 1);
    REQUIRE(wptr->writes == 1);
    REQUIRE(wptr->closes == 1);
}
