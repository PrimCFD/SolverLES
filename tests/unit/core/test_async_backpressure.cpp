#include "master/io/AsyncWriter.hpp"
#include "master/io/IWriter.hpp"
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <memory>

using namespace core::master::io;

struct SlowSink : IWriter
{
    int opened = 0, written = 0, closed = 0;
    void open_case(const std::string&) override { opened++; }
    void write(const WriteRequest&) override
    {
        written++;
        std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }
    void close() override { closed++; }
};

TEST_CASE("AsyncWriter enforces backpressure", "[io][async]")
{
    auto sink = std::make_unique<SlowSink>();
    AsyncWriter::Options opt;
    opt.max_queue = 1;
    opt.drop_to_null_on_overflow = false;
    AsyncWriter W(std::move(sink), opt);

    W.open_case("async_bp");

    WriteRequest r;
    r.step = 0;
    r.time = 0.0;

    auto t0 = std::chrono::steady_clock::now();
    W.write(r); // queued
    W.write(r); // must block until worker drains one
    auto t1 = std::chrono::steady_clock::now();

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    REQUIRE(ms >= 20);

    W.close();
}
