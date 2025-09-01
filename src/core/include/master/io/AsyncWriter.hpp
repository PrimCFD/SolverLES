#pragma once
#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "master/io/IWriter.hpp" // WriteRequest + IWriter (void API)

/**
 * @file AsyncWriter.hpp
 * @brief Threaded writer decorator with bounded queue and backpressure.
 *
 * @details
 * Wraps any :cpp:class:`IWriter`. `open_case()` spawns a worker thread; `write()` enqueues jobs;
 * `close()` joins the thread and flushes remaining work. The queue capacity and overflow policy
 * are configurable via :cpp:struct:`Options`.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   auto sink = std::make_unique<XdmfHdf5Writer>(cfg);
 *   AsyncWriter W(std::move(sink), {.max_queue=4, .drop_to_null_on_overflow=true});
 *   W.open_case("runA");
 *   W.write(req); // returns quickly
 *   W.close();    // flush & join
 * @endrst
 */

namespace core::master::io
{

class AsyncWriter : public IWriter
{
  public:
    /// \cond DOXYGEN_EXCLUDE
    struct Options
    {
        std::size_t max_queue = 2; // 0 = unbounded (not recommended)
        bool drop_to_null_on_overflow = false;
    };

    // No default arg here; provide a delegating ctor instead.
    explicit AsyncWriter(std::unique_ptr<IWriter> sink);
    AsyncWriter(std::unique_ptr<IWriter> sink, Options o);

    /// \endcond

    ~AsyncWriter() override { close(); }

    void open_case(const std::string& case_name) override;
    void write(const WriteRequest& req) override;
    void close() override;

  private:
    struct Job
    {
        WriteRequest req;
    };

    Options opt_;
    std::unique_ptr<IWriter> sink_;
    std::thread worker_;
    std::mutex mtx_;
    std::condition_variable cv_;
    std::queue<Job> q_;
    std::atomic<bool> stop_{false};

    void run_();
};

} // namespace core::master::io
