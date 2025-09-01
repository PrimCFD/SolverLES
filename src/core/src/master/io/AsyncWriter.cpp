#include "master/io/AsyncWriter.hpp"
#include <chrono>

namespace core::master::io
{

AsyncWriter::AsyncWriter(std::unique_ptr<IWriter> sink) : AsyncWriter(std::move(sink), Options{}) {}

AsyncWriter::AsyncWriter(std::unique_ptr<IWriter> sink, Options o) : opt_(o), sink_(std::move(sink))
{
}

void AsyncWriter::open_case(const std::string& case_name)
{
    if (sink_)
        sink_->open_case(case_name);
    stop_ = false;
    worker_ = std::thread([this] { run_(); });
}

void AsyncWriter::write(const WriteRequest& req)
{
    if (!sink_ || stop_)
        return;
    {
        std::unique_lock<std::mutex> lk(mtx_);
        if (opt_.max_queue && q_.size() >= opt_.max_queue)
        {
            if (opt_.drop_to_null_on_overflow)
            {
                return; // drop to protect cadence
            }
            cv_.wait(lk, [&] { return q_.size() < opt_.max_queue || stop_; });
            if (stop_)
                return;
        }
        q_.push(Job{req});
    }
    cv_.notify_one();
}

void AsyncWriter::close()
{
    if (worker_.joinable())
    {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            stop_ = true;
        }
        cv_.notify_all();
        worker_.join();
    }
    if (sink_)
        sink_->close();
}

void AsyncWriter::run_()
{
    while (true)
    {
        Job job{};
        {
            std::unique_lock<std::mutex> lk(mtx_);
            cv_.wait(lk, [&] { return stop_ || !q_.empty(); });
            if (stop_ && q_.empty())
                break;
            job = q_.front();
            q_.pop();
        }
        sink_->write(job.req);
        cv_.notify_all(); // wake producer waiting on queue space
    }
}

} // namespace core::master::io
