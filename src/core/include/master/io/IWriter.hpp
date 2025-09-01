#pragma once
#include "master/Views.hpp"
#include <string>
#include <vector>

/**
 * @file IWriter.hpp
 * @brief Output abstraction used by the scheduler.
 *
 * @details
 * Writers receive a :cpp:struct:`WriteRequest` with **host-side** views and step/time metadata.
 * The base interface is synchronous; use :cpp:class:`AsyncWriter` as a decorator for background
 * I/O.
 *
 * @rst
 * .. code-block:: cpp
 *
 *   using core::master::io::IWriter;
 *   struct CountingWriter : IWriter {
 *     int writes = 0;
 *     void open_case(const std::string&) override {}
 *     void write(const WriteRequest&) override { ++writes; }
 *     void close() override {}
 *   };
 * @endrst
 *
 * @note `open_case()` returns `void` (legacy boolean return was removed).
 */

namespace core::master::io
{
struct WriteRequest
{
    int step{0};
    double time{0.0};
    std::vector<AnyFieldView> fields; // host-side views
};

class IWriter
{
  public:
    virtual ~IWriter() = default;

    // Removed deprecated Layout parameter.
    virtual void open_case(const std::string& case_name) = 0;
    virtual void write(const WriteRequest& req) = 0; // blocking in the base interface
    virtual void close() = 0;
};

} // namespace core::master::io
