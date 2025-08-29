#pragma once
#include "master/Views.hpp"
#include <string>
#include <vector>

/**
 * @file IWriter.hpp
 * @brief Output abstraction used by the scheduler.
 *
 * @details
 * Writers receive a :cpp:struct:`core::master::io::WriteRequest` that contains
 * a set of host-side views and the step/time metadata. The base interface is
 * synchronous; an asynchronous wrapper can be layered later without changing
 * call sites.
 *
 * @rst
 * .. note::
 *    Writers should not mutate field data. Stage/pack into internal buffers
 *    if a format requires reordering.
 * @endrst
 */

namespace core
{
namespace mesh
{
    class Layout;
}

namespace master::io
{

    struct WriteRequest
    {
        int step{0};
        double time{0.0};
        std::vector<AnyFieldView> fields; // host-side views (stage to device if needed elsewhere)
    };

    class IWriter
    {
      public:
        virtual ~IWriter() = default;

        virtual void open_case(const std::string& root, const core::mesh::Layout& layout) = 0;
        virtual void write(const WriteRequest& req) = 0; // blocking in the base interface
        virtual void close() = 0;
    };

} // namespace master::io
} // namespace core
