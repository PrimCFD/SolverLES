#pragma once
#include "master/Views.hpp"
#include "master/plugin/Program.hpp"
#include "mesh/Mesh.hpp"
#include <memory>
#include <string>

/**
 * @file Scheduler.hpp
 * @brief Time loop and phase runner.
 *
 * @details
 * The scheduler advances time, overlaps halo exchanges with interior work, applies boundary
 * conditions, and triggers writer calls at a configured cadence. It executes the
 * :cpp:struct:`core::master::plugin::StepPlan` provided by the selected program.
 *
 * @rst
 * .. graphviz::
 *
 *    digraph S {
 *      rankdir=LR;
 *      node [shape=box];
 *      halos1 [label="halos.start()"];
 *      interior [label="actions(Interior)"];
 *      halos2 [label="halos.finish()"];
 *      bcs [label="bcs.apply()"];
 *      post [label="actions/globals(PostExchange, PostBC, EndStep)"];
 *      io [label="writer.write()"];
 *      halos1 -> interior -> halos2 -> bcs -> post -> io;
 *    }
 * @endrst
 *
 */

namespace core::master
{

struct RunContext;
class FieldCatalog;
namespace io
{
    class IWriter;
}

struct TimeControls
{
    double dt{0.0};
    double t_end{0.0};
    int write_every{0};
    std::string case_name{"case"};
};

class Scheduler
{
  public:
    Scheduler(const RunContext& rc, FieldCatalog& fields, io::IWriter& writer,
              const core::mesh::Mesh& mesh);

    void set_program(std::unique_ptr<plugin::IProgram> program);
    void run(const TimeControls& tc);

  private:
    const RunContext& rc_;
    FieldCatalog& fields_;
    io::IWriter& writer_;
    const core::mesh::Mesh& mesh_;
    std::unique_ptr<plugin::IProgram> program_;
};

} // namespace core::master
