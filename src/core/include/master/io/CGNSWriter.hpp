#pragma once
#include "master/io/IWriter.hpp"
#include "master/io/StagingPool.hpp"
#include "master/io/WritePlan.hpp"
#include "master/io/WriterConfig.hpp"
#include <memory>
#include <string>

/**
 * @file CGNSWriter.hpp
 * @brief Structured CGNS writer (guarded by SOLVERLES_WITH_CGNS).
 *
 * @details
 * Creates a structured zone once, then writes per-step FlowSolutions. Uses a
 * :cpp:struct:`WritePlan` and :cpp:class:`StagingPool` for packing and optional precision
 * conversion.
 *
 * @note Includes CGNS headers only when `SOLVERLES_WITH_CGNS` is defined.
 */

namespace core::master::io
{

class CGNSWriter : public IWriter
{
  public:
    explicit CGNSWriter(WriterConfig cfg);
    ~CGNSWriter() override;

    void open_case(const std::string& case_name) override;
    void write(const WriteRequest& req) override;
    void close() override;

  private:
    WriterConfig cfg_;
    WritePlan plan_;
    StagingPool pool_;

    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace core::master::io
