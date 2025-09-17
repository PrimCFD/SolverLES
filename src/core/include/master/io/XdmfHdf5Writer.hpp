#pragma once
#include "master/io/IWriter.hpp"
#include "master/io/StagingPool.hpp"
#include "master/io/WritePlan.hpp"
#include "master/io/WriterConfig.hpp"
#include <memory>
#include <string>

/**
 * @file XdmfHdf5Writer.hpp
 * @brief HDF5 datasets per step + XDMF XML index.
 *
 * @details
 * On the first write, builds a :cpp:struct:`WritePlan` and a :cpp:class:`StagingPool`. Each
 * step writes datasets into the HDF5 file under `/Step_xxxxxx`, then rewrites an XDMF XML index
 * referencing those datasets. Precision down-cast is applied during packing if configured.
 *
 */

namespace core::master::io
{

class XdmfHdf5Writer : public IWriter
{
  public:
    explicit XdmfHdf5Writer(WriterConfig cfg);
    ~XdmfHdf5Writer() override;

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
