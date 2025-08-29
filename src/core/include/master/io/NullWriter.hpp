#pragma once
#include "master/io/IWriter.hpp"

namespace core::master::io {

/// No-op writer (useful for perf tests and smoke tests)
class NullWriter final : public IWriter {
public:
  void open_case(const std::string&, const core::mesh::Layout&) override {}
  void write(const WriteRequest&) override {}
  void close() override {}
};

} // namespace core::master::io
