#pragma once
// Minimal seam types to satisfy Master/Scheduler construction & calls
namespace core { namespace mesh {
  class Layout { public: Layout() = default; };
  class HaloExchange { public: void start(); void finish(); };
  class Boundary     { public: void apply(); };
}} // namespace core::mesh
