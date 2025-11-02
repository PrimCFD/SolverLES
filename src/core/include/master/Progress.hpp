#pragma once
#include <cstdio>
#include <chrono>
#include <string>
#include <algorithm>
#include <unistd.h>

namespace core::master::prog {
struct Bar {
  int total_steps = 1;
  std::chrono::steady_clock::time_point t0{};
  bool is_tty = false;
  int  last_drawn = -1;   // last integer percent drawn

  void start(int steps) {
    total_steps = std::max(1, steps);
    t0 = std::chrono::steady_clock::now();
    is_tty = ::isatty(::fileno(stderr));
    last_drawn = -1;
  }
  void update(int step) {
    if (!is_tty) return;               // no redraw spam under ctest/redirection
    int pct = int( (100.0 * step) / std::max(1, total_steps) );
    if (pct == last_drawn) return;     // avoid overhead
    last_drawn = pct;
    const int width = 40;
    int fill = (pct * width) / 100;
    auto now = std::chrono::steady_clock::now();
    double dt  = std::chrono::duration<double>(now - t0).count();
    double eta = (pct>0) ? dt*(100.0/pct - 1.0) : 0.0;
    std::fprintf(stderr, "\r[");
    for (int i=0;i<width;i++) std::fputc(i<fill ? '=' : ' ', stderr);
    std::fprintf(stderr, "] %3d%%  t=%.1fs  ETA=%.1fs", pct, dt, eta);
    std::fflush(stderr);
  }
  void finish() {
    if (!is_tty) return;
    update(total_steps);
    std::fprintf(stderr, "\n");
  }
};
} // namespace core::master::prog
