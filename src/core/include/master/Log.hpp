#pragma once
#include <cstdio>
#include <cstdarg>
#include <string>
#include <atomic>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>

namespace core::master::logx {

enum class Level : int { Quiet=0, Error=1, Warn=2, Info=3, Debug=4 };

struct Config {
  Level level = Level::Info;     // default
  bool  color = true;            // ANSI colors if TTY
  bool  rank0_only = false;      // gate INFO/DEBUG to rank0
};

inline int  g_rank = 0;
inline bool g_is_tty = false;
inline std::atomic<Level> g_level{Level::Info};
inline std::atomic<bool>  g_color{true};
inline std::atomic<bool>  g_rank0_only{false};

inline bool is_tty(FILE* f) {
  return ::isatty(::fileno(f)) == 1;
}

inline Level level_from_env() {
  const char* v = std::getenv("SOLVER_LOG");
  if (!v) return Level::Info;
  std::string s(v);
  for (auto& c: s) c = tolower(c);
  if (s=="quiet") return Level::Quiet;
  if (s=="error") return Level::Error;
  if (s=="warn" || s=="warning") return Level::Warn;
  if (s=="info" || s=="mini") return Level::Info;
  if (s=="debug" || s=="full") return Level::Debug;
  return Level::Info;
}

inline void init(const Config& cfg = {}) {
  int inited=0; MPI_Initialized(&inited);
  if (inited) MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  g_is_tty = is_tty(stdout);
  g_level.store(cfg.level == Level::Info ? level_from_env() : cfg.level);
  g_color.store(cfg.color && g_is_tty);
  g_rank0_only.store(cfg.rank0_only);
}

inline const char* c(Level L) {
  static const char* K = "\033[0m";
  static const char* R = "\033[31m";
  static const char* Y = "\033[33m";
  static const char* C = "\033[36m";
  static const char* G = ""; // default
  if (!g_color.load()) return "";
 switch (L) {
    case Level::Error: return R;
    case Level::Warn:  return Y;
    case Level::Info:  return C;
    default:           return G;
  }
}
inline const char* reset() { return g_color.load() ? "\033[0m" : ""; }

inline bool gate(Level L) {
  if (L > g_level.load()) return true;                         // filtered by level
  if (g_rank0_only.load() && g_rank != 0 && L >= Level::Info)  // collapse info/debug to rank0
    return true;
  return false;
}

inline void vprint(Level L, const char* fmt, va_list ap) {
  if (gate(L)) return;
  if (L == Level::Info || L == Level::Debug)
    std::fprintf(stderr, "%s", c(L));
  if (L == Level::Error) std::fprintf(stderr, "%s[error]%s ", c(L), reset());
  else if (L == Level::Warn) std::fprintf(stderr, "%s[warn ]%s ", c(L), reset());
  // Compact rank tag for non-rank0 or for debug
  if (g_rank != 0 && L >= Level::Info)
    std::fprintf(stderr, "[r%d] ", g_rank);
  std::vfprintf(stderr, fmt, ap);
  if (L == Level::Info || L == Level::Warn || L == Level::Error)
    std::fprintf(stderr, "%s", reset());
  std::fflush(stderr);
}
inline void print(Level L, const char* fmt, ...) {
  va_list ap; va_start(ap, fmt); vprint(L, fmt, ap); va_end(ap);
}

// Convenience
#define LOGD(...) ::core::master::logx::print(::core::master::logx::Level::Debug, __VA_ARGS__)
#define LOGI(...) ::core::master::logx::print(::core::master::logx::Level::Info,  __VA_ARGS__)
#define LOGW(...) ::core::master::logx::print(::core::master::logx::Level::Warn,  __VA_ARGS__)
#define LOGE(...) ::core::master::logx::print(::core::master::logx::Level::Error, __VA_ARGS__)

} // namespace core::master::logx
