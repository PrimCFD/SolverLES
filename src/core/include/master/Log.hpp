#pragma once
#include <atomic>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <string>

namespace core::master::logx
{

enum class Level : int
{
    Quiet = 0,
    Error = 1,
    Warn = 2,
    Info = 3,
    Debug = 4
};

struct Config
{
    Level level = Level::Info; // default verbosity (overridden by SOLVER_LOG if level==Info)
    bool rank0_only = false;   // gate INFO/DEBUG to rank 0
};

inline int g_rank = 0;
inline std::atomic<Level> g_level{Level::Info};
inline std::atomic<bool> g_rank0_only{false};

inline Level level_from_env()
{
    const char* v = std::getenv("SOLVER_LOG");
    if (!v)
        return Level::Info;
    std::string s(v);
    for (auto& c : s)
        c = static_cast<char>(::tolower(c));
    if (s == "quiet")
        return Level::Quiet;
    if (s == "error")
        return Level::Error;
    if (s == "warn" || s == "warning")
        return Level::Warn;
    if (s == "info" || s == "mini")
        return Level::Info;
    if (s == "debug" || s == "full")
        return Level::Debug;
    return Level::Info;
}

inline void init(const Config& cfg = {})
{
    int inited = 0;
    MPI_Initialized(&inited);
    if (inited)
        MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);

    // If caller leaves level at Info, allow SOLVER_LOG to override (preserves existing behavior)
    g_level.store(cfg.level == Level::Info ? level_from_env() : cfg.level);
    g_rank0_only.store(cfg.rank0_only);
}

inline const char* level_tag(Level L)
{
    switch (L)
    {
    case Level::Error:
        return "[error] ";
    case Level::Warn:
        return "[warn ] ";
    case Level::Info:
        return "[info ] ";
    case Level::Debug:
        return "[debug] ";
    default:
        return "";
    }
}

inline bool gate(Level L)
{
    if (L > g_level.load())
        return true; // filtered by level
    if (g_rank0_only.load() && g_rank != 0 && L >= Level::Info)
        return true; // collapse info/debug to rank0
    return false;
}

inline void vprint(Level L, const char* fmt, va_list ap)
{
    if (gate(L))
        return;
    // Level tag
    std::fputs(level_tag(L), stderr);
    // Compact rank tag for non-rank0 on info/debug
    if (g_rank != 0 && L >= Level::Info)
        std::fprintf(stderr, "[r%d] ", g_rank);
    std::vfprintf(stderr, fmt, ap);
    std::fflush(stderr);
}

inline void print(Level L, const char* fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprint(L, fmt, ap);
    va_end(ap);
}

// Convenience
#define LOGD(...) ::core::master::logx::print(::core::master::logx::Level::Debug, __VA_ARGS__)
#define LOGI(...) ::core::master::logx::print(::core::master::logx::Level::Info, __VA_ARGS__)
#define LOGW(...) ::core::master::logx::print(::core::master::logx::Level::Warn, __VA_ARGS__)
#define LOGE(...) ::core::master::logx::print(::core::master::logx::Level::Error, __VA_ARGS__)

} // namespace core::master::logx
