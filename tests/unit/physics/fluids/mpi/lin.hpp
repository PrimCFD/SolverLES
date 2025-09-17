#pragma once
#include <cstddef>
// Row-major linear index, matching Fortran kernels’ layout
static inline std::size_t lin(int I, int J, int K, int nx_tot, int ny_tot)
{
    return static_cast<std::size_t>(I + nx_tot * (J + ny_tot * K));
}