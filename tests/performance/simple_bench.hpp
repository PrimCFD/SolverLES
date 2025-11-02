#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <vector>

namespace bench
{

// Run `func()` `nIters` times, return mean & std-dev in microseconds
template <typename F> std::pair<double, double> run(const F& func, int nIters = 30)
{
    using clock = std::chrono::high_resolution_clock;
    std::vector<double> samples;
    samples.reserve(nIters);

    for (int i = 0; i < nIters; ++i)
    {
        auto t0 = clock::now();
        func();
        auto t1 = clock::now();
        samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }

    double sum = 0.0;
    for (double s : samples)
        sum += s;
    double mean = sum / samples.size();

    double var = 0.0;
    for (double s : samples)
        var += (s - mean) * (s - mean);
    double stddev = std::sqrt(var / (samples.size() - 1));

    return {mean, stddev}; // µs
}

// Like run(), but:
//  * Barrier before/after each iteration
//  * Measures wall-time per rank, then reduces to the global MAX (critical path)
//  * Returns mean/stddev of these global max samples (µs)
template <typename F>
std::pair<double, double> run_mpi_max(MPI_Comm comm, const F& func, int nIters = 30)
{
    using clock = std::chrono::high_resolution_clock;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    std::vector<double> samples;
    samples.reserve(nIters);

    for (int i = 0; i < nIters; ++i)
    {
        MPI_Barrier(comm);
        auto t0 = clock::now();
        func();
        MPI_Barrier(comm);
        auto t1 = clock::now();
        double local_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        double global_us = 0.0;
        MPI_Allreduce(&local_us, &global_us, 1, MPI_DOUBLE, MPI_MAX, comm);
        samples.push_back(global_us);
    }
    double sum = 0.0;
    for (double s : samples)
        sum += s;
    const double mean = sum / samples.size();
    double var = 0.0;
    for (double s : samples)
        var += (s - mean) * (s - mean);
    const double stddev = std::sqrt(var / (samples.size() - 1));
    return {mean, stddev};
}

// Print result in a CTest-friendly key=value format
inline void report(const std::string& name, double mean_us, double stddev_us, double bytes = 0)
{
    std::cout << name << " mean_us=" << mean_us << " stddev_us=" << stddev_us;
    if (bytes > 0)
        std::cout << " bytes=" << bytes << " MBps=" << (bytes / 1e6) / (mean_us * 1e-6);
    std::cout << '\n';
}

// Only rank 0 prints (to avoid N duplicate lines in mpirun output)
inline void report_root(MPI_Comm comm, const std::string& name, double mean_us, double stddev_us,
                        double bytes = 0)
{
    int r = 0;
    MPI_Comm_rank(comm, &r);
    if (r == 0)
        report(name, mean_us, stddev_us, bytes);
}

} // namespace bench
