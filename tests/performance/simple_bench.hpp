#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
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

// Print result in a CTest-friendly key=value format
inline void report(const std::string& name, double mean_us, double stddev_us, double bytes = 0)
{
    std::cout << name << " mean_us=" << mean_us << " stddev_us=" << stddev_us;
    if (bytes > 0)
        std::cout << " bytes=" << bytes << " MBps=" << (bytes / 1e6) / (mean_us * 1e-6);
    std::cout << '\n';
}
} // namespace bench
