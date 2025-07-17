#!/bin/bash
# format-all.sh — Format C++, Fortran, and CMake files

set -e

echo "🔧 Formatting C++ with clang-format..."
find src/ tests/ cmake/ -type f \( -name "*.cpp" -o -name "*.hpp" -o -name "*.cc" -o -name "*.h" \) \
  -exec clang-format -i {} +

echo "🧮 Formatting Fortran with fprettify..."
find src/ tests/ -type f -name "*.f90" \
  -exec fprettify --config .fprettifyrc {} +

echo "🛠️  Formatting CMake with cmake-format..."
find . -type f \( -name "CMakeLists.txt" -o -name "*.cmake" \) \
  -exec cmake-format -i {} +

echo "✅ All files formatted!"
