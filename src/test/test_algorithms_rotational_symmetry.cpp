
#include "geoma/geoma.h"

#include <gtest/gtest.h>
#include <numbers>
#include <vector>
#include <chrono>
#include <iostream>


using namespace geoma;

/**
 * Set has rotational symmetry.
 */
TEST(TestAlgorithms, has_rotational_symmetry_01) {
  std::vector<Segment> segments{
    {{0.0, 0.0}, {2.0, 0.0}},
    {{2.0, 0.0}, {2.0, 2.0}},
    {{2.0, 2.0}, {0.0, 2.0}},
    {{0.0, 2.0}, {0.0, 0.0}}
  };

  auto result = has_rotational_symmetry<Segment, Transform>(segments, 1e-8);
  EXPECT_TRUE(result);
}

/**
 * Set doesn't have rotational symmetry.
 */
TEST(TestAlgorithms, has_rotational_symmetry_02) {
  std::vector<Segment> segments{
    {{0.0, 0.0}, {4.0, 0.0}},
    {{4.0, 0.0}, {1.0, 2.0}},
    {{1.0, 2.0}, {0.0, 0.0}}
  };

  EXPECT_FALSE(has_rotational_symmetry(segments, 1e-8));
}

/**
 * Stress test: 100k segments.
 */
TEST(TestAlgorithms, has_rotational_symmetry_03) {
  constexpr std::size_t num_segments = 100000;
  constexpr double radius = 100.0;
  constexpr double tol = 1e-8;

  std::vector<Segment> segments;
  segments.reserve(num_segments);

  // Generate a regular polygon approximating a circle.
  for (std::size_t i = 0; i < num_segments; ++i) {
    const double a1 =
      2.0 * std::numbers::pi * static_cast<double>(i) / num_segments;

    const double a2 =
      2.0 * std::numbers::pi * static_cast<double>(i + 1) / num_segments;

    segments.push_back(Segment{
      {radius * std::cos(a1), radius * std::sin(a1)},
      {radius * std::cos(a2), radius * std::sin(a2)}
    });
  }

  ASSERT_EQ(segments.size(), num_segments);

  const auto start = std::chrono::steady_clock::now();

  const bool result = has_rotational_symmetry(segments, tol);

  const auto elapsed = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start
  ).count();

  std::cout << "has_rotational_symmetry ("
            << num_segments << " segments): "
            << elapsed << " s\n";

  EXPECT_TRUE(result);
}