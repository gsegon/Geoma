
#include "geoma/geoma.h"

#include <gtest/gtest.h>
#include <numbers>
#include <vector>

using namespace geoma;

/**
 * Two sets in canonical form but - equal.
 */
TEST(TestAlgorithms, verify_equal_01) {
  std::vector<Segment> source{
    {{0.0, 0.0}, {3.0, 0.0}},
    {{3.0, 0.0}, {3.0, 2.0}}
  };

  Transform transform{
    .angle = std::numbers::pi / 2.0,
    .translation = {5.0, 3.0}
  };

  std::vector<Segment> target{
    {{5.0, 3.0}, {5.0, 6.0}},
    {{5.0, 6.0}, {3.0, 6.0}}
  };

  EXPECT_TRUE(verify_equal(source, target, transform, 1e-8));
}

/**
 * Two sets in canonical form - not equal.
 */
TEST(TestAlgorithms, verify_equal_02) {
  std::vector<Segment> source{
    {{0.0, 0.0}, {3.0, 0.0}},
    {{3.0, 0.0}, {3.0, 2.0}}
  };

  std::vector<Segment> target{
    {{0.0, 0.0}, {3.0, 0.0}},
    {{3.0, 0.0}, {3.0, 3.0}}
  };

  Transform transform{};

  EXPECT_FALSE(verify_equal(source, target, transform, 1e-8));
}


/**
 * Stress test: 100k segments.
 */
TEST(TestAlgorithms, verify_equal_03) {
  constexpr std::size_t num_segments = 100000;
  constexpr double tol = 1e-8;

  std::vector<Segment> source;
  std::vector<Segment> target;

  source.reserve(num_segments);
  target.reserve(num_segments);

  Transform transform{
    .angle = std::numbers::pi / 2.0,
    .translation = {100.0, 200.0}
  };

  // Generate distinct, already-canonicalized segments.
  for (std::size_t i = 0; i < num_segments; ++i) {
    const double y = static_cast<double>(i);

    Segment segment{
      {0.0, y},
      {10.0, y}
    };

    source.push_back(segment);
    target.push_back(transform.apply(segment));
  }

  // Verify that input ordering does not matter.
  std::ranges::reverse(target);

  const auto start = std::chrono::steady_clock::now();

  const bool equal = verify_equal(source, target, transform, tol);

  const auto elapsed = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start
  ).count();

  std::cout << "verify_equal (" << num_segments
            << " segments): " << elapsed << " s\n";

  EXPECT_TRUE(equal);

  // Modify one segment to make the geometry different.
  target.front().p2.y += 1.0;

  EXPECT_FALSE(verify_equal(source, target, transform, tol));
}