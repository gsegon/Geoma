
#include <geoma/geoma.h>

#include <gtest/gtest.h>
#include <numbers>
#include <vector>

using namespace geoma;

TEST(TestAlgorithms, find_transform_01) {
  std::vector<Segment> source{
    {{0.0, 0.0}, {3.0, 0.0}},
    {{3.0, 0.0}, {3.0, 2.0}}
  };

  Transform expected{
    .angle = std::numbers::pi / 2.0,
    .translation = {5.0, 3.0}
  };

  std::vector<Segment> target;
  target.reserve(source.size());

  for (auto const& segment : source) {
    target.push_back(expected.apply(segment));
  }

  ToleranceConfig config{};
  config.tol = 1e-8;

  auto result = find_transform(source, target, config);

  ASSERT_TRUE(result.has_value());
  EXPECT_TRUE(verify_equal(source, target, *result, config.tol));
}

TEST(TestAlgorithms, find_transform_02) {
  std::vector<Segment> source{
    {{0.0, 0.0}, {3.0, 0.0}},
    {{3.0, 0.0}, {3.0, 2.0}}
  };

  std::vector<Segment> target{
    {{0.0, 0.0}, {3.0, 0.0}},
    {{10.0, 0.0}, {10.0, 2.0}}
  };

  ToleranceConfig config{};
  config.tol = 1e-8;

  auto result = find_transform(source, target, config);

  EXPECT_FALSE(result.has_value());
}


TEST(TestAlgorithms, find_transform_03) {
  constexpr std::size_t num_segments = 100000;
  constexpr double tol = 1e-8;

  std::vector<Segment> source;
  std::vector<Segment> target;

  source.reserve(num_segments);
  target.reserve(num_segments);

  Transform expected{
    .angle = std::numbers::pi / 2.0,
    .translation = {100.0, 200.0}
  };

  // Generate distinct, already-canonicalized segments.
  for (std::size_t i = 0; i < num_segments; ++i) {
    const double x = static_cast<double>(i) * 7.0;
    const double y = static_cast<double>(i) * 5.0;
    const double length = static_cast<double>(i + 1);

    Segment segment{
      {x, y},
      {x + length, y}
    };

    source.push_back(segment);
    target.push_back(expected.apply(segment));
  }

  // Ensure matching does not depend on input order.
  std::ranges::reverse(target);

  ToleranceConfig config{};
  config.tol = tol;

  const auto start = std::chrono::steady_clock::now();

  const auto result = find_transform(source, target, config);

  const auto elapsed = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start
  ).count();

  std::cout << "find_transform (" << num_segments
            << " segments): " << elapsed << " s\n";

  ASSERT_TRUE(result.has_value());

  EXPECT_TRUE(verify_equal(source, target, *result, tol));
}