#include "test_user_customization_point.h"

#include <geoma/geoma.h>
#include <gtest/gtest.h>


using namespace geoma;
using namespace geoma_test;


/**
 * Two co-linear segments: {{0.0, 0.0}, {1.0, 0.0}} and {{1.0, 0.0}, {3.0,
 * 0.0}}. No overlap. Canonicalization must produce a set with only one segment:
 * {{0.0, 0.0}, {3.0, 0.0}}
 */
TEST(TestAlgorithms, canonicalize_01)
{
  std::vector<Segment> segments{{{0.0, 0.0}, {1.0, 0.0}},
                                {{1.0, 0.0}, {3.0, 0.0}}};

  auto const result = canonicalize(segments, ToleranceConfig{});

  ASSERT_EQ(result.size(), 1);

  auto const& segment = *result.begin();

  EXPECT_DOUBLE_EQ(segment.p1.x, 0.0);
  EXPECT_DOUBLE_EQ(segment.p1.y, 0.0);
  EXPECT_DOUBLE_EQ(segment.p2.x, 3.0);
  EXPECT_DOUBLE_EQ(segment.p2.y, 0.0);
}

/**
 * Two non-co-linear segments: {{0.0, 0.0}, {3.0, 0.0}} and {{3.0, 0.0},
 * {3.0, 2.0}}. Not co-linear. Canonicalization must produce a set with two
 * segments.
 */
TEST(TestAlgorithms, canonicalize_02)
{
  std::vector<Segment> segments{{{0.0, 0.0}, {3.0, 0.0}},
                                {{3.0, 0.0}, {3.0, 2.0}}};

  auto const result = canonicalize(segments, ToleranceConfig{});

  EXPECT_EQ(result.size(), 2);
}

/**
 * Two co-linear segments: {{0.0, 0.0}, {2.0, 0.0}} and {{1.0, 0.0}, {3.0, 0.0}}
 * Two segments overlap.
 * Canonicalize must produce a set with only one segment:  {{0.0, 0.0}, {3.0,
 * 0.0}}
 */
TEST(TestAlgorithms, canonicalize_03)
{
  std::vector<Segment> segments{{{0.0, 0.0}, {2.0, 0.0}},
                                {{1.0, 0.0}, {3.0, 0.0}}};

  auto const result = canonicalize(segments, ToleranceConfig{});

  ASSERT_EQ(result.size(), 1);

  auto const& segment = *result.begin();

  EXPECT_DOUBLE_EQ(segment.p1.x, 0.0);
  EXPECT_DOUBLE_EQ(segment.p2.x, 3.0);
}

/**
 * Generate a total of num_input_segments such that after canonicalizing the
 * result has num_non_colinear_segments.
 *
 * Settings:
 *  num_input_segments         = 10M
 *  num_non_colinear_segments  = 1k
 *
 */
TEST(TestAlgorithms, canonicalize_04)
{
  constexpr std::size_t num_input_segments        = 10000000;
  constexpr std::size_t num_non_colinear_segments = 1000;

  ASSERT_GT(num_non_colinear_segments, 0);
  ASSERT_GE(num_input_segments, num_non_colinear_segments);

  std::vector<Segment> segments;

  for (std::size_t i = 0; i < num_input_segments; ++i)
    {
      std::size_t const line  = i % num_non_colinear_segments;
      std::size_t const index = i / num_non_colinear_segments;

      double const x = static_cast<double>(index);
      double const y = static_cast<double>(line);

      segments.emplace_back(Segment{{x, y}, {x + 1.0, y}});
    }

  ToleranceConfig config{};
  config.tolerance = 1e-8;

  auto const start  = std::chrono::steady_clock::now();
  auto const result = canonicalize(segments, ToleranceConfig{});
  auto const elapsed =
    std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();

  std::cout << "canonicalize: " << elapsed << " s\n";
  EXPECT_EQ(segments.size(), num_input_segments);
  EXPECT_EQ(result.size(), num_non_colinear_segments);
}