#ifndef UNTITLED1_NORMALIZE_H
#define UNTITLED1_NORMALIZE_H
#include "config.h"
#include "primitives.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace UngaBunga
{
  // Public API - Start

  /**
   * Check for invalid segments, do any preconditioning that might be beneficial
   * for performance. Can be no-op initially.
   *
   * @param segments
   * @param tolerance_config
   * @return
   */
  [[nodiscard]]
  std::set<Segment>
  normalize(const std::set<Segment> &segments,
            const ToleranceConfig   &tolerance_config);

  /**
   * Convert a given set of segments into a set of segments where neighbouring
   * colinear segments are merged into one. Resulting set must represent the
   * same geometric union as the input set.
   * @param segments
   * @param tolerance_config
   * @return
   */
  [[nodiscard]]
  std::set<Segment>
  canonicalize(const std::set<Segment> &segments,
               const ToleranceConfig   &tolerance_config);

  /**
   * Generates candidates of the transform for two sets of linear segments in
   * canonical form. Assumption is that the two sets represent identical
   * geometric union under some transformation. Candidates are later checked
   * with `verify_equal`.
   *
   * @param source
   * @param target
   * @param tolerance_config
   * @return
   */
  [[nodiscard]]
  std::vector<Transform>
  generate_candidates(const std::set<Segment> &source,
                      const std::set<Segment> &target,
                      const ToleranceConfig   &tolerance_config);

  /**
   * Check if `target = transform(source)`.
   *
   * @param source
   * @param target
   * @param transform
   * @param tol
   * @return
   */
  [[nodiscard]]
  bool
  verify_equal(const std::set<Segment> &source,
               const std::set<Segment> &target,
               const Transform         &transform,
               double                  tol);


  /**
   * Find the transform
   * @param source
   * @param target
   * @param tolerance_config
   * @return
   */
  [[nodiscard]]
  std::optional<Transform>
  find_transform(
    std::set<Segment> const& source,
    std::set<Segment> const& target,
    ToleranceConfig const& tolerance_config
  );

  // Public API - End

  namespace implementation_vdev
  {
    inline std::set<Segment>
    normalize(const std::set<Segment> &segments,
              const ToleranceConfig   &tolerance_config);

    [[nodiscard]]
    inline std::set<Segment>
    canonicalize(const std::set<Segment> &segments,
                 const ToleranceConfig   &tolerance_config);

    [[nodiscard]]
    std::vector<Transform>
    generate_candidates(const std::set<Segment> &source,
                        const std::set<Segment> &target,
                        const ToleranceConfig   &tolerance_config);

    [[nodiscard]]
    bool
    verify_equal(const std::set<Segment> &source,
                 const std::set<Segment> &target,
                 const Transform         &transform,
                 double                  tol);

    [[nodiscard]]
    std::optional<Transform>
    find_transform(
      std::set<Segment> const& source,
      std::set<Segment> const& target,
      ToleranceConfig const& tolerance_config
    );

  } // namespace implementation_vdev

  [[nodiscard]]
  inline std::set<Segment>
  canonicalize(const std::set<Segment> &segments,
               const ToleranceConfig   &tolerance_config)
  {
    return implementation_vdev::canonicalize(segments, tolerance_config);
  }

  [[nodiscard]]
  inline std::vector<Transform>
  generate_candidates(const std::set<Segment> &source,
                      const std::set<Segment> &target,
                      const ToleranceConfig   &tolerance_config)
  {
    return implementation_vdev::generate_candidates(source,
                                                    target,
                                                    tolerance_config);
  }

  [[nodiscard]]
  inline bool
  verify_equal(const std::set<Segment> &source,
               const std::set<Segment> &target,
               const Transform         &transform,
               double                   tol)
    {
      return implementation_vdev::verify_equal(source, target, transform, tol);
    }

  [[nodiscard]]
    inline std::optional<Transform>
    find_transform(
      std::set<Segment> const& source,
      std::set<Segment> const& target,
      ToleranceConfig const& tolerance_config
    )
  {
    return implementation_vdev::find_transform(source, target, tolerance_config);
  }

  namespace implementation_vdev
  {

    namespace detail
    {
      // Signed perpendicular distances to the
      // supporting line of s1.
      inline auto signed_distance = [](const Segment &s,
                                       const Point   &p) -> double {
        const double dx = s.p2.x - s.p1.x;
        const double dy = s.p2.y - s.p1.y;

        const double len = std::hypot(dx, dy);

        if (len == 0.0)
          throw std::invalid_argument("signed_distance: degenerate segment");

        return (dx * (p.y - s.p1.y) - dy * (p.x - s.p1.x)) / len;
      };

      inline auto project = [](const Segment &s, const Point &p) -> double {
        const double dx = s.p2.x - s.p1.x;
        const double dy = s.p2.y - s.p1.y;

        const double len = std::hypot(dx, dy);

        if (len == 0.0)
          throw std::invalid_argument("project: degenerate segment");

        return (dx * (p.x - s.p1.x) + dy * (p.y - s.p1.y)) / len;
      };

      inline auto make_transform = [](const Point &source_p1,
                                      const Point &source_p2,
                                      const Point &target_p1,
                                      const Point &target_p2) -> Transform {
        const double sx = source_p2.x - source_p1.x;
        const double sy = source_p2.y - source_p1.y;

        const double tx = target_p2.x - target_p1.x;
        const double ty = target_p2.y - target_p1.y;

        const double angle = std::atan2(ty, tx) - std::atan2(sy, sx);

        const double c = std::cos(angle);
        const double s = std::sin(angle);

        Point translation{target_p1.x - (c * source_p1.x - s * source_p1.y),
                          target_p1.y - (s * source_p1.x + c * source_p1.y)};

        return Transform{.angle = angle, .translation = translation};
      };
    } // namespace detail

    // No-op for now
    inline std::set<Segment>
    normalize(const std::set<Segment> &segments,
              const ToleranceConfig   &tolerance_config)
    {
      std::set<Segment> result = segments;
      return result;
    }



    [[nodiscard]]
    inline std::set<Segment>
    canonicalize(const std::set<Segment> &segments,
                 const ToleranceConfig   &tolerance_config)
    {
      const double tol = tolerance_config.tol;

      auto distance = [](const auto &a, const auto &b) {
        return std::hypot(b.x - a.x, b.y - a.y);
      };

      auto length = [&](const Segment &s) { return distance(s.p1, s.p2); };

      // Check if two segments are collinear and their
      // projected intervals overlap or touch.
      auto try_merge = [](const Segment &s1,
                          const Segment &s2,
                          double         tol) -> std::optional<Segment> {
        const double dx = s1.p2.x - s1.p1.x;
        const double dy = s1.p2.y - s1.p1.y;

        const double len = std::hypot(dx, dy);

        // Reject segments with length smaller then tolerance.
        // Calling canonicalize is expected with normalized segment set
        // (deduplicated, degenerate removed or otherwise handled, ...)
        if (len <= tol)
          throw std::logic_error{
            "segment merge:"
            " The resulting segment has a length smaller, or equal, to tolerance"
            ". This is not allowed. "};

        // Unit direction along s1
        const double ux = dx / len;
        const double uy = dy / len;

        // Check whether both endpoints of s2 lie on the supporting
        // line of s1, i.e. check if the two given segments are co-linear.
        if (std::abs(detail::signed_distance(s1, s2.p1)) > tol ||
            std::abs(detail::signed_distance(s1, s2.p2)) > tol)
          {
            return std::nullopt;
          }

        // Project s2 endpoints onto s1's direction.
        const double t1 = detail::project(s1, s2.p1);
        const double t2 = detail::project(s1, s2.p2);

        const double lo = std::min(t1, t2);
        const double hi = std::max(t1, t2);

        // s1 occupies [0, len].
        if (lo > len + tol || hi < -tol)
          return std::nullopt;

        // Determine extreme endpoints using the
        // original points, without projection drift.
        struct Candidate
        {
          double                 t;
          decltype(s1.p1) const *point;
        };

        std::array<std::pair<double, Point>, 4> endpoints{
          {{0.0, s1.p1}, {len, s1.p2}, {t1, s2.p1}, {t2, s2.p2}}};

        auto [min_it, max_it] =
          std::minmax_element(endpoints.begin(),
                              endpoints.end(),
                              [](const auto &a, const auto &b) {
                                return a.first < b.first;
                              });

        return Segment{.p1 = min_it->second, .p2 = max_it->second};
      };

      // Copy into a vector so we can merge and erase
      // without depending on std::set ordering.
      std::vector<Segment> work;

      for (const auto &s : segments)
        {
          if (length(s) > tol)
            work.push_back(s);
        }

      bool changed = true;

      while (changed)
        {
          changed = false;

          for (auto i = 0; i < work.size(); ++i)
            {
              for (auto j = i + 1; j < work.size(); ++j)
                {
                  auto merged = try_merge(work[i], work[j], tol);

                  if (!merged)
                    continue;

                  work[i] = *merged;
                  work.erase(work.begin() + j);

                  changed = true;
                  break;
                }

              if (changed)
                break;
            }
        }
      return std::set<Segment>{work.begin(), work.end()};
    }

    [[nodiscard]]
    inline std::vector<Transform>
    generate_candidates(const std::set<Segment> &source,
                        const std::set<Segment> &target,
                        const ToleranceConfig   &tolerance_config)
    {
      std::vector<Transform> candidates;

      const double tol = tolerance_config.tol;

      auto length = [](const Segment &s) -> double {
        return std::hypot(s.p2.x - s.p1.x, s.p2.y - s.p1.y);
      };

      // Generate transformations for every compatible pair.
      for (const Segment &s : source)
        {
          const double source_length = length(s);

          if (source_length <= tol)
            continue;

          for (const Segment &t : target)
            {
              const double target_length = length(t);

              if (std::abs(source_length - target_length) > tol)
                continue;

              // Corresponding endpoint orientations
              candidates.push_back(
                detail::make_transform(s.p1, s.p2, t.p1, t.p2));

              // Reversed endpoint correspondence
              candidates.push_back(
                detail::make_transform(s.p1, s.p2, t.p2, t.p1));
            }
        }

      return candidates;
    }

    [[nodiscard]]
    inline bool
    verify_equal(const std::set<Segment> &source,
                 const std::set<Segment> &target,
                 const Transform         &transform,
                 double                   tol)
    {
      auto covers = [tol](const Segment           &segment,
                          const std::set<Segment> &segments) -> bool {
        const double dx = segment.p2.x - segment.p1.x;
        const double dy = segment.p2.y - segment.p1.y;

        const double len = std::hypot(dx, dy);

        if (len <= tol)
          return false;

        const double ux = dx / len;
        const double uy = dy / len;

        // Collect intervals of segments lying on
        // the supporting line of 'segment'.
        std::vector<std::pair<double, double>> intervals;

        for (const Segment &candidate : segments)
          {
            auto signed_distance = [&](const Point &p) {
              return ux * (p.y - segment.p1.y) - uy * (p.x - segment.p1.x);
            };

            if (std::abs(signed_distance(candidate.p1)) > tol ||
                std::abs(signed_distance(candidate.p2)) > tol)
              {
                continue;
              }

            auto project = [&](const Point &p) {
              return ux * (p.x - segment.p1.x) + uy * (p.y - segment.p1.y);
            };

            double t1 = project(candidate.p1);
            double t2 = project(candidate.p2);

            if (t1 > t2)
              std::swap(t1, t2);

            // Clip the candidate interval to [0, len].
            const double begin = std::max(0.0, t1);
            const double end   = std::min(len, t2);

            if (end >= begin)
              intervals.emplace_back(begin, end);
          }

        if (intervals.empty())
          return false;

        std::ranges::sort(intervals);

        // Check that the union of intervals covers [0, len].
        double covered_until = 0.0;

        for (const auto &[begin, end] : intervals)
          {
            if (begin > covered_until + tol)
              return false;

            covered_until = std::max(covered_until, end);

            if (covered_until >= len - tol)
              return true;
          }

        return false;
      };

      // Transform the source geometry.
      std::set<Segment> transformed_source;

      for (const Segment &segment : source)
        {
          transformed_source.insert(transform.apply(segment));
        }

      // Every transformed source segment must be
      // covered by the target geometry.
      for (const Segment &segment : transformed_source)
        {
          if (!covers(segment, target))
            return false;
        }

      // Every target segment must be covered by
      // the transformed source geometry.
      for (const Segment &segment : target)
        {
          if (!covers(segment, transformed_source))
            return false;
        }

      return true;
    }


    [[nodiscard]]
    inline std::optional<Transform>
    find_transform(
      std::set<Segment> const& source,
      std::set<Segment> const& target,
      ToleranceConfig const& tolerance_config
    ) {
      const double tol = tolerance_config.tol;

      if (source.empty() || target.empty()) {
        if (source.empty() && target.empty()) {
          throw std::invalid_argument("Source and target are empty");
        }
        throw std::invalid_argument("Source or target is empty");
      }

      auto length = [](Segment const& segment) -> double {
        return std::hypot(
          segment.p2.x - segment.p1.x,
          segment.p2.y - segment.p1.y
        );
      };


      // Try each compatible pair, stopping at the first valid transform.
      for (Segment const& s : source) {
        const double source_length = length(s);

        if (source_length <= tol) {
          throw std::invalid_argument("find_transform: Length of a segment in source segments is <= tolerance. This is not allowed.");
        }

        // Compare a segment from source with each segment from target.
        // If length is equal construct the transform and check if the source
        // and target sets are equal under it. If they are, return early with
        // the transform.
        for (Segment const& t : target) {

          // Filter: Check if lengths are the same, they have to be.
          // In canonical form of the segments set, each segment will have a
          // matching segment in the other set, but they are not necessarily
          // ordered the same so we check for each.
          if (std::abs(source_length - length(t)) > tol) {
            continue;
          }

          // First endpoint correspondence.
          const Transform forward = detail::make_transform(
            s.p1, s.p2, t.p1, t.p2
          );

          if (implementation_vdev::verify_equal(source, target, forward, tol)) {
            return forward;
          }

          // Reversed endpoint correspondence.
          const Transform reversed = detail::make_transform(
            s.p1, s.p2, t.p2, t.p1
          );

          if (implementation_vdev::verify_equal(source, target, reversed, tol)) {
            return reversed;
          }
        }
      }
      return std::nullopt;
    }

  } // namespace implementation_vdev

} // namespace UngaBunga

#endif // UNTITLED1_NORMALIZE_H
