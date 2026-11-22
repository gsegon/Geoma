#ifndef GEOMA_ALGORITHMS_H
#define GEOMA_ALGORITHMS_H
#include <algorithm>
#include <array>
#include <cmath>
#include <geoma/config.h>
#include <geoma/primitives.h>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace geoma
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
  std::vector<Segment>
  normalize(std::vector<Segment>   segments,
            ToleranceConfig const& tolerance_config);

  /**
   * Convert a given set of segments into a set of segments where neighbouring
   * colinear segments are merged into one. Resulting set must represent the
   * same geometric union as the input set.
   * @param segments
   * @param tolerance_config
   * @return
   */
  [[nodiscard]]
  std::vector<Segment>
  canonicalize(std::span<Segment const> segments,
               ToleranceConfig const&   tolerance_config);


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
  verify_equal(std::set<Segment> const& source,
               std::set<Segment> const& target,
               Transform const&         transform,
               double                   tol);


  /**
   * Find the transform
   * @param source
   * @param target
   * @param tolerance_config
   * @return
   */
  [[nodiscard]]
  std::optional<Transform>
  find_transform(std::span<Segment const> source,
                 std::span<Segment const> target,
                 ToleranceConfig const&   tolerance_config);

  // Public API - End

  namespace implementation_vdev
  {
    inline std::vector<Segment>
    normalize(std::vector<Segment>   segments,
              ToleranceConfig const& tolerance_config);

    [[nodiscard]]
    std::vector<Segment>
    canonicalize(std::span<Segment const> segments,
                 ToleranceConfig const&   tolerance_config);

    [[nodiscard]]
    std::vector<Transform>
    generate_candidates(std::set<Segment> const& source,
                        std::set<Segment> const& target,
                        ToleranceConfig const&   tolerance_config);

    [[nodiscard]]
    bool
    verify_equal(std::span<Segment const> source,
                 std::span<Segment const> target,
                 Transform const&         transform,
                 double                   tol);

    [[nodiscard]]
    std::optional<Transform>
    find_transform(std::span<Segment const> source,
                   std::span<Segment const> target,
                   ToleranceConfig const&   tolerance_config);

  } // namespace implementation_vdev

  [[nodiscard]]
  inline std::vector<Segment>
  canonicalize(std::span<Segment const> segments,
               ToleranceConfig const&   tolerance_config)
  {
    return implementation_vdev::canonicalize(segments, tolerance_config);
  }

  [[nodiscard]]
  inline std::vector<Transform>
  generate_candidates(std::set<Segment> const& source,
                      std::set<Segment> const& target,
                      ToleranceConfig const&   tolerance_config)
  {
    return implementation_vdev::generate_candidates(source,
                                                    target,
                                                    tolerance_config);
  }

  [[nodiscard]]
  inline bool
  verify_equal(std::span<Segment const>  source,
               std::span<Segment const>  target,
               Transform const&         transform,
               double                   tol)
  {
    return implementation_vdev::verify_equal(source, target, transform, tol);
  }

  [[nodiscard]]
  inline std::optional<Transform>
  find_transform(std::span<Segment const> source,
                 std::span<Segment const> target,
                 ToleranceConfig const&   tolerance_config)
  {
    return implementation_vdev::find_transform(source,
                                               target,
                                               tolerance_config);
  }

  namespace implementation_vdev
  {

    namespace detail
    {
      // Signed perpendicular distances to the
      // supporting line of s1.
      inline auto signed_distance = [](Segment const& s,
                                       Point const&   p) -> double {
        double const dx = s.p2.x - s.p1.x;
        double const dy = s.p2.y - s.p1.y;

        double const len = std::hypot(dx, dy);

        if (len == 0.0)
          throw std::invalid_argument("signed_distance: degenerate segment");

        return (dx * (p.y - s.p1.y) - dy * (p.x - s.p1.x)) / len;
      };

      inline auto project = [](Segment const& s, Point const& p) -> double {
        double const dx = s.p2.x - s.p1.x;
        double const dy = s.p2.y - s.p1.y;

        double const len = std::hypot(dx, dy);

        if (len == 0.0)
          throw std::invalid_argument("project: degenerate segment");

        return (dx * (p.x - s.p1.x) + dy * (p.y - s.p1.y)) / len;
      };

      inline auto make_transform = [](Point const& source_p1,
                                      Point const& source_p2,
                                      Point const& target_p1,
                                      Point const& target_p2) -> Transform {
        double const sx = source_p2.x - source_p1.x;
        double const sy = source_p2.y - source_p1.y;

        double const tx = target_p2.x - target_p1.x;
        double const ty = target_p2.y - target_p1.y;

        double const angle = std::atan2(ty, tx) - std::atan2(sy, sx);

        double const c = std::cos(angle);
        double const s = std::sin(angle);

        Point translation{target_p1.x - (c * source_p1.x - s * source_p1.y),
                          target_p1.y - (s * source_p1.x + c * source_p1.y)};

        return Transform{.angle = angle, .translation = translation};
      };


      struct SegmentBounds
      {
        double min_x;
        double min_y;
        double max_x;
        double max_y;
      };

      inline SegmentBounds
      segment_bounds(Segment const& s)
      {
        return {std::min(s.p1.x, s.p2.x),
                std::min(s.p1.y, s.p2.y),
                std::max(s.p1.x, s.p2.x),
                std::max(s.p1.y, s.p2.y)};
      }

      inline bool
      bounds_intersect(SegmentBounds const& a, SegmentBounds const& b)
      {
        return a.min_x <= b.max_x && b.min_x <= a.max_x && a.min_y <= b.max_y &&
               b.min_y <= a.max_y;
      }

      // Bounding-volume hierarchy for fast spatial queries.
      struct SegmentBvh
      {
        static constexpr std::size_t npos = static_cast<std::size_t>(-1);

        struct Node
        {
          SegmentBounds box;
          std::size_t   begin;
          std::size_t   end;
          std::size_t   left  = npos;
          std::size_t   right = npos;
        };

        std::span<Segment const>   segments;
        std::vector<SegmentBounds> boxes;
        std::vector<std::size_t>   order;
        std::vector<Node>          nodes;

        explicit SegmentBvh(std::span<Segment const> input)
          : segments(input)
        {
          boxes.reserve(segments.size());

          for (auto const& s : segments)
            {
              boxes.push_back(segment_bounds(s));
            }

          order.resize(segments.size());
          std::iota(order.begin(), order.end(), std::size_t{0});

          if (!order.empty())
            {
              build(0, order.size());
            }
        }

        std::size_t
        build(std::size_t begin, std::size_t end)
        {
          auto box = boxes[order[begin]];

          for (std::size_t i = begin + 1; i < end; ++i)
            {
              auto const& b = boxes[order[i]];

              box.min_x = std::min(box.min_x, b.min_x);
              box.min_y = std::min(box.min_y, b.min_y);
              box.max_x = std::max(box.max_x, b.max_x);
              box.max_y = std::max(box.max_y, b.max_y);
            }

          auto const node = nodes.size();
          nodes.push_back({box, begin, end});

          // Small leaves avoid unnecessary tree traversal.
          if (end - begin <= 8)
            {
              return node;
            }

          bool const use_x = box.max_x - box.min_x >= box.max_y - box.min_y;

          auto const mid = begin + (end - begin) / 2;

          std::nth_element(order.begin() + begin,
                           order.begin() + mid,
                           order.begin() + end,
                           [&](std::size_t a, std::size_t b) {
                             auto const& x = boxes[a];
                             auto const& y = boxes[b];

                             return use_x ? std::midpoint(x.min_x, x.max_x) <
                                              std::midpoint(y.min_x, y.max_x) :
                                            std::midpoint(x.min_y, x.max_y) <
                                              std::midpoint(y.min_y, y.max_y);
                           });

          auto const left  = build(begin, mid);
          auto const right = build(mid, end);

          nodes[node].left  = left;
          nodes[node].right = right;

          return node;
        }

        template <class Fn>
        void
        query(SegmentBounds const& box, Fn&& fn) const
        {
          if (!nodes.empty())
            {
              visit(0, box, fn);
            }
        }

        template <class Fn>
        void
        visit(std::size_t node_id, SegmentBounds const& box, Fn& fn) const
        {
          auto const& node = nodes[node_id];

          if (!bounds_intersect(node.box, box))
            {
              return;
            }

          if (node.left == npos)
            {
              for (std::size_t i = node.begin; i < node.end; ++i)
                {
                  auto const index = order[i];

                  if (bounds_intersect(boxes[index], box))
                    {
                      fn(segments[index]);
                    }
                }
              return;
            }

          visit(node.left, box, fn);
          visit(node.right, box, fn);
        }
      };

      // Check whether every segment in source is covered by
      // a union of approximately collinear target segments.
      inline bool
      covered_by(std::span<Segment const> source,
                 SegmentBvh const&        target,
                 double                   tol)
      {
        std::vector<std::pair<double, double>> intervals;

        for (auto const& s : source)
          {
            double const dx  = s.p2.x - s.p1.x;
            double const dy  = s.p2.y - s.p1.y;
            double const len = std::hypot(dx, dy);

            if (len <= tol)
              {
                continue;
              }

            double const ux = dx / len;
            double const uy = dy / len;

            intervals.clear();

            // Query segments whose bounding boxes overlap ours.
            auto box = segment_bounds(s);

            box.min_x -= tol;
            box.min_y -= tol;
            box.max_x += tol;
            box.max_y += tol;

            target.query(box, [&](Segment const& t) {
              double const ax = t.p1.x - s.p1.x;
              double const ay = t.p1.y - s.p1.y;
              double const bx = t.p2.x - s.p1.x;
              double const by = t.p2.y - s.p1.y;

              // Both endpoints must lie within tolerance
              // of the supporting line of s.
              if (std::abs(ux * ay - uy * ax) > tol ||
                  std::abs(ux * by - uy * bx) > tol)
                {
                  return;
                }

              // Project onto s, whose interval is [0, len].
              double const t1 = ux * ax + uy * ay;
              double const t2 = ux * bx + uy * by;

              double const lo = std::min(t1, t2);
              double const hi = std::max(t1, t2);

              if (hi >= -tol && lo <= len + tol)
                {
                  intervals.emplace_back(lo, hi);
                }
            });

            std::ranges::sort(intervals);

            double covered_end = 0.0;

            for (auto const& [lo, hi] : intervals)
              {
                if (lo > covered_end + tol)
                  {
                    return false;
                  }

                covered_end = std::max(covered_end, hi);

                if (covered_end >= len - tol)
                  {
                    break;
                  }
              }

            if (covered_end < len - tol)
              {
                return false;
              }
          }

        return true;
      }


    } // namespace detail

    // No-op for now
    inline std::vector<Segment>
    normalize(std::vector<Segment>   segments,
              ToleranceConfig const& tolerance_config)
    {
      return segments;
    }



    [[nodiscard]]
    inline std::vector<Segment>
    canonicalize(std::span<Segment const> segments,
                 ToleranceConfig const&   tolerance_config)
    {
      double const tol = tolerance_config.tol;

      if (!std::isfinite(tol) || tol < 0.0)
        {
          throw std::invalid_argument("canonicalize: invalid tolerance");
        }

      struct LineKey
      {
        double ux;
        double uy;
        double offset;

        bool
        operator==(LineKey const&) const = default;
      };

      struct LineHash
      {
        std::size_t
        operator()(LineKey const& key) const noexcept
        {
          std::size_t h = std::hash<double>{}(key.ux);

          auto combine = [&](double value) {
            h ^= std::hash<double>{}(value) + 0x9e3779b97f4a7c15ULL + (h << 6) +
                 (h >> 2);
          };

          combine(key.uy);
          combine(key.offset);

          return h;
        }
      };

      struct Interval
      {
        double lo;
        double hi;
        Point  p_lo;
        Point  p_hi;
      };

      std::unordered_map<LineKey, std::vector<Interval>, LineHash> groups;

      // Group segments by supporting line.
      for (auto const& s : segments)
        {
          double const dx  = s.p2.x - s.p1.x;
          double const dy  = s.p2.y - s.p1.y;
          double const len = std::hypot(dx, dy);

          if (!std::isfinite(len))
            {
              throw std::invalid_argument("canonicalize: invalid segment");
            }

          if (len <= tol)
            {
              continue;
            }

          double ux = dx / len;
          double uy = dy / len;

          if (ux < 0.0 || (ux == 0.0 && uy < 0.0))
            {
              ux = -ux;
              uy = -uy;
            }

          double const offset = -uy * s.p1.x + ux * s.p1.y;
          bool const   use_x  = std::abs(ux) >= std::abs(uy);

          double const t1 = use_x ? s.p1.x : s.p1.y;
          double const t2 = use_x ? s.p2.x : s.p2.y;

          if (!std::isfinite(offset) || !std::isfinite(t1) ||
              !std::isfinite(t2))
            {
              throw std::invalid_argument("canonicalize: invalid coordinates");
            }

          auto& intervals = groups[LineKey{ux, uy, offset}];

          if (t1 <= t2)
            {
              intervals.push_back({t1, t2, s.p1, s.p2});
            }
          else
            {
              intervals.push_back({t2, t1, s.p2, s.p1});
            }
        }

      std::vector<Segment> result;
      result.reserve(groups.size());

      // Sort and merge intervals within each supporting line.
      for (auto& [line, intervals] : groups)
        {
          std::ranges::sort(intervals, {}, &Interval::lo);

          double const axis_tol =
            tol * std::max(std::abs(line.ux), std::abs(line.uy));

          Interval merged = intervals.front();

          for (std::size_t i = 1; i < intervals.size(); ++i)
            {
              auto const& next = intervals[i];

              if (next.lo <= merged.hi + axis_tol)
                {
                  if (next.hi > merged.hi)
                    {
                      merged.hi   = next.hi;
                      merged.p_hi = next.p_hi;
                    }
                }
              else
                {
                  result.push_back({merged.p_lo, merged.p_hi});
                  merged = next;
                }
            }

          result.push_back({merged.p_lo, merged.p_hi});
        }

      return result;
    }



    [[nodiscard]]
    inline std::vector<Transform>
    generate_candidates(std::set<Segment> const& source,
                        std::set<Segment> const& target,
                        ToleranceConfig const&   tolerance_config)
    {
      std::vector<Transform> candidates;

      double const tol = tolerance_config.tol;

      auto length = [](Segment const& s) -> double {
        return std::hypot(s.p2.x - s.p1.x, s.p2.y - s.p1.y);
      };

      // Generate transformations for every compatible pair.
      for (Segment const& s : source)
        {
          double const source_length = length(s);

          if (source_length <= tol)
            continue;

          for (Segment const& t : target)
            {
              double const target_length = length(t);

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
  inline bool verify_equal(
    std::span<Segment const> source,
    std::span<Segment const> target,
    Transform const& transform,
    double tol
  ) {
      if (!std::isfinite(tol) || tol < 0.0) {
        throw std::invalid_argument("verify_equal: invalid tolerance");
      }

      if (source.empty() || target.empty()) {
        return source.empty() && target.empty();
      }

      // Transform the already canonicalized source.
      std::vector<Segment> transformed;
      transformed.reserve(source.size());

      for (auto const& segment : source) {
        transformed.push_back(transform.apply(segment));
      }

      // Check whether the transformed source is covered by target.
      detail::SegmentBvh target_index(target);

      if (!detail::covered_by(transformed, target_index, tol)) {
        return false;
      }

      // Check whether target is covered by transformed source.
      detail::SegmentBvh source_index(transformed);

      return detail::covered_by(target, source_index, tol);
    }


    [[nodiscard]]
    inline std::optional<Transform>
    find_transform(std::span<Segment const> source,
                   std::span<Segment const> target,
                   ToleranceConfig const&   tolerance_config)
    {
      double const tol = tolerance_config.tol;

      if (source.empty() || target.empty())
        {
          if (source.empty() && target.empty())
            {
              throw std::invalid_argument("Source and target are empty");
            }
          throw std::invalid_argument("Source or target is empty");
        }

      auto length = [](Segment const& segment) -> double {
        return std::hypot(segment.p2.x - segment.p1.x,
                          segment.p2.y - segment.p1.y);
      };


      // Try each compatible pair, stopping at the first valid transform.
      for (Segment const& s : source)
        {
          double const source_length = length(s);

          if (source_length <= tol)
            {
              throw std::invalid_argument(
                "find_transform: Length of a segment in source segments is <= tolerance. This is not allowed.");
            }

          // Compare a segment from source with each segment from target.
          // If length is equal construct the transform and check if the source
          // and target sets are equal under it. If they are, return early with
          // the transform.
          for (Segment const& t : target)
            {
              // Filter: Check if lengths are the same, they have to be.
              // In canonical form of the segments set, each segment will have a
              // matching segment in the other set, but they are not necessarily
              // ordered the same so we check for each.
              if (std::abs(source_length - length(t)) > tol)
                {
                  continue;
                }

              // First endpoint correspondence.
              Transform const forward =
                detail::make_transform(s.p1, s.p2, t.p1, t.p2);

              if (implementation_vdev::verify_equal(
                    source, target, forward, tol))
                {
                  return forward;
                }

              // Reversed endpoint correspondence.
              Transform const reversed =
                detail::make_transform(s.p1, s.p2, t.p2, t.p1);

              if (implementation_vdev::verify_equal(
                    source, target, reversed, tol))
                {
                  return reversed;
                }
            }
        }
      return std::nullopt;
    }

  } // namespace implementation_vdev

} // namespace geoma

#endif // GEOMA_ALGORITHMS_H
