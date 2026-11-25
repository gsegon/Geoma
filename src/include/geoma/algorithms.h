#ifndef GEOMA_ALGORITHMS_H
#define GEOMA_ALGORITHMS_H
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace geoma
{
  // Public API - Start

  // User customization - Point
  template <class T>
  struct point_traits;

  // User customization - Segment
  template <class T>
  struct segment_traits;

  // User customization - Transform
  template <class T>
  struct transform_traits;



  template <class T>
  concept PointLike = requires(T const& p, double x, double y) {
    { point_traits<T>::x(p) } -> std::convertible_to<double>;
    { point_traits<T>::y(p) } -> std::convertible_to<double>;
    { point_traits<T>::make(x, y) } -> std::same_as<T>;
  };

  /**
   * SegmentLike concept.
   *
   *
   */
  template <class T>
  concept SegmentLike =
    requires(T const& s, typename segment_traits<T>::point_type p) {
      typename segment_traits<T>::point_type;

      requires PointLike<typename segment_traits<T>::point_type>;

      {
        segment_traits<T>::p1(s)
      } -> std::convertible_to<typename segment_traits<T>::point_type>;

      {
        segment_traits<T>::p2(s)
      } -> std::convertible_to<typename segment_traits<T>::point_type>;

      { segment_traits<T>::make(p, p) } -> std::same_as<T>;
    };

  /**
   * TransformLike concept.
   *
   * Must provide `apply` trait for PointLike and SegmentLike.
   * Must provide `make` trait.
   *
   * SegmentLike points must be PointLike.
   */
  template <class T>
  concept TransformLike =
    requires(T const&                                          t,
             typename transform_traits<T>::point_type const&   p,
             typename transform_traits<T>::segment_type const& s,
             double                                            angle) {
      // Check that transform traits actually define both of these
      typename transform_traits<T>::point_type;
      typename transform_traits<T>::segment_type;

      // Make sure they satisfy the PointLike and SegmentLike concepts
      requires PointLike<typename transform_traits<T>::point_type>;
      requires SegmentLike<typename transform_traits<T>::segment_type>;

      // Make sure that point type of segment type of transform
      // is same as tranform point type
      requires std::same_as<typename segment_traits<typename transform_traits<
                              T>::segment_type>::point_type,
                            typename transform_traits<T>::point_type>;

      // Check the transformation can transform points
      {
        transform_traits<T>::apply(t, p)
      } -> std::same_as<typename transform_traits<T>::point_type>;

      // Check the transformation can transform segments
      {
        transform_traits<T>::apply(t, s)
      } -> std::same_as<typename transform_traits<T>::segment_type>;

      // Check that the transform can be constructed from angle and point.
      { transform_traits<T>::make(angle, p) } -> std::same_as<T>;
    };

  /**
   * Config parameters.
   *
   */
  struct ToleranceConfig
  {
    double tolerance = 1e-10;
  };



  /**
   * Check for invalid segments, do any preconditioning that might be beneficial
   * for performance. Can be no-op initially.
   *
   * @param segments
   * @param tolerance_config
   * @return
   */
  // [[nodiscard]]
  // std::vector<Segment>
  // normalize(std::vector<Segment>   segments,
  //           ToleranceConfig const& tolerance_config);

  /**
   * Convert a given set of segments into a set of segments where neighbouring
   * colinear segments are merged into one. Resulting set must represent the
   * same geometric union as the input set.
   * @param segments
   * @param tolerance_config
   * @return
   */
  template <SegmentLike S>
  [[nodiscard]]
  std::vector<S>
  canonicalize(std::span<S const>     segments,
               ToleranceConfig const& tolerance_config);


  /**
   * Check if `target = transform(source)`.
   *
   * @param source
   * @param target
   * @param transform
   * @param tol
   * @return
   */
  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  bool
  verify_equal(std::span<S const> source,
               std::span<S const> target,
               TR const&          transform,
               double             tol);



  /**
   *
   * @tparam S
   * @tparam TR
   * @param source
   * @param target
   * @param tolerance_config
   * @return
   */
  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  std::optional<TR>
  find_transform(std::span<S const>     source,
                 std::span<S const>     target,
                 ToleranceConfig const& tolerance_config);


  /**
   *
   * @tparam S
   * @tparam TR
   * @param source
   * @param target
   * @param tolerance_config
   * @return
   */
  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  std::optional<TR>
  find_transform(std::vector<S> const&  source,
                 std::vector<S> const&  target,
                 ToleranceConfig const& tolerance_config);


  /**
   * Detect if the set of segments has a rotational symmetry.
   * @param segments
   * @param tol
   * @return
   */
  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  inline bool
  has_rotational_symmetry(std::span<S const> segments, double tol);

  // Public API - End

  /**
   * Implementation declaration sections
   */
  namespace implementation_vdev
  {
    // inline std::vector<Segment>
    // normalize(std::vector<Segment>   segments,
    //           ToleranceConfig const& tolerance_config);

    template <SegmentLike S>
    [[nodiscard]]
    std::vector<S>
    canonicalize(std::span<S const>     segments,
                 ToleranceConfig const& tolerance_config);


    template <SegmentLike S, TransformLike TR>
    [[nodiscard]]
    bool
    verify_equal(std::span<S const> source,
                 std::span<S const> target,
                 TR const&          transform,
                 double             tol);

    template <SegmentLike S, TransformLike TR>
    [[nodiscard]]
    std::optional<TR>
    find_transform(std::span<S const>     source,
                   std::span<S const>     target,
                   ToleranceConfig const& tolerance_config);

    template <SegmentLike S, TransformLike TR>
    [[nodiscard]]
    inline bool
    has_rotational_symmetry(std::span<S const> segments, double tol);


  } // namespace implementation_vdev



  // Public API definitions - Start

  template <SegmentLike S>
  [[nodiscard]]
  inline std::vector<S>
  canonicalize(std::span<S const>     segments,
               ToleranceConfig const& tolerance_config)
  {
    return implementation_vdev::canonicalize<S>(segments, tolerance_config);
  }

  template <SegmentLike S>
  [[nodiscard]]
  inline std::vector<S>
  canonicalize(std::vector<S> const&  segments,
               ToleranceConfig const& tolerance_config)
  {
    return canonicalize<S>(std::span<S const>{segments}, tolerance_config);
  }


  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  inline bool
  verify_equal(std::span<S const> source,
               std::span<S const> target,
               TR const&          transform,
               double             tol)
  {
    return implementation_vdev::verify_equal<S, TR>(source,
                                                    target,
                                                    transform,
                                                    tol);
  }

  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  inline bool
  verify_equal(std::vector<S> const& source,
               std::vector<S> const& target,
               TR const&             transform,
               double                tolerance)
  {
    return verify_equal<S, TR>(std::span<S const>{source},
                               std::span<S const>{target},
                               transform,
                               tolerance);
  }

  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  inline std::optional<TR>
  find_transform(std::span<S const>     source,
                 std::span<S const>     target,
                 ToleranceConfig const& tolerance_config)
  {
    return implementation_vdev::find_transform<S, TR>(source,
                                                      target,
                                                      tolerance_config);
  }

  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  inline std::optional<TR>
  find_transform(std::vector<S> const&  source,
                 std::vector<S> const&  target,
                 ToleranceConfig const& tolerance_config)
  {
    return find_transform<S, TR>(std::span<S const>{source},
                                 std::span<S const>{target},
                                 tolerance_config);
  }


  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  inline bool
  has_rotational_symmetry(std::span<S const> segments, double tol)
  {
    return implementation_vdev::has_rotational_symmetry<S, TR>(segments, tol);
  }

  template <SegmentLike S, TransformLike TR>
  [[nodiscard]]
  inline bool
  has_rotational_symmetry(std::vector<S> const& segments, double tol)
  {
    return implementation_vdev::has_rotational_symmetry<S, TR>(
      std::span<S const>{segments}, tol);
  }


  // Public API definitions - End



  // Implementation
  namespace implementation_vdev
  {

    namespace detail
    {

      /**
       * Calculate the parametric value of the point P on segment S.
       * 0 - if P is at the P1 of the segment
       * 1 - if P is at the P2 of the segment
       * Linear in between.
       * @tparam S Segment
       * @tparam P Point
       */
      template <SegmentLike S, PointLike P>
      inline auto project = [](S const& s, P const& p) -> double {
        using ST = segment_traits<S>;
        using PT = point_traits<P>;

        double const dx = PT::x(ST::p2(s)) - PT::x(ST::p1(s));
        double const dy = PT::y(ST::p2(s)) - PT::y(ST::p1(s));

        double const len = std::hypot(dx, dy);

        if (len == 0.0)
          throw std::invalid_argument("project: degenerate segment");

        return (dx * (PT::x(p) - PT::x(ST::p1(s))) +
                dy * (PT::y(p) - PT::y(ST::p1(s)))) /
               len;
      };

      /**
       * make a transform between the would be segments
       *  {source_p1, source_p2}
       *  {target_p1, target_p2}
       *
       * @tparam P
       */
      template <PointLike P, TransformLike TR>
      inline auto make_transform = [](P const& source_p1,
                                      P const& source_p2,
                                      P const& target_p1,
                                      P const& target_p2) -> TR {
        using PT  = point_traits<P>;
        using TRT = transform_traits<TR>;

        double const sx = PT::x(source_p2) - PT::x(source_p1);
        double const sy = PT::y(source_p2) - PT::y(source_p1);

        double const tx = PT::x(target_p2) - PT::x(target_p1);
        double const ty = PT::y(target_p2) - PT::y(target_p1);

        double const angle = std::atan2(ty, tx) - std::atan2(sy, sx);

        double const c = std::cos(angle);
        double const s = std::sin(angle);

        P const translation{.x = PT::x(target_p1) -
                                 (c * PT::x(source_p1) - s * PT::y(source_p1)),
                            .y = PT::y(target_p1) -
                                 (s * PT::x(source_p1) + c * PT::y(source_p1))};

        return TRT::make(angle, translation);
      };


      /**
       * Struct to hold segment bounds for BVH
       */
      struct SegmentBounds
      {
        double min_x;
        double min_y;
        double max_x;
        double max_y;
      };

      template <SegmentLike S>
      inline SegmentBounds
      segment_bounds(S const& s)
      {
        using ST = segment_traits<S>;
        using PT = point_traits<typename ST::point_type>;

        return {.min_x = std::min(PT::x(ST::p1(s)), PT::x(ST::p2(s))),
                .min_y = std::min(PT::y(ST::p1(s)), PT::y(ST::p2(s))),
                .max_x = std::max(PT::x(ST::p1(s)), PT::x(ST::p2(s))),
                .max_y = std::max(PT::y(ST::p1(s)), PT::y(ST::p2(s)))};
      }

      /**
       * Checks if there is an overlap of 2 boxes in 2D.
       * @param a
       * @param b
       * @return
       */
      inline bool
      bounds_intersect(SegmentBounds const& a, SegmentBounds const& b)
      {
        return

          // True if either of these:
          // a----b--a--b
          // a--b--b----A
          // i.e. checks if there is an intersection
          a.min_x <= b.max_x && b.min_x <= a.max_x &&

          // Same as above but for y
          a.min_y <= b.max_y && b.min_y <= a.max_y;
      }

      /**
       *  Bounding volume hierarchy for fast spatial checks.
       *
       *
       *
       * @tparam S
       */
      template <SegmentLike S>
      struct SegmentBvh
      {
        static constexpr std::size_t npos = static_cast<std::size_t>(-1);

        /**
         * Node struct
         */
        struct Node
        {
          SegmentBounds box{}; // 2D box containing all segments of this node
          std::size_t   begin_seg_index{};       // Begin index  (inclusive)
          std::size_t   end_seg_index{};         // End index    (exclusive)
          std::size_t   left_child_node  = npos; // left child node index
          std::size_t   right_child_node = npos; // right child node index
        };

        std::span<S const>         segments;
        std::vector<SegmentBounds> boxes;
        std::vector<std::size_t>   order;
        std::vector<Node>          nodes;

        explicit SegmentBvh(std::vector<S> const& input)
        : SegmentBvh(std::span<S const>{input})
        {}

        explicit SegmentBvh(std::span<S const> input)
          : segments(input)
        {
          // Resrve to avoid reallocations.
          boxes.reserve(segments.size());

          // Construct the boxes.
          for (auto const& s : segments)
            {
              boxes.push_back(segment_bounds(s));
            }

          order.resize(segments.size());
          std::iota(order.begin(), order.end(), std::size_t{0});

          if (!order.empty())
            {
              // Build the BVH
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

          auto const width  = box.max_x - box.min_x;
          auto const height = box.max_y - box.min_y;

          bool const use_x = width >= height;

          auto const mid = std::midpoint(begin, end); // Median center

          // Puts the value at mid to whatever it would be if the order was
          // fully ordered and also all elements before mid are less then or
          // equal to elements after n-th, as given by the predicate. Complexity
          // should be O(n), so O(n log(n)) for recursive call.
          std::nth_element(
            order.begin() + begin, // Beginning of range
            order.begin() + mid,   // Position of the nth element
            order.begin() + end,   // End of range

            // Returns true if the center of bounding box `a` precedes
            // the center of bounding box `b` along the selected axis.
            [&](std::size_t a, std::size_t b) {
              auto const& box_a_bounds = boxes[a];
              auto const& box_b_bounds = boxes[b];

              if (use_x)
                return std::midpoint(box_a_bounds.min_x, box_a_bounds.max_x) <
                       std::midpoint(box_b_bounds.min_x, box_b_bounds.max_x);
              else
                return std::midpoint(box_a_bounds.min_y, box_a_bounds.max_y) <
                       std::midpoint(box_b_bounds.min_y, box_b_bounds.max_y);
            });

          auto const left  = build(begin, mid); // build on left
          auto const right = build(mid, end);  // build on right

          nodes[node].left_child_node  = left;
          nodes[node].right_child_node = right;

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

          // Filter through full groups of nodes.
          if (!bounds_intersect(node.box, box))
            {
              return;
            }

          // Once here we are at the leaf node.
          // Run the function on any segments of this box.
          if (node.left_child_node == npos)
            {
              for (std::size_t i = node.begin_seg_index; i < node.end_seg_index;
                   ++i)
                {
                  auto const index = order[i];

                  if (bounds_intersect(boxes[index], box))
                    {
                      fn(segments[index]);
                    }
                }
              return;
            }

          visit(node.left_child_node, box, fn);
          visit(node.right_child_node, box, fn);
        }
      };

      // Check whether every segment in source is covered by
      // a union of approximately collinear target segments.
      template <SegmentLike S>
      inline bool
      covered_by(std::span<S const>   source,
                 SegmentBvh<S> const& target,
                 double               tol)
      {
        std::vector<std::pair<double, double>> intervals;

        for (auto const& s : source)
          {
            using ST = segment_traits<S>;
            using P  = ST::point_type;
            using PT = point_traits<P>;

            P const p1 = ST::p1(s);
            P const p2 = ST::p2(s);

            double const x1 = PT::x(p1);
            double const y1 = PT::y(p1);
            double const x2 = PT::x(p2);
            double const y2 = PT::y(p2);

            double const dx  = x2 - x1;
            double const dy  = y2 - y1;
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

            // Enlarge by tolerance
            box.min_x -= tol;
            box.min_y -= tol;
            box.max_x += tol;
            box.max_y += tol;

            // The function will be called on every segment within the box.
            target.query(box, [&](S const& t) {
              P const pt1 = ST::p1(t);
              P const pt2 = ST::p2(t);

              double const t_x1 = PT::x(pt1);
              double const t_y1 = PT::y(pt1);
              double const t_x2 = PT::x(pt2);
              double const t_y2 = PT::y(pt2);

              double const ax = t_x1 - x1;
              double const ay = t_y1 - y1;
              double const bx = t_x2 - x1;
              double const by = t_y2 - y1;

              // Both endpoints must lie within tolerance
              // of the supporting line of s.
              double const distance_a = std::abs(ux * ay - uy * ax);
              double const distance_b = std::abs(ux * by - uy * bx);

              if (distance_a > tol || distance_b > tol)
                {
                  return;   // reject
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


            // The source segment is project on [0, len], and target segments
            // are projected on the same coordinate system and stored as [lo, hi]
            // intervals. Intervals are sorted.
            // Accumulate covered_end
            double covered_end = 0.0;
            for (auto const& [lo, hi] : intervals)
              {

                // New lo must not be bigger the current covered_end (within tol)
                if (lo > covered_end + tol)
                  {
                    return false;
                  }

                covered_end = std::max(covered_end, hi);

                // If covered_end is bigger then len (within tolerance)
                // we exit early and go to next source segment.
                if (covered_end >= len - tol)
                  {
                    break;
                  }
              }

            // After accumulation, if end length of intervals is lower then len
            // of the source segment (within tol), we return early with false.
            if (covered_end < len - tol)
              {
                return false;
              }
          }

        // All source segments have target segments whose intervals cover them
        // fully.
        return true;
      }

      template <SegmentLike S>
      inline bool
      covered_by(std::vector<S> const& source,
                 SegmentBvh<S> const&  target,
                 double                tol)
      {
        return covered_by(std::span<S const>(source), target, tol);
      }


    } // namespace detail

    // // No-op for now
    // inline std::vector<Segment>
    // normalize(std::vector<Segment>   segments,
    //           ToleranceConfig const& tolerance_config)
    // {
    //   return segments;
    // }



    template <SegmentLike S>
    [[nodiscard]]
    std::vector<S>
    canonicalize(std::span<S const>     segments,
                 ToleranceConfig const& tolerance_config)
    {
      using ST = segment_traits<S>;
      using P  = ST::point_type;
      using PT = point_traits<P>;

      double const tol = tolerance_config.tolerance;

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
        P      p_lo;
        P      p_hi;
      };

      std::unordered_map<LineKey, std::vector<Interval>, LineHash> groups;

      // Group segments by supporting line.
      for (auto const& s : segments)
        {
          P const p1 = ST::p1(s);
          P const p2 = ST::p2(s);

          double const x1 = PT::x(p1);
          double const y1 = PT::y(p1);
          double const x2 = PT::x(p2);
          double const y2 = PT::y(p2);

          double const dx  = x2 - x1;
          double const dy  = y2 - y1;
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

          // Canonicalize the supporting line direction.
          if (ux < 0.0 || (ux == 0.0 && uy < 0.0))
            {
              ux = -ux;
              uy = -uy;
            }

          double const offset = -uy * x1 + ux * y1;

          // Use the dominant coordinate as the interval axis.
          bool const use_x = std::abs(ux) >= std::abs(uy);

          double const t1 = use_x ? x1 : y1;
          double const t2 = use_x ? x2 : y2;

          if (!std::isfinite(offset) || !std::isfinite(t1) ||
              !std::isfinite(t2))
            {
              throw std::invalid_argument("canonicalize: invalid coordinates");
            }

          auto& intervals = groups[LineKey{ux, uy, offset}];

          if (t1 <= t2)
            {
              intervals.push_back({t1, t2, p1, p2});
            }
          else
            {
              intervals.push_back({t2, t1, p2, p1});
            }
        }

      std::vector<S> result;
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
                  result.push_back(ST::make(merged.p_lo, merged.p_hi));

                  merged = next;
                }
            }

          result.push_back(ST::make(merged.p_lo, merged.p_hi));
        }

      return result;
    }



    template <SegmentLike S, TransformLike TR>
    [[nodiscard]]
    inline bool
    verify_equal(std::span<S const> source,
                 std::span<S const> target,
                 TR const&          transform,
                 double             tol)
    {
      if (!std::isfinite(tol) || tol < 0.0)
        {
          throw std::invalid_argument("verify_equal: invalid tolerance");
        }

      // If they are both empty, they are equal.
      // Otherwise, they are not.
      if (source.empty() || target.empty())
        {
          return source.empty() && target.empty();
        }

      // Transform the already canonicalized source.
      std::vector<S> transformed;
      transformed.reserve(source.size());

      // Transform the source segments
      for (auto const& segment : source)
        {
          transformed.push_back(transform.apply(segment));
        }

      // Construct the BVH for the target segments.
      detail::SegmentBvh target_index(target);

      // Check that evey segment of `transformed` is covered by `target`
      if (!detail::covered_by(transformed, target_index, tol))
        {
          return false;
        }

      // Construct the BVH for the transformed segments.
      detail::SegmentBvh transformed_index(transformed);

      // Check that evey segment of `target` is covered by `transformed`
      return detail::covered_by(target, transformed_index, tol);
    }


    template <SegmentLike S, TransformLike TR>
    [[nodiscard]]
    inline std::optional<TR>
    find_transform(std::span<S const>     source,
                   std::span<S const>     target,
                   ToleranceConfig const& tolerance_config)
    {
      using ST = segment_traits<S>;
      using P  = ST::point_type;
      using PT = point_traits<P>;

      auto const canonical_source =
        implementation_vdev::canonicalize(source, tolerance_config);

      auto const canonical_target =
        implementation_vdev::canonicalize(target, tolerance_config);

      double const tol = tolerance_config.tolerance;

      if (source.empty() || target.empty())
        {
          if (source.empty() && target.empty())
            {
              throw std::invalid_argument("Source and target are empty");
            }
          throw std::invalid_argument("Source or target is empty");
        }

      auto length = [](S const& s) -> double {
        P const p1 = ST::p1(s);
        P const p2 = ST::p2(s);

        double const x1 = PT::x(p1);
        double const y1 = PT::y(p1);
        double const x2 = PT::x(p2);
        double const y2 = PT::y(p2);

        double const dx  = x2 - x1;
        double const dy  = y2 - y1;
        double const len = std::hypot(dx, dy);

        return len;
      };


      // Try each compatible pair, stopping at the first valid transform.
      for (S const& s : source)
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
          for (S const& t : target)
            {
              // Filter: Check if lengths are the same, they have to be.
              // In canonical form of the segments set, each segment will have a
              // matching segment in the other set, but they are not necessarily
              // ordered the same so we check for each.
              if (std::abs(source_length - length(t)) > tol)
                {
                  continue;
                }

              // Points
              P const ps1 = ST::p1(s);
              P const ps2 = ST::p2(s);

              P const pt1 = ST::p1(t);
              P const pt2 = ST::p2(t);

              // First endpoint correspondence.
              TR const forward =
                detail::make_transform<P, TR>(ps1, ps2, pt1, pt2);

              if (implementation_vdev::verify_equal(
                    source, target, forward, tol))
                {
                  return forward;
                }

              // Reversed endpoint correspondence.

              TR const reversed =
                detail::make_transform<P, TR>(ps1, ps2, pt2, pt1);

              if (implementation_vdev::verify_equal(
                    source, target, reversed, tol))
                {
                  return reversed;
                }
            }
        }
      return std::nullopt;
    }

    template <SegmentLike S, TransformLike TR>
    [[nodiscard]]
    inline bool
    has_rotational_symmetry(std::span<S const> segments, double tol)
    {
      if (!std::isfinite(tol) || tol < 0.0)
        {
          throw std::invalid_argument("invalid tolerance");
        }

      if (segments.empty())
        {
          return false;
        }

      // Compute the length-weighted geometric center.
      using ST = segment_traits<S>;
      using P  = ST::point_type;
      using PT = point_traits<P>;

      P center = PT::make(0, 0);

      double total_length = 0.0;

      for (auto const& s : segments)
        {
          double const len = std::hypot(PT::x(ST::p2(s)) - PT::x(ST::p1(s)),
                                        PT::y(ST::p2(s)) - PT::y(ST::p1(s)));

          if (!std::isfinite(len))
            {
              throw std::invalid_argument("invalid segment");
            }

          if (len <= tol)
            {
              continue;
            }

          center.x += (PT::x(ST::p1(s)) + PT::x(ST::p2(s))) * 0.5 * len;
          center.y += (PT::y(ST::p1(s)) + PT::y(ST::p2(s))) * 0.5 * len;
          total_length += len;
        }

      if (total_length <= 0.0)
        {
          return false;
        }

      center.x /= total_length;
      center.y /= total_length;

      // Find an endpoint farthest from the center.
      P      anchor = PT::make(0, 0);
      double radius = 0.0;

      for (auto const& s : segments)
        {
          for (auto const& p : {ST::p1(s), ST::p2(s)})
            {
              double const r = std::hypot(p.x - center.x, p.y - center.y);

              if (r > radius)
                {
                  radius = r;
                  anchor = p;
                }
            }
        }

      if (radius <= tol)
        {
          return false;
        }

      double const ax = anchor.x - center.x;
      double const ay = anchor.y - center.y;

      // Test rotations mapping the anchor to other endpoints.
      for (auto const& [p1, p2] : segments)
        {
          for (auto const& p : {p1, p2})
            {
              double const bx = p.x - center.x;
              double const by = p.y - center.y;

              if (std::abs(std::hypot(bx, by) - radius) > tol)
                {
                  continue;
                }

              double const angle =
                std::atan2(ax * by - ay * bx, ax * bx + ay * by);

              // Ignore the identity rotation.
              if (std::abs(angle) <= std::max(1e-12, tol / radius))
                {
                  continue;
                }

              double const c  = std::cos(angle);
              double const sn = std::sin(angle);

              TR transform{
                .angle       = angle,
                .translation = {center.x - (c * center.x - sn * center.y),
                                center.y - (sn * center.x + c * center.y)}};

              if (implementation_vdev::verify_equal(
                    segments, segments, transform, tol))
                {
                  return true;
                }
            }
        }

      return false;
    }


  } // namespace implementation_vdev

} // namespace geoma

#endif // GEOMA_ALGORITHMS_H
