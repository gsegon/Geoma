#ifndef GEOMA_TEST_PRIMITIVE_TYPES_H
#define GEOMA_TEST_PRIMITIVE_TYPES_H


#include <geoma/algorithms.h>


namespace geoma_test
{

  struct Point
  {
    double x{};
    double y{};
  };

  struct Segment
  {
    Point p1;
    Point p2;

    bool
    operator<(Segment const& other) const
    {
      return std::tie(p1.x, p1.y, p2.x, p2.y) <
             std::tie(other.p1.x, other.p1.y, other.p2.x, other.p2.y);
    }
  };

  struct Transform
  {
    double angle{};
    Point  translation;

    [[nodiscard]] Point
    apply(Point const& p) const;

    [[nodiscard]] Segment
    apply(Segment const& s) const;
  };

  [[nodiscard]]
  inline Point
  Transform::apply(Point const& p) const
  {
    double const c = std::cos(angle);
    double const s = std::sin(angle);

    return Point{.x = c * p.x - s * p.y + translation.x,
                 .y = s * p.x + c * p.y + translation.y};
  }

  [[nodiscard]]
  inline Segment
  Transform::apply(Segment const& segment) const
  {
    return Segment{.p1 = apply(segment.p1), .p2 = apply(segment.p2)};
  }


} // namespace geoma_test

namespace geoma
{


  template <>
  struct point_traits<geoma_test::Point>
  {
    using point_type = geoma_test::Point;

    static double
    x(point_type const& p)
    {
      return p.x;
    }

    static double
    y(point_type const& p)
    {
      return p.y;
    }

    static point_type
    make(double x, double y)
    {
      return {x, y};
    }
  };

  template <>
  struct segment_traits<geoma_test::Segment>
  {
    using segment_type = geoma_test::Segment;
    using point_type   = geoma_test::Point;

    static point_type
    p1(segment_type const& s)
    {
      return s.p1;
    }

    static point_type
    p2(segment_type const& s)
    {
      return s.p2;
    }

    static segment_type
    make(point_type p1, point_type p2)
    {
      return {.p1 = p1, .p2 = p2};
    }
  };


  template <>
  struct transform_traits<geoma_test::Transform>
  {
    using transform_type = geoma_test::Transform;
    using point_type     = geoma_test::Point;
    using segment_type   = geoma_test::Segment;

    [[nodiscard]]
    static point_type
    apply(transform_type const& t, point_type const& p)
    {
      return t.apply(p);
    }

    [[nodiscard]]
    static segment_type
    apply(transform_type const& t, segment_type const& s)
    {
      return t.apply(s);
    }

    [[nodiscard]]
    static transform_type
    make(double angle, point_type const& translation)
    {
      return geoma_test::Transform{.angle = angle, .translation = translation};
    }
  };

} // namespace geoma



#endif // GEOMA_TEST_PRIMITIVE_TYPES_H
