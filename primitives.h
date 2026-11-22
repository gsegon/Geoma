#ifndef UNTITLED1_PRIMITIVES_H
#define UNTITLED1_PRIMITIVES_H
#include <cmath>
#include <tuple>

// UngaBunga!!! Better not end up in production code, gonna Unga all the Bunga!!
namespace UngaBunga
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
    operator<(const Segment &other) const
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
    apply(const Point &p) const;

    [[nodiscard]] Segment
    apply(const Segment &s) const;
  };

  [[nodiscard]]
  inline Point
  Transform::apply(const Point &p) const
  {
    const double c = std::cos(angle);
    const double s = std::sin(angle);

    return Point{c * p.x - s * p.y + translation.x,
                 s * p.x + c * p.y + translation.y};
  }

  [[nodiscard]]
  inline Segment
  Transform::apply(const Segment &segment) const
  {
    return Segment{apply(segment.p1), apply(segment.p2)};
  }


} // namespace UngaBunga

#endif // UNTITLED1_PRIMITIVES_H
