#ifndef GEOMA_PRIMITIVES_H
#define GEOMA_PRIMITIVES_H
#include <cmath>
#include <tuple>

namespace geoma
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

    return Point{.x = c * p.x - s * p.y + translation.x,
                 .y = s * p.x + c * p.y + translation.y};
  }

  [[nodiscard]]
  inline Segment
  Transform::apply(const Segment &segment) const
  {
    return Segment
    {
      .p1 = apply(segment.p1),
      .p2 = apply(segment.p2)
    };
  }


} // namespace geoma

#endif // GEOMA_PRIMITIVES_H
