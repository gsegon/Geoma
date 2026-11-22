#include <geoma/geoma.h>
#include <iostream>
#include <numbers>
#include <set>

int
main()
{
  using namespace geoma;

  std::set<Segment> source{{{0, 0}, {3, 0}}, {{3, 0}, {3, 2}}};

  Transform expected{.angle = std::numbers::pi / 2.0, .translation = {5, 3}};

  std::set<Segment> target;
  for (auto const &segment : source)
    {
      target.insert(expected.apply(segment));
    }

  auto result = find_transform(source, target, ToleranceConfig{});

  if (!result)
    {
      std::cout << "No transformation found\n";
      return 1;
    }

  std::cout << "Angle: " << result->angle << '\n';
  std::cout << "Translation: " << result->translation.x << ", "
            << result->translation.y << '\n';

  return 0;
}