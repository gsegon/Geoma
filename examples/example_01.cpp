#include <geoma/geoma.h>
#include <iostream>
#include <numbers>



int
main()
{
  using namespace geoma;

  std::vector<Segment> source{{{0, 0}, {3, 0}}, {{3, 0}, {3, 2}}};

  Transform expected{.angle = std::numbers::pi / 2.0, .translation = {5, 3}};

  std::vector<Segment> target;
  target.reserve(source.size());
for (auto const &segment : source)
    {
      target.push_back(expected.apply(segment));
    }

  auto result = find_transform<Segment>(source, target, ToleranceConfig{});

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