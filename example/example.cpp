#include "geoma-cp.h"

#include <iostream>

int
main()
{
  std::vector<MyTypes::Segment> source{{{0, 0}, {3, 0}}, {{3, 0}, {3, 2}}};
  std::vector<MyTypes::Segment> target{{{5, 3}, {5, 6}}, {{5, 6}, {3, 6}}};

  auto result = find_transform<MyTypes::Segment, MyTypes::Transform>(
    source, target, geoma::ToleranceConfig{});

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