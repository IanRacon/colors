#pragma once
#include <map>

struct Range{
    unsigned min;
    unsigned max;
};
enum class Color{BLACK, PURPLE, RED, YELLOW};
const std::map<Color, Range> colorRanges = {{Color::BLACK, {0u, 1u}}, {Color::PURPLE, {1u, 2u}}, 
                                            {Color::RED, {2u, 3u}}, {Color::YELLOW, {3u, 4u}}};

class Dye{
public:
    Dye(Color color);
private:
    Color color;
    double **densityMap;
};