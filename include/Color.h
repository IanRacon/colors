#pragma once
#include <map>

enum class Color{BLACK, PURPLE, RED, YELLOW, GREEN};

class Dye{
public:
    Dye(Color color);
private:
    Color color;
    double **densityMap;
};