#pragma once

class Tank{
public:
    Tank(int sizeX, int sizeY);
    int getSizeX() const;
    int getSizeY() const;
private:
    const int sizeX;
    const int sizeY;
    double** densityMap;
};