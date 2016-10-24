#include "function.h"
#include <cmath>
#include <iostream>
#include "imnmath.hpp"

namespace function
{
double clockwiseX(double angle)
{
    return sin(angle);
}
double clockwiseY(double angle)
{
    return -cos(angle);
}
double counterClockwiseX(double angle)
{
    return -sin(angle);
}
double counterClockwiseY(double angle)
{
    return cos(angle);
}
double calculateAngle(int x, int y, double radius)
{
    if (radius == 0)
        return 0;
    if (y < 0)
        return M_PI * 2 - acos(x / radius);
    else if (y >= 0)
        return acos(x / radius);
}
double calculateSpinV(double spinRange, double radius, double speedFactor, double angle, double (*dirfunc)(const double))
{
    return speedFactor * sin(M_PI * radius / spinRange) * dirfunc(angle);
}
double **fillVDistrib(double spinCenterX, double spinCenterY, double range,
                      double speedFactor, int size, double (*dirfunc)(const double))
{
    double **vDistribution = imn<double>::matrix(size, size);
    double radius = 0;
    double x2;
    double y2;
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
        {
            y2 = i - spinCenterY;
            x2 = j - spinCenterX;
            radius = sqrt(pow(x2, 2) + pow(y2, 2));
            if (radius <= range)
            {
                vDistribution[i][j] = calculateSpinV(range, radius, speedFactor, calculateAngle(x2, y2, radius), dirfunc);
                // std::cout << "vDistribution[" << i << "][" << j << "]"
                //           << " = " << vDistribution[i][j] << std::endl;
            }
        }
    return vDistribution;
}
}