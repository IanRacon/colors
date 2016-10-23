#include "function.h"
#include <cmath>

namespace function
{
double calculateAngle(int x, int y, double radius)
{
    if (y < 0)
        return M_PI * 2 - acos(x / radius);
    else if (y >= 0)
        return acos(x / radius);
}
double calculateSpinVX(double spinRange, double radius, double speedFactor, double angle)
{
    return speedFactor * sin(M_PI * radius / spinRange) * sin(angle);
}
double calculateSpinVY(double spinRange, double radius, double speedFactor, double angle)
{
    return speedFactor * sin(M_PI * radius / spinRange) * (-cos(angle));
}
double **fillVDistributionX(double spinCenterX, double spinCenterY, double range,
                            double speedFactor, const int size)
{
    double **vDistributionX = new double *[size];
    for (int i = 0; i < size; ++i)
        vDistributionX[i] = new double[size];
    double radius = 0;
    int x2;
    int y2;
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
        {
            x2 = i - spinCenterY;
            y2 = j - spinCenterX;
            radius = sqrt(pow(x2, 2) + pow(y2, 2));
            if (radius <= range)
                vDistributionX[i][j] = calculateSpinVX(range, radius, speedFactor, calculateAngle(x2, y2, radius));
        }
}
double **fillVDistributionYdouble(double spinCenterX, double spinCenterY, double range,
                                  double speedFactor, int size)
{
    double **vDistributionY = new double *[size];
    for (int i = 0; i < size; ++i)
        vDistributionY[i] = new double[size];
    double radius = 0;
    int x2;
    int y2;
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
        {
            x2 = i - spinCenterY;
            y2 = j - spinCenterX;
            radius = sqrt(pow(x2, 2) + pow(y2, 2));
            if (radius <= range)
                vDistributionY[i][j] = calculateSpinVY(range, radius, speedFactor, calculateAngle(x2, y2, radius));
        }
}
}