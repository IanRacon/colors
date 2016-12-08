#include "numerical.h"
#include <cmath>
#include <iostream>
#include "imnmath.hpp"

namespace numerical
{
double clockwiseY(const double angle)
{
    return sin(angle);
}
double clockwiseX(const double angle)
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
                      double speedFactor, int containerSize, double (*dirfunc)(const double))
{
    double **vDistribution = imn<double>::matrix(containerSize, containerSize);
    double radius = 0;
    double x2;
    double y2;
    for (int i = 0; i < containerSize; ++i)
        for (int j = 0; j < containerSize; ++j)
        {
            y2 = i - spinCenterY;
            x2 = j - spinCenterX;
            radius = sqrt(pow(x2, 2) + pow(y2, 2));
            if (radius <= range)
                vDistribution[i][j] = calculateSpinV(range, radius, speedFactor, calculateAngle(x2, y2, radius), dirfunc);
        }
    return vDistribution;
}
CSR fillAlphaMatrix(double **velocityX, double **velocityY, int nx, int ny, double dt, double dx)
{
    CSR alphaMatrix(nx * nx, ny * ny, nx * ny);
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            double l = i * nx + j;
            double alpha1 = -(dt * velocityX[i][j]) / (4.0 * dx);
            double alpha2 = -(dt * velocityY[i][j]) / (4.0 * dx);
            double alpha3 = 1;
            double alpha4 = -alpha2;
            double alpha5 = -alpha1;
            // double alpha1 = 1;
            // double alpha2 = 2;
            // double alpha3 = 3;
            // double alpha4 = 4;
            // double alpha5 = 5;
            if (l - nx >= 0)
                alphaMatrix.setValue(l, l - nx, alpha1);
            if (l - 1 >= 0)
                alphaMatrix.setValue(l, l - 1, alpha2);

            alphaMatrix.setValue(l, l, alpha3);

            if (l + 1 < nx * nx)
                alphaMatrix.setValue(l, l + 1, alpha4);
            if (l + nx < nx * nx)
                alphaMatrix.setValue(l, l + nx, alpha5);
        }
    }
    alphaMatrix.setEndIndicator();
    return alphaMatrix;
}
CSR fillBetaMatrix(double **velocityX, double **velocityY, int nx, int ny, double dt, double dx)
{
    CSR betaMatrix(nx * nx, ny * ny, nx * ny);
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            double l = i * nx + j;
            double beta1 = (dt * velocityX[i][j]) / (4.0 * dx);
            double beta2 = (dt * velocityY[i][j]) / (4.0 * dx);
            double beta3 = 1;
            double beta4 = -beta2;
            double beta5 = -beta1;
            if (l - nx >= 0)
                betaMatrix.setValue(l, l - nx, beta1);
            if (l - 1 >= 0)
                betaMatrix.setValue(l, l - 1, beta2);

            betaMatrix.setValue(l, l, beta3);

            if (l + 1 < nx * nx)
                betaMatrix.setValue(l, l + 1, beta4);
            if (l + nx < nx * nx)
                betaMatrix.setValue(l, l + nx, beta5);
        }
    }
    betaMatrix.setEndIndicator();
    return betaMatrix;
}
}