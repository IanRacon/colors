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
CSR fillAlphaMatrix(double **velocityArrayX, double **velocityArrayY, int matrixSize, double timeStep, double moveStep)
{
    double scatterMatrixSize = matrixSize * matrixSize;
    CSR scatterMatrix(scatterMatrixSize, scatterMatrixSize, scatterMatrixSize);
    for (int i = 0; i < matrixSize; ++i)
    {
        for (int j = 0; j < matrixSize; ++j)
        {
            double l = i * matrixSize + j;
            double alpha1 = +calculateModulus(timeStep, moveStep, velocityArrayX[i][j]);
            double alpha2 = +calculateModulus(timeStep, moveStep, velocityArrayY[i][j]);
            double alpha3 = 1;
            double alpha4 = -calculateModulus(timeStep, moveStep, velocityArrayY[i][j]);
            double alpha5 = -calculateModulus(timeStep, moveStep, velocityArrayX[i][j]);
            if (l - matrixSize >= 0)
                scatterMatrix.setValue(l, l - matrixSize, alpha1);
            if (l - 1 >= 0)
                scatterMatrix.setValue(l, l - 1, alpha2);

            scatterMatrix.setValue(l, l, alpha3);

            if (l + 1 < scatterMatrixSize)
                scatterMatrix.setValue(l, l + 1, alpha4);
            if (l + matrixSize < scatterMatrixSize)
                scatterMatrix.setValue(l, l + matrixSize, alpha5);
        }
    }
    scatterMatrix.setEndIndicator();
    return scatterMatrix;
}
CSR fillBetaMatrix(double **velocityArrayX, double **velocityArrayY, int matrixSize, double timeStep, double moveStep)
{
    double scatterMatrixSize = matrixSize * matrixSize;
    CSR scatterMatrix(scatterMatrixSize, scatterMatrixSize, scatterMatrixSize);
    for (int i = 0; i < matrixSize; ++i)
    {
        for (int j = 0; j < matrixSize; ++j)
        {
            double l = i * matrixSize + j;
            double beta1 = -calculateModulus(timeStep, moveStep, velocityArrayX[i][j]);
            double beta2 = -calculateModulus(timeStep, moveStep, velocityArrayY[i][j]);
            double beta3 = 1;
            double beta4 = +calculateModulus(timeStep, moveStep, velocityArrayY[i][j]);
            double beta5 = +calculateModulus(timeStep, moveStep, velocityArrayX[i][j]);
            if (l - matrixSize >= 0)
                scatterMatrix.setValue(l, l - matrixSize, beta1);
            if (l - 1 >= 0)
                scatterMatrix.setValue(l, l - 1, beta2);

            scatterMatrix.setValue(l, l, beta3);

            if (l + 1 < scatterMatrixSize)
                scatterMatrix.setValue(l, l + 1, beta4);
            if (l + matrixSize < scatterMatrixSize)
                scatterMatrix.setValue(l, l + matrixSize, beta5);
        }
    }
    scatterMatrix.setEndIndicator();
    return scatterMatrix;
}
double calculateModulus(double timeStep, double moveStep, double velocity)
{
    return (timeStep * velocity) / (4.0 * moveStep);
}
}