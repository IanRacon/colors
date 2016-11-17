#pragma once
#include "CSR.h"
namespace numerical
{
double clockwiseX(const double angle);
double clockwiseY(const double angle);
double counterClockwiseX(double angle);
double counterClockwiseY(double angle);
double calculateAngle(int x, int y, double radius);
double calculateSpinV(double range, double radius, double speedFactor, double angle, double (*dirfunc)(const double));
double **fillVDistrib(double spinCenterX, double spinCenterY, double range,
                      double speedFactor, int containerSize, double (*dirfunc)(const double));
CSR fillRoDistrib(double **velocityArrayX, double **velocityArrayY, int matrixSize, double timeStep, double moveStep);
CSR initialfillRoDistrib(double **velocityArrayX, double **velocityArrayY, int matrixSize, double timeStep, double moveStep);
double calculateModulus(double timeStep, double moveStep, double velocity);
}