#pragma once
namespace function
{
double clockwiseX(double angle);
double clockwiseY(double angle);
double counterClockwiseX(double angle);
double counterClockwiseY(double angle);
double calculateAngle(int x, int y, double radius);
double calculateSpinV(double range, double radius, double speedFactor, double angle, double (*dirfunc)(const double));
double **fillVDistrib(double spinCenterX, double spinCenterY, double range,
                      double speedFactor, int size, double (*dirfunc)(const double));
}