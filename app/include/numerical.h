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
void fillVDistribVortex(double spinCenterX, double spinCenterY, double range, double moveStep,
                        double speedFactor, int cols, double **vDistributionX, double **vDistributionY);
CSR fillAlphaMatrix(double **velocityArrayX, double **velocityArrayY, int nrows, int ncols, double timeStep, double moveStep);
CSR fillBetaMatrix(double **velocityArrayX, double **velocityArrayY, int nrows, int ncols, double timeStep, double moveStep);
double calculateModulus(double timeStep, double moveStep, double velocity);
double rminus(int l, int ny, double dx, double *ax, double *ay, double *qmin, double *qtd);
double rplus(int l, int ny, double dx, double *ax, double *ay, double *qmax, double *qtd);
void licz_krok_fct(double dt, double dx, int n, int nx, int ny, double *vxn, double *vyn, double *un, double *un1);
}