#include "numerical.h"
#include <cmath>
#include <iostream>
#include "imnmath.hpp"
#include <fftw3.h>

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
                      double speedFactor, int cols, double (*dirfunc)(const double))
{
    double **vDistribution = imn<double>::matrix(cols, cols);
    double radius = 0;
    double x2;
    double y2;
    for (int i = 0; i < cols; ++i)
        for (int j = 0; j < cols; ++j)
        {
            y2 = i - spinCenterY;
            x2 = j - spinCenterX;
            radius = sqrt(pow(x2, 2) + pow(y2, 2));
            if (radius <= range)
                vDistribution[i][j] = calculateSpinV(range, radius, speedFactor, calculateAngle(x2, y2, radius), dirfunc);
        }
    return vDistribution;
}
void fillVDistribVortex(double spinCenterX, double spinCenterY, double range, double moveStep,
                        double speedFactor, int cols, double **vDistributionX, double **vDistributionY)
{
    double *zeta = new double[cols * cols];
    double *cos_psi = new double[cols * cols];
    double *psi = new double[cols * cols];
    fftw_plan plan_1, plan_2;
    double radius = 0;
    double x2, y2, gamma;
    for (int i = 0; i < cols; ++i)
        for (int j = 0; j < cols; ++j)
        {
            y2 = i - spinCenterY;
            x2 = j - spinCenterX;
            radius = sqrt(pow(x2, 2) + pow(y2, 2));
            gamma = 1.1;
            zeta[i * cols + j] = -(gamma * pow(radius / range, 2) - 1) * exp(-pow(radius / range, 2));
            // if (radius <= range)
            //     vDistribution[i][j] = calculateSpinV(range, radius, speedFactor, calculateAngle(x2, y2, radius), dirfunc);
        }
    plan_1 = fftw_plan_r2r_2d(cols, cols, zeta, cos_psi, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(plan_1);
    double dkx = 2. * M_PI / cols;
    double dky = 2. * M_PI / cols;
    double norma = cols * cols;
    double skx, sky, sk;
    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int l = j + i * cols;
            skx = dkx * i;
            sky = dky * j;
            sk = skx * skx + sky * sky;
            if (i == 0 && j == 0)
            {
                cos_psi[l] = 0.;
            }
            else
            {
                cos_psi[l] = cos_psi[l] / sk / norma;
            }
        }
    }
    plan_2 = fftw_plan_r2r_2d(cols, cols, cos_psi, psi, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(plan_2);

    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int l = i * cols + j;
            if ((j != 0) && (j != (cols - 1)))
                vDistributionX[i][j] = (psi[l + 1] - psi[l - 1]) / 2. / moveStep;
            if ((i != 0) && (i != (cols - 1)))
                vDistributionY[i][j] = -(psi[l + cols] - psi[l - cols]) / 2. / moveStep;
        }
    }
    fftw_destroy_plan(plan_1);
    fftw_destroy_plan(plan_2);
    fftw_cleanup();
    delete[] psi;
    delete[] zeta;
    delete[] cos_psi;
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
void licz_krok_fct(double dt, double dx, int n, int nx, int ny, double *vxn, double *vyn, double *un, double *un1)
{

    double q1, q2, q3, q4, q5, c;
    double vxp, vyp;
    double *flx, *gly, *fhx, *ghy;
    double *qlx, *qly;
    double qhx, qhy;

    double *ax, *ay;
    double *cx, *cy;
    double *qmin, *qmax;
    double *qtd;

    int i, j, l;

    flx = new double[n];
    gly = new double[n];
    fhx = new double[n];
    ghy = new double[n];

    qlx = new double[n]; //q^L_(i+1/2,j)
    qly = new double[n]; //q^L_(i,j+1/2)

    ax = new double[n]; //a_(i+1/2,j)
    ay = new double[n]; //a_(i,j+1/2)

    cx = new double[n]; //c_(i+1/2,j)
    cy = new double[n]; //c_(i,j+1/2)

    qtd = new double[n];  //q_td
    qmin = new double[n]; //q_min
    qmax = new double[n]; //q_max

    for (l = 0; l < n; l++)
    {
        flx[l] = 0.;
        gly[l] = 0.;
        fhx[l] = 0.;
        ghy[l] = 0.;
        qlx[l] = 0.;
        qly[l] = 0.;
        ax[l] = 0.;
        ay[l] = 0.;
        qtd[l] = 0.;
        cx[l] = 0.;
        cy[l] = 0.;
        qmin[l] = 0.;
        qmax[l] = 0.;
    }
    //lower order fluxes

    for (i = 0; i < (nx - 1); i++)
    {
        for (j = 0; j < (ny - 1); j++)
        {
            l = j + i * ny;

            vxp = (vxn[l] + vxn[l + ny]) / 2.;
            if (vxp > 0)
            {
                qlx[l] = un[l];
            }
            else
            {
                qlx[l] = un[l + ny];
            }

            vyp = (vyn[l] + vyn[l + 1]) / 2.;
            if (vyp > 0)
            {
                qly[l] = un[l];
            }
            else
            {
                qly[l] = un[l + 1];
            }

            flx[l] = qlx[l] * vxp * dx * dt;
            gly[l] = qly[l] * vyp * dx * dt;
        }
    }

    //higher order fluxes

    for (i = 1; i < (nx - 2); i++)
    {
        for (j = 1; j < (ny - 2); j++)
        {
            l = j + i * ny;

            qhx = 7. / 12. * (un[l + ny] + un[l]) - 1. / 12. * (un[l + 2 * ny] + un[l - ny]);
            qhy = 7. / 12. * (un[l + 1] + un[l]) - 1. / 12. * (un[l + 2] + un[l - 1]);

            vxp = (vxn[l] + vxn[l + ny]) / 2.;
            vyp = (vyn[l] + vyn[l + 1]) / 2.;

            fhx[l] = qhx * vxp * dx * dt;
            ghy[l] = qhy * vyp * dx * dt;
        }
    }

    //antidiffusive fluxes

    for (i = 1; i < (nx - 2); i++)
    {
        for (j = 1; j < (ny - 2); j++)
        {
            l = j + i * ny;
            ax[l] = fhx[l] - flx[l];
            ay[l] = ghy[l] - gly[l];
        }
    }

    //low-order transported-diffusedsolution

    for (i = 1; i < (nx - 1); i++)
    {
        for (j = 1; j < (ny - 1); j++)
        {
            l = j + i * ny;
            qtd[l] = un[l] - (flx[l] - flx[l - ny] + gly[l] - gly[l - 1]) / dx / dx;
        }
    }

    //qmin,qmax

    for (i = 1; i < (nx - 1); i++)
    {
        for (j = 1; j < (ny - 1); j++)
        {
            l = j + i * ny;

            q1 = fmax(un[l], qtd[l]);
            q2 = fmax(un[l - ny], qtd[l - ny]);
            q3 = fmax(un[l + ny], qtd[l + ny]);
            q4 = fmax(un[l - 1], qtd[l - 1]);
            q5 = fmax(un[l + 1], qtd[l + 1]);
            qmax[l] = fmax(q1, fmax(q2, fmax(q3, fmax(q4, q5))));

            q1 = fmin(un[l], qtd[l]);
            q2 = fmin(un[l - ny], qtd[l - ny]);
            q3 = fmin(un[l + ny], qtd[l + ny]);
            q4 = fmin(un[l - 1], qtd[l - 1]);
            q5 = fmin(un[l + 1], qtd[l + 1]);
            qmin[l] = fmin(q1, fmin(q2, fmin(q3, fmin(q4, q5))));
        }
    }

    //Cx,Cy

    for (i = 1; i < (nx - 2); i++)
    {
        for (j = 1; j < (ny - 2); j++)
        {
            l = j + i * ny;

            if (ax[l] > 0)
            {
                cx[l] = fmin(rplus(l + ny, ny, dx, ax, ay, qmax, qtd), rminus(l, ny, dx, ax, ay, qmin, qtd));
            }
            else
            {
                cx[l] = fmin(rplus(l, ny, dx, ax, ay, qmax, qtd), rminus(l + ny, ny, dx, ax, ay, qmin, qtd));
            }

            if (ay[l] > 0)
            {
                cy[l] = fmin(rplus(l + 1, ny, dx, ax, ay, qmax, qtd), rminus(l, ny, dx, ax, ay, qmin, qtd));
            }
            else
            {
                cy[l] = fmin(rplus(l, ny, dx, ax, ay, qmax, qtd), rminus(l + 1, ny, dx, ax, ay, qmin, qtd));
            }
        }
    }

    //modyfikacja: Ax,Ay
    for (i = 1; i < (nx - 2); i++)
    {
        for (j = 1; j < (ny - 2); j++)
        {
            l = j + i * ny;

            c = ax[l] * (qtd[l + ny] - qtd[l]);
            if (c < 0)
                ax[l] = 0.;

            c = ay[l] * (qtd[l + 1] - qtd[l]);
            if (c < 0)
                ay[l] = 0.;
        }
    }

    //Acx,Acy
    for (i = 1; i < (nx - 2); i++)
    {
        for (j = 1; j < (ny - 2); j++)
        {
            l = j + i * ny;

            if (cx[l] > 0 && cx[l] <= 1.)
                ax[l] = ax[l] * cx[l];
            if (cy[l] > 0 && cy[l] <= 1.)
                ay[l] = ay[l] * cy[l];
        }
    }

    //un1, prawa strona rowania rozniczkowego (ff)
    for (i = 1; i < (nx - 1); i++)
    {
        for (j = 1; j < (ny - 1); j++)
        {
            l = j + i * ny;
            un1[l] = qtd[l] - (ax[l] - ax[l - ny] + ay[l] - ay[l - 1]) / dx / dx;
        }
    }

    delete[] flx;
    delete[] gly;
    delete[] fhx;
    delete[] ghy;
    delete[] qlx;
    delete[] qly;

    delete[] ax;
    delete[] ay;
    delete[] cx;
    delete[] cy;
    delete[] qmin;
    delete[] qmax;
    delete[] qtd;

    return;
}

double rplus(int l, int ny, double dx, double *ax, double *ay, double *qmax, double *qtd)
{
    double pplus, qplus, rp;
    pplus = fmax(ax[l - ny], 0) - fmin(ax[l], 0) + fmax(ay[l - 1], 0) - fmin(ay[l], 0);
    qplus = (qmax[l] - qtd[l]) * dx * dx;
    if (pplus > 0)
    {
        rp = fmin(1., qplus / pplus);
    }
    else
    {
        rp = 0.;
    }
    return rp;
}

double rminus(int l, int ny, double dx, double *ax, double *ay, double *qmin, double *qtd)
{
    double pminus, qminus, rm;
    pminus = fmax(ax[l], 0) - fmin(ax[l - ny], 0) + fmax(ay[l], 0) - fmin(ay[l - 1], 0);
    qminus = (qtd[l] - qmin[l]) * dx * dx;
    if (pminus > 0)
    {
        rm = fmin(1., qminus / pminus);
    }
    else
    {
        rm = 0.;
    }
    return rm;
}
}