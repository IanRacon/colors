#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include "imnmath.hpp"
//#include "easylogging++.h"
#include "INIReader.h"
#include "Configuration.h"
#include "utils.h"
#include "Tank.h"
#include "CSR.h"
#include "numerical.h"
#include "procedures.h"
#include <numeric>

//INITIALIZE_EASYLOGGINGPP

void mainRoutine()
{
    int cols = Configuration::getInstance().getInt("sizeX");
    int rows = Configuration::getInstance().getInt("sizeY");
    double time = Configuration::getInstance().getDouble("time");
    double timeStep = Configuration::getInstance().getDouble("timeStep");
    double moveStep = Configuration::getInstance().getDouble("moveStep");
    int colorPlacementX = Configuration::getInstance().getDouble("colorPlacementX");
    int colorPlacementY = Configuration::getInstance().getDouble("colorPlacementY");
    double densityValue = Configuration::getInstance().getDouble("densityValue");
    double mixerX = Configuration::getInstance().getDouble("mixerX");
    double mixerY = Configuration::getInstance().getDouble("mixerY");
    double mixerSize = Configuration::getInstance().getDouble("mixerSize");
    double mixerSpeed = Configuration::getInstance().getDouble("mixerSpeed");
    int frames = Configuration::getInstance().getInt("frames");
    bool printGraphs = Configuration::getInstance().getBool("printGraphs");

    double **velocityArrayX = numerical::fillVDistrib(mixerX, mixerY, mixerSize, mixerSpeed, rows, numerical::clockwiseX);
    double **velocityArrayY = numerical::fillVDistrib(mixerX, mixerY, mixerSize, mixerSpeed, rows, numerical::clockwiseY);
    if (printGraphs)
    {
        imnd::plot_2d_system("velocityX.png", velocityArrayX, rows, cols, 1, 1);
        imnd::plot_2d_system("velocityY.png", velocityArrayY, rows, cols, 1, 1);
    }

    // CSR betaMatrix = numerical::initialfillRoDistrib(velocityArrayX, velocityArrayY, rows, timeStep, moveStep);
    CSR betaMatrix = numerical::fillBetaMatrix(velocityArrayX, velocityArrayY, rows, timeStep, moveStep);
    CSR alphaMatrix(rows * rows, cols * cols, rows * cols);
    std::vector<double> initialDensityMatrix(cols * rows);
    initialDensityMatrix[colorPlacementY * rows + colorPlacementX] = densityValue;

    std::vector<double> rhsVector = procedures::product(betaMatrix, initialDensityMatrix);
    std::vector<double> result = initialDensityMatrix;
    std::vector<double> x0(cols * rows);
    int counter = 0;
    ofstream file;
    file.open("total_sum.dat");
    for (double i = 0; i < time; i += timeStep)
    {
        //std::cout << "time: " << i << std::endl;
        alphaMatrix = numerical::fillAlphaMatrix(velocityArrayX, velocityArrayY, rows, timeStep, moveStep);
        result = procedures::conjugateGradient(alphaMatrix, rhsVector, result);
        betaMatrix = numerical::fillBetaMatrix(velocityArrayX, velocityArrayY, rows, timeStep, moveStep);
        rhsVector = procedures::product(betaMatrix, result);
        if (counter++ % int(time / timeStep / frames) == 0)
        {
            file << std::accumulate(result.begin(), result.end(), 0.0) << std::endl;
            imnd::push_vector_data2D(result, rows, cols);
        }
    }
    file.close();
    imnd::write_data2D("density.dat", rows, cols, moveStep, moveStep);
    imnd::free_data2D(rows, cols);
}
int main(void)
{
    mainRoutine();
    return 0;
}
