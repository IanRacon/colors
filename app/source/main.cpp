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
#include <stdexcept>

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
    std::string dataName = Configuration::getInstance().getString("outputDataName");
    std::string graphsName = Configuration::getInstance().getString("outputGraphsName");

    // double **velocityArrayX = numerical::fillVDistrib(mixerX, mixerY, mixerSize, mixerSpeed, rows, numerical::clockwiseX);
    // double **velocityArrayY = numerical::fillVDistrib(mixerX, mixerY, mixerSize, mixerSpeed, rows, numerical::clockwiseY);
    double **velocityArrayX = imn<double>::matrix(cols, cols);
    double **velocityArrayY = imn<double>::matrix(cols, cols);
    numerical::fillVDistribVortex(mixerX, mixerY, mixerSize, moveStep, mixerSpeed, cols, velocityArrayX, velocityArrayY);
    if (printGraphs)
    {
        imnd::plot_2d_system("velocityX.png", velocityArrayX, rows, cols, 1, 1);
        imnd::plot_2d_system("velocityY.png", velocityArrayY, rows, cols, 1, 1);
    }

    CSR betaMatrix = numerical::fillBetaMatrix(velocityArrayX, velocityArrayY, rows, cols, timeStep, moveStep);
    CSR alphaMatrix(cols * cols, rows * rows, cols * rows);
    std::vector<double> initialDensityMatrix = procedures::getInitialDensity(colorPlacementX, colorPlacementY, rows, cols, moveStep, 8 * moveStep);

    std::vector<double> rhsVector = procedures::product(betaMatrix, initialDensityMatrix);
    std::vector<double> result = initialDensityMatrix;
    std::vector<double> x0(cols * rows);
    int counter = 0;
    int snapshotCounter = 0;
    for (double i = 0; i < time; i += timeStep)
    {
        alphaMatrix = numerical::fillAlphaMatrix(velocityArrayX, velocityArrayY, rows, cols, timeStep, moveStep);
        result = procedures::bicgstab(alphaMatrix, rhsVector, x0);

        //result = procedures::conjugateGradient(alphaMatrix, rhsVector, result);
        betaMatrix = numerical::fillBetaMatrix(velocityArrayX, velocityArrayY, rows, cols, timeStep, moveStep);

        rhsVector = procedures::product(betaMatrix, result);

        if (counter++ % int(time / timeStep / 10) == 0)
        {
            ++snapshotCounter;
            std::string filename = dataName + "_" + std::to_string(snapshotCounter) + ".dat";
            std::cout << "Zapisywanie danych do pliku: " << filename << std::endl;
            utils::saveData2D(filename, result, rows, cols);
        }
        if (counter % int(time / timeStep / frames) == 0)
        {
            std::cout << "Czas programu: " << i << "s, czas całkowity: " << time << "s" << std::endl;
            imnd::push_vector_data2D(result, rows, cols);
        }
    }
    imnd::write_data2D("density.dat", rows, cols, moveStep, moveStep);
    imnd::free_data2D(rows, cols);
    delete[] velocityArrayX;
    delete[] velocityArrayY;
}
int main(void)
{
    mainRoutine();
    return 0;
}
