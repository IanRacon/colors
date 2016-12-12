#pragma once
#include <string>
#include <vector>
#include <fstream>

namespace utils
{
void fill2Darray(int x1, int y1, int x2, int y2, double **array, double val);
void logArray(int x1, int y1, int x2, int y2, double **array);
constexpr int to1D(int cols, int x, int y)
{
    return y * cols + x;
}
void saveData2D(const std::string &filename, const std::vector<double> &data, int rows, int cols);
void saveData2D(const std::string &filename, double *data, int rows, int cols);
void saveData2DRGB(const std::string &filename, double *r, double *g, double *b, int rows, int cols);
}