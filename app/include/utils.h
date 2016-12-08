#pragma once

namespace utils
{
void fill2Darray(int x1, int y1, int x2, int y2, double **array, double val);
void logArray(int x1, int y1, int x2, int y2, double **array);
constexpr int to1D(int cols, int x, int y)
{
    return y*cols + x;
}
}