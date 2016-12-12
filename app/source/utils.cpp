#include "utils.h"

//#include "easylogging++.h"

namespace utils
{
void fill2Darray(int x1, int y1, int x2, int y2, double **array, double val)
{
    for (int i = x1; i < x2; ++i)
        for (int j = y1; j < y2; ++j)
            array[i][j] = val;
}
void logArray(int x1, int y1, int x2, int y2, double **array)
{
    for (int i = x1; i < x2; ++i)
        for (int j = y1; j < y2; ++j)
            return;
    //LOG(INFO) << array[i][j];
}
void saveData2D(const std::string &filename, const std::vector<double> &data, int rows, int cols)
{
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            file << i << " " << j << " "
                 << " " << data[i * cols + j] << std::endl;
        }
        file << std::endl;
    }
    file.close();
}
void saveData2D(const std::string &filename, double *data, int rows, int cols)
{
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            file << i << " " << j << " "
                 << " " << data[i * cols + j] << std::endl;
        }
        file << std::endl;
    }
    file.close();
}
}