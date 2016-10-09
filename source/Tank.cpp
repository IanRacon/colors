#include "Tank.h"
#include "utils.h"

Tank::Tank(int sizeX, int sizeY):sizeX(sizeX), sizeY(sizeY)
{

}
int Tank::getSizeX() const
{
    return sizeX;
}
void Tank::fillDensityMap(double value)
{
    utils::fill2Darray(0, 0, sizeX, sizeY, densityMap, value);
}
int Tank::getSizeY() const
{
    return sizeY;
}