#include "Tank.h"

Tank::Tank(int sizeX, int sizeY):sizeX(sizeX), sizeY(sizeY)
{

}
int Tank::getSizeX() const
{
    return sizeX;
}
int Tank::getSizeY() const
{
    return sizeY;
}