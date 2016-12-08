#pragma once
#include "numerical.h"
#include <math.h>

namespace mixertraits
{
double middleFastMixer(int x, int y, double range, double speed, double (*dirfunc)(const double))
{
    double radius = sqrt(pow(x, 2) + pow(y, 2));
    double angle = numerical::calculateAngle(x, y, radius);
    if (radius <= range)
        return speed * sin(M_PI * radius / range) * dirfunc(angle);
    else
        return 0;
}
}