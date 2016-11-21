#include "Mixer.h"
#include "Configuration.h"
#include "mixertraits.h"
#include "numerical.h"

Mixer::Mixer()
{
    size = Configuration::getInstance().getDouble("mixerSize");
    speed = Configuration::getInstance().getDouble("mixerSpeed");
}
double Mixer::velocityX(int mixerX, int mixerY) const
{
    return mixertraits::middleFastMixer(mixerX, mixerY, size, speed, numerical::clockwiseX);
}
double Mixer::velocityY(int mixerX, int mixerY) const
{
    return mixertraits::middleFastMixer(mixerX, mixerY, size, speed, numerical::clockwiseY);
}