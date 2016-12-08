#pragma once

enum class Rotation
{
    CLOCKWISE,
    COUNTERCLOCKWISE
};
class Mixer
{
  public:
    Mixer();
    Mixer(double speed, double size, Rotation direction);
    double velocityX(int mixerX, int mixerY) const;
    double velocityY(int mixerX, int mixerY) const;

  private:
    double size;
    double speed;
};