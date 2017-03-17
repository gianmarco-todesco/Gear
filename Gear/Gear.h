#ifndef GEAR_H
#define GEAR_H

#include "Mesh.h"
#include "Curve.h"

class Gear : public Mesh
{
protected:
  Curve *m_curve;
  double m_rotation;
  QVector3D m_position;
public:

  Gear(Curve *curve);
  virtual ~Gear();

  double getRotation() const { return m_rotation; }
  void setRotation(double angle) { m_rotation = angle; }

  QVector3D getPosition() const { return m_position; }
  void setPosition(const QVector3D &pos) { m_position = pos; }

  void build();
  Curve *getCurve() { return m_curve; }
};


Gear *makeConjugate(Gear *g, double distance);


class GearLink {
  Gear *m_driver, *m_driven;
public:
  GearLink(Gear *driver, Gear *driven) : m_driver(driver), m_driven(driven) {}

  void update();
};



#endif

