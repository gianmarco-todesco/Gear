#ifndef GEAR_H
#define GEAR_H

#include "Mesh.h"
#include "Curve.h"

class Gear : public Mesh
{
protected:
  Curve *m_curve;
public:
  double m_rotation;
  QVector3D m_position;

  AbstractGear(AbstractCurve *curve);
  virtual ~AbstractGear();

  void build();
  AbstractCurve *getCurve() { return m_curve; }
};


AbstractGear *makeConjugate(AbstractGear *g, double distance);


class GearLink {
  AbstractGear *m_driver, *m_driven;
public:
  GearLink(AbstractGear *driver, AbstractGear *driven) : m_driver(driver), m_driven(driven) {}

  void update();
};


#endif

