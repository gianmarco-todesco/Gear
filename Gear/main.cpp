#include "Viewer.h"
#include <QtGui/QApplication>
#include <qdebug.h>

#include "Curve.h"
#include <qmath.h>
#include <assert.h>

class EllipsePolarFunction : public PolarFunction {
public:
  double getRadius(double phi) const {
    return 1.0/(1.0 + 0.6*cos(phi));
  }
};

/*
class SquarePolarFunction : public PolarFunction {
  double m_edgeLength, m_cornerRadius;
public:
  SquarePolarFunction() : m_edgeLength(2), m_cornerRadius(0.5) {}
  double getRadius(double phi) const {

    return 1.0/(1.0 + 0.6*cos(phi));
  }
};
*/


void foo()
{
  EllipsePolarFunction f;
  Curve crv;
  crv.build(f);
  int m = 1000000;
  QVector2D oldpos;
  double length = 0.0;
  for(int i=0;i<m;i++)
  {
    double phi = 2*M_PI*i/(m-1);
    QVector2D posref = f.getPoint(phi);
    QVector2D normref = f.getNormal(phi);
    
    QVector2D pos = crv.getPointFromPhi(phi);

    assert((pos-posref).length()<1e-6);
    if(i>0) length += (pos-oldpos).length();
    oldpos = pos;
    double s = crv.getSfromPhi(phi);
    assert(fabs(s-length)<1.0e-3);

    double phi2 = crv.getPhifromS(s);
    assert(fabs(phi2-phi)<1.0e-7);

    QVector2D norm = crv.getNormalFromPhi(phi);
    assert((norm-normref).length()<1e-2);
  }

}


void foo2()
{
  Curve *crv = makeEllipse(0.5,0.7);
  crv->getSfromPhi(6.2046454908398410);
  
}


void foo3()
{
  Curve *crv1 = makeSquare(0.5,0.1);
  Curve *crv2 = makeConjugate(crv1, 1.0);
  delete crv1;
  delete crv2;


}


int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  foo3();
  // foo();
  Viewer w;
  w.show();
  return a.exec();
}
