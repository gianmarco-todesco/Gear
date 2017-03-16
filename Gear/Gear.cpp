#include "Gear.h"
#include "Curve.h"
#include <qmath.h>

Gear::Gear(Curve *curve)
  : m_curve(curve)
  , m_rotation(0)
{
}

Gear::~Gear()
{
  delete m_curve;
}

void Gear::build()
{
  int m = 20;
  int toothCount = 64;
  int n = toothCount * m;
  struct ShapePoint { QVector3D p0,p1,nrm0,nrm1; }; 
  QVector<ShapePoint> shape(n);
  //Curves::Ellipse crv(1.0,0.7);
  //crv.build(1000);
  Curve &crv = *m_curve;

  for(int i=0;i<n;i++) 
  {
    double phi = crv.getPhifromS(crv.getLength()*i/n);
    QVector2D pos1 = crv.getPointFromPhi(phi);
    QVector2D nrm1 = crv.getNormalFromPhi(phi);
    QVector2D nrm0 = pos1.normalized();

    shape[i].p0 = nrm0 * 0.1;
    shape[i].nrm0 = -nrm0;

    double t = (double)(i%m)/(double)(m/2);
    if(t>=1.0) t=0.0;
    else t = 4*t*(1-t);
    
    shape[i].p1 = pos1 + nrm1*((t-0.5)*0.05);
  }
  for(int i=0;i<n;i++) 
  {
    shape[i].nrm1 = QVector3D::crossProduct(shape[(i+1)%n].p1 - shape[(i+n-1)%n].p1, QVector3D(0,0,1)).normalized();
  }

  QVector3D dz(0,0,0.05);
  for(int i=0;i<n;i++) 
  {
    addVertex(shape[i].p0+dz, QVector3D(0,0,1));
    addVertex(shape[i].p1+dz, QVector3D(0,0,1));

    addVertex(shape[i].p1+dz, shape[i].nrm1);
    addVertex(shape[i].p1-dz, shape[i].nrm1);

    addVertex(shape[i].p1-dz, QVector3D(0,0,-1));
    addVertex(shape[i].p0-dz, QVector3D(0,0,-1));

    addVertex(shape[i].p0-dz, shape[i].nrm0);
    addVertex(shape[i].p0+dz, shape[i].nrm0);
  }

  for(int i=0;i<n;i++) 
  {
    int k0 = i*8, k1 = ((i+1)%n)*8;
    // addFace(k0,k0+1,k1+1,k1);
    addFace(k0+2,k0+3,k1+3,k1+2);
    // addFace(k0+4,k0+5,k1+5,k1+4);
    addFace(k0+6,k0+7,k1+7,k1+6);
  }

  for(int q=0; q<toothCount; q++) 
  {
    int ia = q * m, ib = ia + m/2;
    int ja=ia; 
    int jb=ib;
    while(ja+1<jb-1)
    {
      addFace(ja*8+1, (ja+1)*8+1, (jb-1)*8+1, jb*8+1);
      addFace(ja*8+4, jb*8+4, (jb-1)*8+4, (ja+1)*8+4);
      ja++;
      jb--;
    }
    addFace(ia*8+1, ib*8+1, ib*8, ia*8);
    addFace(ia*8+5, ib*8+5, ib*8+4, ia*8+4);
    while(ib+1<=ia+m)
    {
      int ib1 = (ib+1)%n;
      addFace(ib*8+1, ib1*8+1, ib1*8, ib*8);
      addFace(ib*8+5, ib1*8+5, ib1*8+4, ib*8+4);
      ib++;
    }
  }

}


void GearLink::update()
{
  double s = m_driver->getCurve()->getSfromPhi(M_PI * m_driver->getRotation() / 180.0);
  s = m_driven->getCurve()->getLength()*0.5 - s;
  m_driven->setRotation(180 + 180 * m_driven->getCurve()->getPhifromS(s) / M_PI);
}



/*
Gear1::Gear1()
  : AbstractGear(makeEllipse(256,0.5,0.7))
{
  build();
  QVector2D qqq = m_curve->getPointFromT(0) - m_curve->getPointFromT(0.5);
  m_distance = sqrt(qqq.x()*qqq.x()+qqq.y()*qqq.y());
}
*/

double getNorm(const QVector2D &p) { return sqrt(p.x()*p.x()+p.y()*p.y()); }

/*
AbstractGear *makeConjugate(AbstractGear *g, double distance)
{
  AbstractCurve *srcCrv = g->getCurve();
  QVector<QVector2D> pts;
  double oldr = distance - getNorm(srcCrv->getPoint(0));
  pts.append(QVector2D(oldr, 0.0));
  QVector2D oldSrcP = srcCrv->getPoint(0);
  double phi = 0.0;
  for(int i=0;i<srcCrv->getPointCount();i++)
  {
    QVector2D srcP = srcCrv->getPoint(i);
    double r = distance - getNorm(srcP);
    double d = getNorm(srcP - oldSrcP);
    double dphi = acos((r*r+oldr*oldr-d*d)/(2*oldr*r));
    phi += dphi;
    pts.append(r*QVector2D(cos(phi), sin(phi)));
    oldr = r;
    oldSrcP = srcP;
  }
  AbstractCurve *crv = new AbstractCurve();
  crv->build(pts);
  return new AbstractGear(crv);
}
*/


