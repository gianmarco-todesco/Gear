#include "Curve.h"
#include <qmath.h>
#include <QDebug>
#include <assert.h>

const int TickCount = 20000;

void normalizePeriodicValue(double &x, int &q, double period)
{
  double oldx = x;
  q = 0;
  if(x>=period) 
  { 
    q = (int)floor(x/period); 
    x -= q*period;
  }
  else if(x<0.0)
  {
    q = - (int)ceil(-x/period); 
    x -= q*period;
  }
  assert(0.0<=x && x<period);
  assert(fabs(x + period*q - oldx)<1.0e-6);  
}

QVector2D PolarFunction::getPoint(double phi) const
{
  double r = getRadius(phi);
  return r*QVector2D(cos(phi), sin(phi));
}

QVector2D PolarFunction::getNormal(double phi) const
{
  double h = 0.0001;
  QVector2D v = (getPoint(phi+h) - getPoint(phi-h)).normalized();
  return QVector2D(v.y(),-v.x());
}

//=============================================================================


Curve::Curve()
  : m_isOpen(true)
  , m_length(0)
{
}

Curve::~Curve()
{
}

void Curve::build(const PolarFunction &f)
{
  m_pts.clear();
  m_length = 0;  

  int n = TickCount;
  m_pts.resize(n);
  
  m_pts[0].s = 0;
  m_pts[0].phi = 0.0;
  m_pts[0].radius = f.getRadius(0.0);
  m_pts[0].pos = QVector2D(m_pts[0].radius, 0.0);
  double totPhi = 2*M_PI;

  for(int i=1; i<n; i++)
  {
    Point &pt = m_pts[i];
    double phi = pt.phi = totPhi * i/n;
    pt.radius = f.getRadius(phi);
    pt.pos = pt.radius * QVector2D(cos(phi),sin(phi));
    m_length += (pt.pos - m_pts[i-1].pos).length();
    pt.s = m_length;
  }
  m_length += (m_pts.back().pos - m_pts.front().pos).length();
  
  // compute normals
  for(int i=0; i<n; i++)
    m_pts[i].norm = f.getNormal(m_pts[i].phi);

  postBuild();
}

void Curve::build(const QVector<QVector2D> &pts)
{
  assert(pts[0].x()>0.0 && fabs(pts[0].y())<1.0e-7);
  m_pts.clear();
  m_length = 0;  

  int n = pts.count();
  m_pts.resize(n);
  
  m_pts[0].s = 0;
  m_pts[0].phi = 0.0;
  m_pts[0].pos = pts[0];
  m_pts[0].radius = pts[0].x();

  for(int i=1; i<n; i++)
  {
    Point &pt = m_pts[i];
    pt.pos = pts[i];
    double phi = atan2(pt.pos.y(), pt.pos.x());
    if(phi<0.0)phi+=2*M_PI;
    pt.phi = phi;
    pt.radius = pt.pos.length();
    m_length += (pt.pos - m_pts[i-1].pos).length();
    pt.s = m_length;
  }
  m_length += (m_pts.back().pos - m_pts.front().pos).length();
  
  // compute normals
  for(int i=0; i<n; i++)
  {
    QVector2D v = m_pts[(i+1)%n].pos - m_pts[(i+n-1)%n].pos;
    m_pts[i].norm = QVector2D(v.y(),-v.x()).normalized();
  }

  postBuild();
}


void Curve::postBuild()
{
  int n = m_pts.count();
  m_sToIndex.resize(n);
  int j = 0;
  for(int i=0;i<n;i++)
  {
    double s = m_length * i/n;
    while(j+1<n && m_pts[j+1].s<=s) j++;
    m_sToIndex[i] = j; //  = max j : m_pts[j+1].s <= s
  }

  check1();
}

void Curve::check1()
{
  int m = 10000000;
  int n = m_sToIndex.count();
  for(int i=0;i<m;i++)
  {
    double s = m_length * i / m;
    assert(s<m_length);
    int j = (int)floor(s*n/m_length);
    assert(0<=j && j<n);
    int a = m_sToIndex[j];
    assert(0<=a && a<n);
    double s0 = m_pts[a].s;
    assert(s0<=s);
    double s1 = m_length;
    if(j+1<n)
    {
      int b = m_sToIndex[j+1];
      assert(0<=b && b<n);
      if(b+1<n)
        s1 = m_pts[b+1].s;
    }
    assert(s<=s1);
  }

  n = m_pts.count();
  for(int i=0;i<n;i++)
  {
    double phi = 2*M_PI * (i + 0.5) / n;
    assert(m_pts[i].phi < phi && phi < (i+1<n ? m_pts[i+1].phi : 2*M_PI));
    int j = (int)floor(phi*n/(2*M_PI));
    assert(j==i);
  }

  qDebug() << "Check1 ok";
}


int Curve::getIndexFromPhi(double phi, double *off) const
{
  const double roundAngle = 2*M_PI;
  assert(0<=phi && phi< roundAngle);
  int n = m_pts.count();
  int j = (int)floor(phi*n/roundAngle);
  assert(0<=j && j<n);
  double phi0 = m_pts[j].phi;
  if(phi0>phi && j>0) {j--; phi0 = m_pts[j].phi;}
  if(j+1<n && m_pts[j+1].phi<=phi) {j++; phi0 = m_pts[j].phi;}
  double phi1 = j+1<n ? m_pts[j+1].phi : roundAngle;
  assert(phi0<=phi && phi<=phi1);
  if(off) *off = (phi-phi0)/(phi1-phi0);
  return j;
}

QVector2D Curve::getPointFromPhi(double phi) const
{
  int q;
  normalizePeriodicValue(phi, q, 2*M_PI);
  double t;
  int j = getIndexFromPhi(phi, &t);
  QVector2D p0 = m_pts[j].pos;
  QVector2D p1 = j+1<m_pts.count() ? m_pts[j+1].pos : m_pts[0].pos;
  return p0*(1-t) + p1*t;
}

QVector2D Curve::getNormalFromPhi(double phi) const
{
  int q;
  normalizePeriodicValue(phi, q, 2*M_PI);
  double t;
  int j = getIndexFromPhi(phi, &t);
  QVector2D p0 = m_pts[j].norm;
  QVector2D p1 = j+1<m_pts.count() ? m_pts[j+1].norm : m_pts[0].norm;
  return (p0*(1-t) + p1*t).normalized();
}


double Curve::getRadiusFromPhi(double phi) const
{
  int q;
  normalizePeriodicValue(phi, q, 2*M_PI);
  double t;
  int j = getIndexFromPhi(phi, &t);
  double r0 = m_pts[j].radius;
  double r1 = j+1<m_pts.count() ? m_pts[j+1].radius : m_pts[0].radius;
  return r0*(1-t) + r1*t;
}



double Curve::getSfromPhi(double phi) const
{
  int q;
  normalizePeriodicValue(phi, q, 2*M_PI);
  double t;
  int j = getIndexFromPhi(phi, &t);
  double s0 = m_pts[j].s;
  double s1 = j+1<m_pts.count() ? m_pts[j+1].s : m_length;
  return s0*(1-t) + s1*t + q * m_length;
}


double Curve::getPhifromS(double s) const
{
  int q;
  normalizePeriodicValue(s, q, m_length);
  const double roundAngle = 2*M_PI;
  double s0,s1,phi0,phi1;
  if(s>=m_pts.back().s)
  {
    s0 = m_pts.back().s;
    s1 = m_length;
    phi0 = m_pts.back().phi;
    phi1 = roundAngle;
  }
  else
  {
    int n = m_pts.count(); 
    int m = m_sToIndex.count();
    int j = (int)floor(m*s/m_length);
    assert(0<=j && j<m);
    int a = m_sToIndex[j];
    assert(0<=a && a<n);
    int b = m_pts.count()-1;
    if(j+1<m)
    {
      b = m_sToIndex[j+1];
      if(b+1<n) b++;
    }

    assert(m_pts[a].s<=s && s<m_pts[b].s);
    while(b-a>1)
    {
      int c=(a+b)/2;
      if(m_pts[c].s<=s)a=c; else b=c;
    }
    s0 = m_pts[a].s;
    s1 = m_pts[b].s;
    phi0 = m_pts[a].phi;
    phi1 = m_pts[b].phi;
  }
  assert(s0<=s && s<s1);
  double t = (s-s0)/(s1-s0);
  return phi0 * (1-t) + phi1 * t + q * roundAngle;
}




/*
void build(const QVector<QVector2D> &pts);

  virtual QVector2D getPointFromT(double t) const;
  virtual QVector2D getNormalFromT(double t) const;

  double getLength() const { return m_length; }
  double getTfromS(double s) const;

  double getSfromPhi(double phi) const;
  double getPhifromS(double s) const;

  int getPointCount() const { return m_pts.count(); }
  QVector2D getPoint(int index) const { return m_pts[index].pos; }

  */

//------------------------------------------------


AbstractCurve::AbstractCurve()
  : m_length(0)
  , m_phi0(0)
{
}

AbstractCurve::~AbstractCurve()
{
}

double normalizeAngle(double phi)
{
  const double roundAngle = 2*M_PI;
  if(phi<0.0)
  {
    while(phi<0.0) phi += roundAngle;
  }
  else if(phi>=roundAngle)
  {
    while(phi>=0.0) phi -= roundAngle;
  }
  return phi;
}


double normalizeAngle(double phi, double oldPhi)
{
  const double roundAngle = 2*M_PI;
  phi = normalizeAngle(phi);
  if(phi>oldPhi) {if(phi-oldPhi>M_PI) phi-=roundAngle; }
  else {if(oldPhi-phi>M_PI) phi += roundAngle; }
  return phi;
}

void AbstractCurve::build(const CurveFunction &f, int n)
{
  m_pts.clear();
  m_length = 0;
  
  Point p;
  p.t = p.s = 0.0;
  p.pos = f.getPoint(0.0);
  m_pts.append(p);
  for(int i=1; i<n; i++)
  {
    p.t = (double)i/(double)(n-1);
    p.pos = f.getPoint(p.t);
    QVector2D d = p.pos - m_pts.back().pos;
    m_length += sqrt(d.x()*d.x()+d.y()*d.y());
    p.s = m_length;
    m_pts.append(p);
  }
  m_phi0 = atan2(m_pts[0].pos.y(),m_pts[0].pos.x());
  m_pts[0].phi = 0.0;
  for(int i=1; i<n; i++)
  {
    QVector2D pos = m_pts[i].pos;
    m_pts[i].phi = normalizeAngle(atan2(pos.y(),pos.x()) - m_phi0, m_pts[i-1].phi);
  }
}

void AbstractCurve::build(const QVector<QVector2D> &pts)
{
  m_pts.clear();
  m_length = 0;
  Point p;
  p.t = p.s = 0.0;
  p.pos = pts[0];
  m_pts.append(p);
  for(int i=1; i<pts.count(); i++)
  {
    p.t = (double)i/(double)(pts.count()-1);
    p.pos = pts[i];
    QVector2D d = p.pos - m_pts.back().pos;
    m_length += sqrt(d.x()*d.x()+d.y()*d.y());
    p.s = m_length;
    m_pts.append(p);
  }
  m_phi0 = atan2( m_pts[0].pos.y(), m_pts[0].pos.x());
  m_pts[0].phi = 0.0;
  for(int i=1; i<pts.count(); i++)
  {
    QVector2D pos = m_pts[i].pos;
    m_pts[i].phi = normalizeAngle(atan2(pos.y(),pos.x()) - m_phi0, m_pts[i-1].phi);
  }
}


QVector2D AbstractCurve::getPointFromT(double t) const
{
  if(t<=0.0) return m_pts[0].pos;
  else if(t>=1.0) return m_pts.back().pos;
  int a=0,b=m_pts.count()-1;
  while(b-a>1)
  {
    int c = (a+b)/2;
    if(m_pts[c].t<=t) a=c;
    else b=c;
  }
  double tt = (t - m_pts[a].t) / (m_pts[b].t-m_pts[a].t);
  return m_pts[a].pos * (1-tt) + m_pts[b].pos * tt;
}


QVector2D AbstractCurve::getNormalFromT(double t) const
{
  double h = 1.0e-4;
  QVector2D v = getPointFromT(t+h)-getPointFromT(t-h);
  return QVector2D(v.y(),-v.x()).normalized();
}

double AbstractCurve::getTfromS(double s) const
{
  if(s<=0.0) return 0.0;
  else if(s>=m_length) return 1.0;
  int a=0,b=m_pts.count()-1;
  while(b-a>1)
  {
    int c = (a+b)/2;
    if(m_pts[c].s<=s) a=c;
    else b=c;
  }
  double t = m_pts[a].t + (s-m_pts[a].s)*(m_pts[b].t-m_pts[a].t)/(m_pts[b].s-m_pts[a].s);
  return t;
}

double AbstractCurve::getSfromPhi(double phi) const
{
  double originalPhi = phi;
  double roundangle = 2*M_PI;
  int q = 0;
  
  if(phi>roundangle) { q = (int)floor(phi/roundangle); phi -= roundangle*q; }
  else if(phi<0) { q = -1-(int)ceil(-phi/roundangle); phi -= roundangle*q; }
  if(phi>=roundangle) { phi-=roundangle; q--; }
  assert(0<=phi && phi<=roundangle);
  //assert(fabs( originalPhi - (phi + q*roundangle)) < 1.0e-6);
  
  phi += m_phi0;
  double s = 0.0;
  if(phi<=m_pts.front().phi) s = m_pts.front().s;
  else if(phi>=m_pts.back().phi) s = m_pts.back().s;
  else
  {
    int a=0, b=m_pts.count()-1;
    while(b-a>1)
    {
      int c = (a+b)/2;
      if(m_pts[c].phi<=phi) a=c; else b=c;
    }
    s = m_pts[a].s + (m_pts[b].s-m_pts[a].s)*(phi-m_pts[a].phi)/ (m_pts[b].phi-m_pts[a].phi);
  }
  return s + q * m_length;
}

double AbstractCurve::getPhifromS(double s) const
{
  double originalS = s;
  int q = 0;
  
  if(s>m_length) { q = (int)floor(s/m_length); s -= m_length*q; }
  else if(s<0) { q = -(int)ceil(-s/m_length); s -= m_length*q; }
  assert(0<=s && s<m_length);
  //assert(fabs(originalS - ( s + q*m_length) ) < 1.0e-6) ;
  
  double phi = 0.0;
  if(s<=m_pts.front().s) phi = m_pts.front().phi;
  else if(s>=m_pts.back().s) s = m_pts.back().phi;
  else
  {
    int a=0, b=m_pts.count()-1;
    while(b-a>1)
    {
      int c = (a+b)/2;
      if(m_pts[c].s<=s) a=c; else b=c;
    }
    s = m_pts[a].phi + (m_pts[b].phi-m_pts[a].phi)*(s-m_pts[a].s)/ (m_pts[b].s-m_pts[a].s);
  }
  return (s + q * (2*M_PI)) + m_phi0;
}


//=============================================================================


class EllipseFunction : public CurveFunction {
  double m_r, m_e;
public:
  EllipseFunction(double r, double e) : m_r(r), m_e(e) {}
  QVector2D getPoint(double t) const {
    double phi = t*M_PI*2;
    double r = m_r/(1.0 + m_e*cos(phi));
    return r*QVector2D(cos(phi), sin(phi));
  }

};

/*

AbstractCurve *makeEllipse(int n, double r, double e)
{  
  AbstractCurve *curve = new AbstractCurve();
  curve->build(EllipseFunction(r,e), n);
  return curve;
}

*/

/*

class SquareFunction : public CurveFunction {
  double m_edgeLength, m_cornerRadius;
public:
  SquareFunction(double edgeLength, double cornerRadius) : m_edgeLength(edgeLength), m_cornerRadius(cornerRadius) {}
  QVector2D getPoint(double t) const {
    if(t<0.0)t=0.0; else if(t>1.0)t=1.0;
    t = 1-t;
    t *= 4.0;
    int i = (int)floor(t);
    t -= i;

    double s0 = m_cornerRadius * M_PI*0.5, s1 = s0 + m_edgeLength - 2*m_cornerRadius;
    t *= s1;

    double q = -m_edgeLength*0.5 + m_cornerRadius;
    QVector2D p;
    if(t<s0) { double phi = M_PI*0.5*t/s0; p = QVector2D(q - m_cornerRadius*cos(phi), -q +  m_cornerRadius*sin(phi)); }
    else { p = QVector2D(q + t-s0, m_edgeLength*0.5); }

    if(i==1) p = QVector2D(p.y(),-p.x());  
    else if(i==2) p = -p;
    else if(i==3) p = QVector2D(-p.y(),p.x()); 
    return p;
  }

};


AbstractCurve *makeSquare(int n, double edgeLength, double cornerRadius)
{
  AbstractCurve *curve = new AbstractCurve();
  curve->build(SquareFunction(edgeLength,cornerRadius), n);
  return curve;
}


*/

namespace {


class EllipsePolarFunction : public PolarFunction {
  double m_r,m_e;
public:
  EllipsePolarFunction(double r, double e) : m_r(r), m_e(e) {}
  double getRadius(double phi) const {
    return m_r/(1.0 + m_e*cos(phi));
  }
};

}


Curve *makeEllipse(double r, double e)
{  
  Curve *curve = new Curve();
  curve->build(EllipsePolarFunction(r,e));
  return curve;
}

class SquarePolarFunction : public PolarFunction {
  double m_radius, m_cornerRadius;
  double m_a,m_b;
public:
  SquarePolarFunction(double radius, double cornerRadius) 
    : m_radius(radius)
    , m_cornerRadius(cornerRadius)
  {
    m_a = m_radius - m_cornerRadius;
    m_b = sqrt(m_a*m_a+m_radius*m_radius);
  }
  double getRadius(double phi) const {
    double cs = cos(phi), sn = sin(phi);
    double absCs = fabs(cs), absSn = fabs(sn); 
    double r = m_radius/qMax(absCs, absSn);
    if(r > m_b)
    {
      QVector2D q(cs>0 ? m_a : -m_a, sn>0 ? m_a : -m_a);
      double cu = cs*q.x() + sn*q.y(), cc = q.lengthSquared();
      double dsc = cu*cu-cc+m_cornerRadius*m_cornerRadius;
      assert(dsc>0.0);
      r = cu+sqrt(dsc);
    }
    return r;
  }
};


Curve *makeSquare(double radius, double cornerRadius)
{
  Curve *curve = new Curve();
  curve->build(SquarePolarFunction(radius,cornerRadius));
  return curve;
}


struct MyPoint { double s, phi, r; };
struct MyPointPair { MyPoint src, dst; };

Curve *makeConjugate(const Curve *srcCurve, double dist)
{

  QVector<MyPointPair> zip;
  MyPointPair oldg, curg;


  curg.src.phi = M_PI;
  curg.src.s = srcCurve->getSfromPhi(M_PI);
  curg.src.r = srcCurve->getRadiusFromPhi(M_PI);

  curg.dst.s = curg.dst.phi = 0.0;
  curg.dst.r = dist - curg.src.r;

  oldg = curg;
  zip.append(curg);

  const double ds = 0.001;
  const double ds2 = ds*ds;
  for(;;)
  {
    curg.src.s = oldg.src.s - ds;
    curg.src.phi = srcCurve->getPhifromS(curg.src.s);
    curg.src.r = srcCurve->getRadiusFromPhi(curg.src.phi);

    curg.dst.s = oldg.dst.s + ds;
    curg.dst.r = dist - curg.src.r;
    double oldr = oldg.dst.r;
    double curr = curg.dst.r;
    double dphi = acos((curr*curr+oldr*oldr-ds2)/(2*oldr*curr)); 
    assert(dphi);
    curg.dst.phi = oldg.dst.phi + dphi; 
    zip.push_back(curg);
    oldg = curg;
    if(curg.dst.phi>2*M_PI) break;
    if((zip.count()%1000)==0)
    {
      qDebug() << curg.dst.phi;
    }
  }

  QVector<QVector2D> pts;
  int n = TickCount;
  pts.resize(n);
  int j = 0;
  for(int i=0;i<n;i++)
  {
    double phi = 2*M_PI*i/n;
    while(zip[j+1].dst.phi<=phi) j++;
    assert(zip[j].dst.phi<=phi && phi<zip[j+1].dst.phi);
    double t = (phi - zip[j].dst.phi)/(zip[j+1].dst.phi - zip[j].dst.phi);
    double r = (1-t)*zip[j].dst.r + t*zip[j+1].dst.r;
    pts[i] = r*QVector2D(cos(phi), sin(phi));
  }

  Curve *curve = new Curve();
  curve->build(pts);
  return curve;
}
