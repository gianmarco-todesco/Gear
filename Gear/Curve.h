#ifndef CURVE_H
#define CURVE_H

#include <QVector>
#include <QPair>
#include <QVector2D>


class CurveFunction {
public:
  virtual QVector2D getPoint(double phi) const = 0;
};

class PolarFunction {
public:
  virtual double getRadius(double phi) const = 0;
  QVector2D getPoint(double phi) const;
  QVector2D getNormal(double phi) const;
};


class Curve {
  struct Point {
    double s;
    double phi;
    double radius;
    QVector2D pos, norm;
  };
  bool m_isOpen;
  QVector<Point> m_pts;
  QVector<int> m_sToIndex;
  
  double m_length;  

public:

  Curve();
  virtual ~Curve();

  bool isOpen() const { return m_isOpen; }

  void build(const PolarFunction &f);
  void build(const QVector<QVector2D> &pts);

  QVector2D getPointFromPhi(double phi) const;
  QVector2D getNormalFromPhi(double phi) const;

  double getRadiusFromPhi(double phi) const;

  double getLength() const { return m_length; }
  double getTfromS(double s) const;

  double getSfromPhi(double phi) const;
  double getPhifromS(double s) const;

  int getPointCount() const { return m_pts.count(); }
  QVector2D getPoint(int index) const { return m_pts[index].pos; }

private:
  void postBuild();
  int getIndexFromPhi(double phi, double *off = 0) const;
  void check1();
};




class AbstractCurve {
  struct Point {
    double s,t;
    QVector2D pos, norm;
    double phi;
  };
  bool m_isOpen;
  QVector<Point> m_pts;

  double m_length;
  double m_phi0;

public:

  AbstractCurve();
  virtual ~AbstractCurve();

  void build(const CurveFunction &f, int n);
  void build(const QVector<QVector2D> &pts);

  virtual QVector2D getPointFromT(double t) const;
  virtual QVector2D getNormalFromT(double t) const;

  double getLength() const { return m_length; }
  double getTfromS(double s) const;

  double getSfromPhi(double phi) const;
  double getPhifromS(double s) const;

  int getPointCount() const { return m_pts.count(); }
  QVector2D getPoint(int index) const { return m_pts[index].pos; }
};


Curve *makeEllipse(double r, double e);
Curve *makeSquare(double radius, double cornerRadius);

Curve *makeConjugate(const Curve *curve, double dist);

#endif

