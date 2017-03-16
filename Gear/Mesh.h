#ifndef MESH_H
#define MESH_H

#include <gl/glew.h>
#include <QVector>
#include <QVector3D>
#include <QGLWidget>
#include <QMatrix4x4>

class QGLShaderProgram;

struct VertexData
{
    QVector3D position;
    QVector3D normal;    
};


class Mesh 
{
  GLuint m_vboIds[2];
  QVector<VertexData> m_vertices;
  QVector<GLushort> m_indices;  
  QMatrix4x4 m_matrix;
  
public:


  Mesh();
  virtual ~Mesh();
  
  int addVertex(const VertexData &vData) { int index = m_vertices.count(); m_vertices << vData; return index; }
  int addVertex(const QVector3D &position, const QVector3D &normal) { VertexData vd; vd.position = position; vd.normal = normal; return addVertex(vd); }
  
  void addFace(int a, int b, int c) { m_indices << a << b << c; }
  void addFace(int a, int b, int c, int d) { m_indices << a << b << c << a << c << d; }

  void initialize();

  void setMatrix(const QMatrix4x4 &matrix) { m_matrix = matrix; }
  const QMatrix4x4 &getMatrix() const { return m_matrix; }

  void draw(QGLShaderProgram *program);
};


#endif
