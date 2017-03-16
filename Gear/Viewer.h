#ifndef VIEWER_H
#define VIEWER_H

#include <gl/glew.h>
#include <QGLWidget>
// #include <QGLFunctions>

#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>

QT_BEGIN_NAMESPACE
class QBasicTimer;
class QGLShaderProgram;
QT_END_NAMESPACE

class Scene;
class Mesh;


class Viewer : public QGLWidget // , protected QGLFunctions
{
  Q_OBJECT
  Scene *m_scene;
  
public:
  explicit Viewer(QWidget *parent = 0);
  virtual ~Viewer();

  void mousePressEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void mouseReleaseEvent(QMouseEvent *e);
  void timerEvent(QTimerEvent *e);

  void initializeGL();
  void resizeGL(int w, int h);
  void paintGL();

  void initShaders();
  void initTextures();
  void initGeometry();
  void draw(Mesh *mesh);

private:
  int m_buttonId;
    QBasicTimer *timer;
    QGLShaderProgram *program;
    // GeometryEngine *geometries;

    GLuint texture;

    QMatrix4x4 projection;

    QVector2D mousePressPosition, lastMousePosition;
    QVector3D rotationAxis;
    qreal angularSpeed;
    QQuaternion rotation;

    // GLuint *vboIds;
    int indicesCount;


};

#endif // VIEWER_H
