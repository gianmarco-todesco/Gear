#include "Viewer.h"

#include <QBasicTimer>
#include <QMouseEvent>
#include <QDebug>

#include <QGLShaderProgram>


#include <qmath.h>
#include <locale.h>

#include "Scene.h"
#include "Mesh.h"

// see http://stackoverflow.com/questions/11845230/glgenbuffers-crashes-in-release-build

//------------------------------------------------



//------------------------------------------------
/*
Gear1(int n)  { 
  tooth.count = n; 
  tooth.d = 0.2; 
  tooth.h = 0.1; 
  radius = tooth.d * tooth.count / (2* M_PI); 
  build();
};

  float getRadius() const { return radius; }
};

void Gear1::buildGeometry(QVector<VertexData> &vertices, QVector<GLushort> &indices)
{
  QVector<QPair<float, float> > ppts;
  for(int i=0; i<tooth.count; i++)
  {
    float dphi = 2 * M_PI / tooth.count; 
    float phi = i * dphi;
    ppts.append(qMakePair(phi, radius));
    ppts.append(qMakePair(phi+dphi*0.5f*0.05f, radius));
    ppts.append(qMakePair(phi+dphi*0.5f*0.3f, radius+tooth.h));
    ppts.append(qMakePair(phi+dphi*0.5f*0.7f, radius+tooth.h));
    ppts.append(qMakePair(phi+dphi*0.5f*0.95f, radius));
    ppts.append(qMakePair(phi+dphi*0.5f, radius));
  }

  int m = ppts.size();
  QVector<QVector3D> pts;
  for(int i=0; i<m; i++)
  {
    double phi = ppts[i].first;
    double csPhi = cos(phi), snPhi = sin(phi);
    pts.append(QVector3D(cos(phi), sin(phi),0));
  }

  QVector3D dz(0,0,0.05);
  for(int i=0; i<m; i++)
  {
    VertexData v;
    v.normal = QVector3D(0,0,1);
    v.position = 0.2*pts[i] + dz;
    vertices.push_back(v);

    v.position = ppts[i].second*pts[i] + dz;
    vertices.push_back(v);

    QVector3D w = (pts[(i+1)%m] - pts[(i+m-1)%m]);
    v.normal = -QVector3D::crossProduct(QVector3D(0,0,1),w).normalized();
    vertices.push_back(v);
    
    v.position -= 2*dz; 
    vertices.push_back(v);
  }

  
  for(int i=0;i<=m;i++)
  {
    int j = i%m;
    indices.push_back(j*4);
    indices.push_back(j*4+1);
  }
  indices.push_back(0);

  indices.push_back(3);
  for(int i=0;i<=m;i++)
  {
    int j = i%m;
    indices.push_back(j*4+2);
    indices.push_back(j*4+3);
  }
}
*/

//------------------------------------------------


Viewer::Viewer(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
    , timer(new QBasicTimer)
    , program(new QGLShaderProgram)
    , indicesCount(0)
    , m_scene(new Scene())
    , m_buttonId(Qt::NoButton)
{
  setAutoFillBackground(false);
}

Viewer::~Viewer()
{
  delete timer; timer = 0;
  delete program; program = 0;
  //delete geometries; geometries = 0;

  //deleteTexture(texture);
  delete m_scene;

  
}

void Viewer::mousePressEvent(QMouseEvent *e)
{
    // Saving mouse press position
    lastMousePosition = mousePressPosition = QVector2D(e->posF());
    m_buttonId = e->button();
}

void Viewer::mouseMoveEvent(QMouseEvent *e)
{
    // Mouse release position - mouse press position
    QVector2D diff = QVector2D(e->posF()) - lastMousePosition;
    lastMousePosition = QVector2D(e->posF());

    if(m_buttonId == Qt::LeftButton)
    {
      m_scene->drag(diff.x(), diff.y());
    }
    else
    {

      // Rotation axis is perpendicular to the mouse position difference 
      // vector
      QVector3D n = QVector3D(diff.y(), diff.x(), 0.0).normalized();
      m_scene->m_rotation = (QQuaternion::fromAxisAndAngle(n, diff.length()*1.0) * m_scene->m_rotation).normalized(); 
    }

    
}

void Viewer::mouseReleaseEvent(QMouseEvent *e)
{
  return;
    // Mouse release position - mouse press position
    QVector2D diff = QVector2D(e->posF()) - mousePressPosition;

    // Rotation axis is perpendicular to the mouse position difference 
    // vector
    QVector3D n = QVector3D(diff.y(), diff.x(), 0.0).normalized();

    // Accelerate angular speed relative to the length of the mouse sweep
    qreal acc = diff.length() / 100.0;

    // Calculate new rotation axis as weighted sum
    rotationAxis = (rotationAxis * angularSpeed + n * acc).normalized();

    // Increase angular speed
    angularSpeed += acc;
}



void Viewer::timerEvent(QTimerEvent *e)
{
    Q_UNUSED(e);
    //gearRotation+= 0.2;
    updateGL();
    return;

    // Decrease angular speed (friction)
    angularSpeed *= 0.99;

    // Stop rotation when speed goes below threshold
    if (angularSpeed < 0.01)
        angularSpeed = 0.0;
    else {
        // Update rotation
        rotation = QQuaternion::fromAxisAndAngle(rotationAxis, angularSpeed) * rotation;

        // Update scene
        updateGL();
    }
}

void Viewer::initializeGL()
{
  glewInit();
    // initializeGLFunctions();
    // setFormat(QGLFormat(QGL::SampleBuffers));

    glClearColor(0.8,0.8,0.9,1.0);

    qDebug() << "Initializing shaders...";
    initShaders();

    // qDebug() << "Initializing textures...";
    // initTextures();

//! [2]
    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);

    // Enable back face culling
    glEnable(GL_CULL_FACE);
//! [2]

    qDebug() << "Initializing geometries...";
    // geometries->init();

    // glGenBuffers(2, vboIds);
    
    m_scene->initialize();

    // using QBasicTimer because its faster that QTimer
    timer->start(12, this);
    qDebug() << "started";
}

//! [3]
void Viewer::initShaders()
{
      qDebug() << "i0";

    // Overriding system locale until shaders are compiled
    setlocale(LC_NUMERIC, "C");

    // Compiling vertex shader
    if (!program->addShaderFromSourceFile(QGLShader::Vertex, "./vshader.glsl")) // :/
    {
      
      qDebug() << program->log();
      abort();
    }
    qDebug() << "vshader ok";      

    // Compiling fragment shader
    if (!program->addShaderFromSourceFile(QGLShader::Fragment, "./fshader.glsl"))
    {
      qDebug() << program->log();
      abort();
    }
    qDebug() << "fshader ok";      
    // Linking shader pipeline
    if (!program->link())
    {
      qDebug() << program->log();
      abort();
    }
    qDebug() << "program ok";      

    // Binding shader pipeline for use
    if (!program->bind())
        close();
    int offset = 0;

    
    // Restore system locale
    setlocale(LC_ALL, "");
}


void Viewer::resizeGL(int w, int h)
{
    // Set OpenGL viewport to cover whole widget
    glViewport(0, 0, w, h);

    // Calculate aspect ratio
    qreal aspect = (qreal)w / ((qreal)h?h:1);

    // Set near plane to 3.0, far plane to 7.0, field of view 45 degrees
    const qreal zNear = 0.5, zFar = 30.0, fov = 15.0;

    // Reset projection
    projection.setToIdentity();

    // Set perspective projection
    projection.perspective(fov, aspect, zNear, zFar);
}
//! [5]


void Viewer::draw(Mesh *mesh)
{

  const QMatrix4x4 &matrix = mesh->getMatrix();

  program->setUniformValue("mv_matrix", matrix);
  program->setUniformValue("normal_matrix", matrix.normalMatrix());
  program->setUniformValue("mvp_matrix", projection * matrix);
  mesh->draw(program);
}

void Viewer::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    m_scene->draw(this);

    //glColor3d(1,1,1);
    //renderText(50,50,"Ciao");
}
