#include "Scene.h"
#include "Gear.h"
#include "Viewer.h"
#include <qmath.h>

Scene::Scene()
{
  /*
  double r = 0.5, e = 0.7;
  m_gears.append(new AbstractGear(makeEllipse(256,r,e)));
  m_gears.append(new AbstractGear(makeEllipse(256,r,e)));
  m_gears.append(new AbstractGear(makeEllipse(256,r,e)));
  m_gears.append(new AbstractGear(makeEllipse(256,r,e)));

  double dist = 0.5/(1.0+e)+0.5/(1.0-e);

  m_gears[0]->m_position = QVector3D(-dist,0,0);
  m_gears[1]->m_position = QVector3D(    0,0,0);
  m_gears[2]->m_position = QVector3D( dist,0,0);
  m_gears[3]->m_position = QVector3D( dist*2,0,0);

  m_links.append(new GearLink(m_gears[0], m_gears[1]));
  m_links.append(new GearLink(m_gears[1], m_gears[2]));
  m_links.append(new GearLink(m_gears[2], m_gears[3]));
  */
  m_gears.append(new AbstractGear(makeSquare(256,1.0,0.1)));

  double dist = 1.2;

  m_gears.append(makeConjugate(m_gears[0], dist));

  m_gears[0]->m_position = QVector3D(-dist/2,0,0);
  m_gears[1]->m_position = QVector3D( dist/2,0,0);

  m_links.append(new GearLink(m_gears[0], m_gears[1]));

  for(int i=0;i<m_gears.count();i++)
  {
    m_gears[i]->build();
  }
  for(int i=0;i<m_links.count();i++)
  {
    m_links[i]->update();
  }


  // m_rotation;
}


Scene::~Scene()
{
  for(int i=0;i<m_gears.count();i++) delete m_gears[i];
}

void Scene::initialize()
{
  for(int i=0;i<m_gears.count();i++) m_gears[i]->initialize();
}

void Scene::draw(Viewer *viewer)
{
  double cameraDist = 10;
  for(int i=0;i<m_gears.count();i++)
  {
    AbstractGear *gear = m_gears[i];
    QMatrix4x4 matrix;
    matrix.translate(0.0, 0.0, -cameraDist);
    matrix.rotate(m_rotation);
    matrix.translate(gear->m_position);
    matrix.rotate(gear->m_rotation, 0,0,1);
    gear->setMatrix(matrix);
    viewer->draw(gear);
  }
}

void Scene::drag(int dx, int dy)
{
  m_gears[0]->m_rotation += dx*0.5;
  for(int i=0;i<m_links.count();i++)
  {
    m_links[i]->update();
  }

}

