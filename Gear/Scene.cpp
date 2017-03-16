#include "Scene.h"
#include "Gear.h"
#include "Viewer.h"
#include <qmath.h>

Scene::Scene()
{
  
  if(false)
  {
    double r = 0.5, e = 0.7;
    int m = 4;
    double dist = r/(1.0+e)+r/(1.0-e);

    for(int i=0;i<m;i++)
    {
      Gear *gear = new Gear(makeEllipse(r,e));
      gear->setPosition( QVector3D(-dist+dist*i,0,0));
      m_gears.append(gear);
    }
  
  

    for(int i=1;i<4;i++)
      m_links.append(new GearLink(m_gears[i-1], m_gears[i]));
  }
  else
  {
    m_gears.append(new Gear(makeEllipse(0.5,0.6)));   
    m_gears.append(new Gear(makeSquare(0.5,0.3)));
    m_gears[1]->setPosition(QVector3D(1,0,0));

  }
  // m_gears.append(new Gear(makeSquare(256,1.0,0.1)));
  // m_gears.append(new Gear(makeEllipse(0.5,0.6)));

  // double dist = 1.2;

 // m_gears.append(makeConjugate(m_gears[0], dist));

  //m_gears[0]->setPosition(QVector3D(-dist/2,0,0));
  //m_gears[1]->m_position = QVector3D( dist/2,0,0);

  //m_links.append(new GearLink(m_gears[0], m_gears[1]));

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
    Gear *gear = m_gears[i];
    QMatrix4x4 matrix;
    matrix.translate(0.0, 0.0, -cameraDist);
    matrix.rotate(m_rotation);
    matrix.translate(gear->getPosition());
    matrix.rotate(gear->getRotation(), 0,0,1);
    gear->setMatrix(matrix);
    viewer->draw(gear);
  }
}

void Scene::drag(int dx, int dy)
{
  m_gears[0]->setRotation(m_gears[0]->getRotation() + dx*0.5);
  for(int i=0;i<m_links.count();i++)
  {
    m_links[i]->update();
  }

}

