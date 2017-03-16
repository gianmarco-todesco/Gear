#ifndef SCENE_H
#define SCENE_H


class AbstractGear;
class GearLink;
class Viewer;
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector>


class Scene
{
  QVector<AbstractGear*> m_gears;
  QVector<GearLink*> m_links;

public:
  QQuaternion m_rotation;

  Scene();
  ~Scene();

  void drag(int dx, int dy);
  void initialize();
  void draw(Viewer *viewer);
};

#endif

