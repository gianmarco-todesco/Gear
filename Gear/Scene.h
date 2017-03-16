#ifndef SCENE_H
#define SCENE_H


class Gear;
class GearLink;
class Viewer;
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector>


class Scene
{
  QVector<Gear*> m_gears;
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

