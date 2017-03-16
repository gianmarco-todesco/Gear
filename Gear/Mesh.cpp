#include "Mesh.h"
#include <QGLShaderProgram>

Mesh::Mesh()
{
  m_vboIds[0] = 0;
  m_vboIds[1] = 0;
  m_matrix.setToIdentity();
}

Mesh::~Mesh()
{
  glDeleteBuffers(2, m_vboIds); 
  for(int i=0;i<2;i++) m_vboIds[i] = 0;
}


void Mesh::initialize()
{  
  // initializeGLFunctions();
  glGenBuffers(2, m_vboIds);
  
 
  glBindBuffer(GL_ARRAY_BUFFER, m_vboIds[0]);
  glBufferData(GL_ARRAY_BUFFER, m_vertices.count() * sizeof(VertexData), m_vertices.data(), GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_vboIds[1]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.count() * sizeof(GLushort), m_indices.data(), GL_STATIC_DRAW);
}

void Mesh::draw(QGLShaderProgram *program)
{
  glBindBuffer(GL_ARRAY_BUFFER, m_vboIds[0]);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_vboIds[1]);

  int offset = 0;

  int vertexLocation = program->attributeLocation("a_position");
  program->enableAttributeArray(vertexLocation);
  glVertexAttribPointer(vertexLocation, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), (const void *)offset);
    
  offset += sizeof(QVector3D);
    
  vertexLocation = program->attributeLocation("a_normal");
  program->enableAttributeArray(vertexLocation);
  glVertexAttribPointer(vertexLocation, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), (const void *)offset);

  glDrawElements(GL_TRIANGLES, m_indices.count(), GL_UNSIGNED_SHORT, 0);
}
