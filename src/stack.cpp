#include "stack.h"

Stack::Stack(int *_orbits_used, bool *vertices_used, int *vertex_alloc) :
  orbits_used(_orbits_used), vertices_used(vertices_used),
  vertex_alloc(vertex_alloc)
{
}


Stack::~Stack()
{
  delete[] orbits_used;
  delete[] vertices_used;
  delete[] vertex_alloc;
}
