#ifndef STACK_H
#define STACK_H


class Stack
{
  public:
    Stack(int *orbits_used, bool *vertices_used, int *vertex_alloc);
    //Stack(Stack copy);
    ~Stack();
    int *orbits_used;
    bool *vertices_used;
    int *vertex_alloc;
};

#endif /* STACK_H */
