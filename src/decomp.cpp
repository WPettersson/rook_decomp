#include "stack.h"

#include <algorithm> // std::copy
#include <iostream>
#include <iomanip>  // std::setw
#include <stdlib.h> // atoi(), exit()

#define GRP_SIZE 7
#define GRP_COUNT 5
#define GROUP(i) (((i + (GRP_SIZE-1)) / GRP_SIZE))
#define INDEX(i) (( (i % GRP_SIZE)==0 ? GRP_SIZE : (i%GRP_SIZE)))

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

int orbit_count;

const int WIDTH = 3;

int k; // Rook graph is R_{k,k}

int *adjacencies; // adjacencies[i*k*k + j] gives the orbit of the edge between
                   // vertices i and j


void print_stack(Stack *stack)
{
  for( int i=0; i < orbit_count; i++)
    std::cout << std::setw(WIDTH) << i << ":" << stack->orbits_used[i] << std::endl;
  for( int i=0; i < (k*k); i++)
  {
    int v = stack->vertex_alloc[i];
    if (v == 0)
    {
      std::cout << std::setw(WIDTH) << " ∞ ";
    } else {
      int index = INDEX(v);
      int grp = GROUP(v);
      std::cout << grp << "," << index;
    }
    if (( i > 0 ) && ( i % k ) == k-1 )
      std::cout << std::endl;
    else
      std::cout << " ";
  }
  std::cout << std::endl;
}


void finish(Stack *stack)
{
  print_stack(stack);
  exit(0);
}

// Attempts to place vert_K in spot_R, using Stack *old.
// Returns a new Stack ptr if successful
// Returns NULL if this vertex cannot be placed here.
Stack* use(int spot_R, int vert_K, Stack *old)
{
  // Create new orbits_used
  int *orbits_used = new int[orbit_count];
  std::copy(old->orbits_used, old->orbits_used + orbit_count, orbits_used);
  
  // Use old vertex_alloc for now.
  int *vertex_alloc = old->vertex_alloc;

  int above_R = spot_R % k;
  int left_R = spot_R - above_R;
  for (int j_R = above_R; j_R < spot_R; j_R+=k)
  {
    int nextVert_K = vertex_alloc[j_R];
    int orbit = adjacencies[vert_K*k*k + nextVert_K];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      return NULL;
    }
  }
  for (int j_R = left_R; j_R < spot_R; j_R++)
  {
    int nextVert_K = vertex_alloc[j_R];
    int orbit = adjacencies[vert_K*k*k + nextVert_K];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      return NULL;
    }
  }
  vertex_alloc = new int[k*k];
  std::copy(old->vertex_alloc, old->vertex_alloc + (k*k), vertex_alloc);
  vertex_alloc[spot_R] = vert_K;
  
  bool *vertices_used= new bool[k*k];
  std::copy(old->vertices_used, old->vertices_used + k*k, vertices_used);
  vertices_used[vert_K] = TRUE;

  return new Stack(orbits_used, vertices_used, vertex_alloc);
}


void fill(int spot_R, Stack *stack)
{
  if (spot_R == k*k)
    finish(stack);
  for (int i = 0; i < k*k; i++)
  {
    if (stack->vertices_used[i])
      continue;
    Stack *newStack = use(spot_R, i, stack);
    if (newStack != NULL)
    {
      fill(spot_R+1, newStack);
      delete newStack;
    }
  }
}

Stack* init()
{
  orbit_count = k*k + GRP_COUNT*(GRP_COUNT-1)*GRP_SIZE/2;
  adjacencies = new int[k*k*k*k];
  for (int i_K = 0 ; i_K < k*k; i_K++)
  {
    for (int j_K = i_K+1 ; j_K < k*k; j_K++)
    {
      int grp_i = GROUP(i_K);
      int grp_j = GROUP(j_K);
      int ind_i = INDEX(i_K);
      int ind_j = INDEX(j_K);
      int orbit;
      if (i_K == 0)
      {
        orbit = grp_j * GRP_SIZE;;
      } 
      else if (grp_i == grp_j) 
      {
        orbit = grp_i * GRP_SIZE + (ind_j - ind_i);
      } 
      else
      {
        int index_diff = (ind_j - ind_i)%7;
        if (index_diff < 0) index_diff += 7;
        orbit = k*k-1   // Skip orbits including infinity or inside a group
                + GRP_SIZE*(GRP_COUNT*(grp_i-1) - ((grp_i-1)*(grp_i-1) + (grp_i-1))/2)
                      // Skip orbits from groups before grp_i
                + (grp_j-grp_i-1)*GRP_SIZE  // Skip orbits from grp_i to groups before grp_j
                + index_diff; // Find right orbit
      }
      //std::cout << "Orbit between " << std::setw(3) << i_K << " and " <<
      //  std::setw(3) << j_K << " : " << orbit << std::endl;
      //std::cout << "Orbit between " << grp_i <<","<<ind_i << " and " <<
      //  grp_j<<","<<ind_j << " : " << orbit << std::endl;
      if (orbit >= orbit_count)
        std::cout << "Bad orbit between " << i_K << " and " << j_K << " : " <<
          orbit << std::endl;
      adjacencies[i_K*k*k + j_K] = orbit;
      adjacencies[j_K*k*k + i_K] = orbit;
    }
  }
//  for (int i = 0 ; i < k*k; i++)
//  {
//    for (int j = 0 ; j < k*k; j++)
//    {
//      std::cout << std::setw(4) << adjacencies[i*k*k + j] << " ";
//    }
//    std::cout << std::endl;
//  }





  int *orbits_used = new int[orbit_count];
  bool *vertices_used = new bool[k*k];
  int *vertex_alloc = new int[k*k];

  for( int i=0; i< k*k; i++)
  {
    vertices_used[i] = FALSE;
  }
  for (int i=0; i < orbit_count; i++)
    orbits_used[i] = 0;

  Stack *s = new Stack(orbits_used, vertices_used, vertex_alloc);
  return s;
}

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    std::cout << "Usage: decomp k [alloc1 [alloc2 ... ]]" << std::endl;
    std::cout << " where vertex alloc1 will be put in spot 0 and so-on" << std::endl;
    return 0;
  }
  k = atoi(argv[1]);

  Stack *stack = init();
  int spot_R = 0;
  for (int i = 2; i < argc; i++) 
  {
    Stack *newStack = use(spot_R, atoi(argv[i]), stack);
    if ( newStack == NULL)
    {
      std::cout << "Your allocation failed when putting vert_K " <<
        atoi(argv[i]) <<
        " into position_R " << spot_R << std::endl;
      return -1;
    }
    spot_R++;
    stack = newStack;
  }

  fill(spot_R, stack);

  return -1; // Did not find decomp.
}
