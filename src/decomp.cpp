#include "stack.h"

#include <algorithm> // std::copy
#include <iostream>
#include <iomanip>  // std::setw
#include <cstdlib> // atoi(), exit()
#include <ctime> // clock()
#include <cstring> // strcmp()

int orbit_count;

const int WIDTH = 3;

int k; // Rook graph is R_{k,k}
int grp_size;
int grp_count;

int *adjacencies; // adjacencies[i*k*k + j] gives the orbit of the edge between
                   // vertices i and j

int count; // Count how many such decompositions we find

bool *prealloc; // Which vertices were assigned from the command line.
int *progress; // Track how we are faring.
bool first; // Is this the first time we're printing the progress (in which
            // case, don't scroll up

time_t last_progress;
const time_t progress_interval = 30*60*CLOCKS_PER_SEC; // Only print progress every "progress_interval"

inline int group(int i)
{
  return (i + grp_size-1)/grp_size;
}

inline int index(int i)
{
  return ( (i % grp_size == 0) ? grp_size: i % grp_size);
}

void print_progress()
{
  if (!first)
    std::cout << "\033[" << (k+1) << "F" << std::endl;
  first = false;
  for (int i = 0; i < k*k; i++)
  {
    if (prealloc[i])
      std::cout << std::setw(WIDTH) << "  X";
    else
      std::cout << std::setw(WIDTH) << progress[i];
    if (( i > 0 ) && ( i % k ) == k-1 )
      std::cout << std::endl;
    else
      std::cout << " ";
  }
}

void print_orbits(Stack *stack)
{
  for( int i=0; i < orbit_count; i++)
    std::cout << std::setw(WIDTH) << i << ":" << stack->orbits_used[i] << std::endl;
}

void print_forced(Stack *stack)
{
  for( int i=0; i < (k*k); i++)
  {
    int v = stack->vertex_alloc[i];
    if (! prealloc[i])
    {
      std::cout << std::setw(WIDTH) << " _ ";
    }
    else if (v == 0)
    {
      std::cout << std::setw(WIDTH) << " ∞ ";
    }
    else
    {
      int ind = index(v);
      int grp = group(v);
      std::cout << grp << "," << ind;
    }
    if (( i > 0 ) && ( i % k ) == k-1 )
      std::cout << std::endl;
    else
      std::cout << " ";
  }
  std::cout << std::endl;
}


void print_stack(Stack *stack)
{
  for( int i=0; i < (k*k); i++)
  {
    int v = stack->vertex_alloc[i];
    if (v == 0)
    {
      std::cout << std::setw(WIDTH) << " ∞ ";
    } else {
      int ind = index(v);
      int grp = group(v);
      std::cout << grp << "," << ind;
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
  count+=1;
  //print_orbits(stack);
  print_stack(stack);

  // Only exit if count = 0, which means we are not looking for all
  if (count == 0)
    exit(0);

  // If we aren't exiting, set first to true so that progress prints
  // don't overwrite the outputs
  first = true;
  return;
}

// Attempts to place vert_K in spot_R, using Stack *old.
// Returns a new Stack ptr if successful
// Returns NULL if this vertex cannot be placed here.
Stack* use(int spot_R, int vert_K, Stack *old)
{
  if (old->vertices_used[vert_K])
    return NULL;
  // Create new orbits_used
  int *orbits_used = new int[orbit_count];
  std::copy(old->orbits_used, old->orbits_used + orbit_count, orbits_used);
  
  // Use old vertex_alloc for now.
  int *vertex_alloc = old->vertex_alloc;

  int col = spot_R % k;
  int row = (spot_R - col) / k;
  int vert_to_check;
  int orbit;
  if (col > 0) // Check to left
  {
    vert_to_check = vertex_alloc[spot_R - 1];
    orbit = adjacencies[vert_K*k*k + vert_to_check];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
      return NULL;
    }
  }
  if (col == k-1) // Check wrap around if at last column
  {
    vert_to_check = vertex_alloc[spot_R - k + 1];
    orbit = adjacencies[vert_K*k*k + vert_to_check];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
      return NULL;
    }
  }

  if (row > 0) // Check above, and above-left.
  {
    // Above
    vert_to_check = vertex_alloc[spot_R - k];
    orbit = adjacencies[vert_K*k*k + vert_to_check];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
      return NULL;
    }
    if (col == 0) // Above-left is above and far-right
      vert_to_check = vertex_alloc[spot_R - 1];
    else // Above-left is just above-left.
      vert_to_check = vertex_alloc[spot_R - k - 1];
    orbit = adjacencies[vert_K*k*k + vert_to_check];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
      return NULL;
    }
  }
  if (row == k-1) // Check bottom-to-top, and bottom-to-top-and-left
  {
    // Bottom-to-top
    vert_to_check = vertex_alloc[col];
    orbit = adjacencies[vert_K*k*k + vert_to_check];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
      return NULL;
    }
    // Bottom-to-top-and-left
    if (col == 0) // Top-and-above-left is top and far-right
      vert_to_check = vertex_alloc[k - 1];
    else // Top-and-above-left is just top-and-above-left.
      vert_to_check = vertex_alloc[col - 1];
    orbit = adjacencies[vert_K*k*k + vert_to_check];
    if (++orbits_used[orbit] > 2)
    {
      delete[] orbits_used;
      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
      return NULL;
    }
  }

//  int above_R = spot_R % k;
//  int left_R = spot_R - above_R;
//  for (int j_R = above_R; j_R < spot_R; j_R+=k)
//  {
//    int nextVert_K = vertex_alloc[j_R];
//    int orbit = adjacencies[vert_K*k*k + nextVert_K];
//    if (++orbits_used[orbit] > 2)
//    {
//      delete[] orbits_used;
//      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
//      return NULL;
//    }
//  }
//  for (int j_R = left_R; j_R < spot_R; j_R++)
//  {
//    int nextVert_K = vertex_alloc[j_R];
//    int orbit = adjacencies[vert_K*k*k + nextVert_K];
//    if (++orbits_used[orbit] > 2)
//    {
//      delete[] orbits_used;
//      //std::cout << "Orbit " << orbit << " used too many times" << std::endl;
//      return NULL;
//    }
//  }
  vertex_alloc = new int[k*k];
  std::copy(old->vertex_alloc, old->vertex_alloc + (k*k), vertex_alloc);
  vertex_alloc[spot_R] = vert_K;
  
  bool *vertices_used= new bool[k*k];
  std::copy(old->vertices_used, old->vertices_used + k*k, vertices_used);
  vertices_used[vert_K] = true;

  return new Stack(orbits_used, vertices_used, vertex_alloc);
}


void fill(int spot_R, Stack *stack)
{
#ifndef QUIET
  if (clock() - last_progress > progress_interval)
  {
    last_progress = clock();
    print_progress();
  }
  progress[spot_R] = 0;
#endif
  if (spot_R == k*k)
    finish(stack);
  for (int i = 0; i < k*k; i++)
  {
#ifndef QUIET
    progress[spot_R]=i;
#endif
    if (stack->vertices_used[i])
      continue;
    Stack *newStack = use(spot_R, i, stack);
    if (newStack != NULL)
    {
      //std::cout << "Put " << group(i) << "," << index(i) << " into " << spot_R << std::endl;
      int newSpot = spot_R+1;
      while (prealloc[newSpot])
        newSpot++; // While loop to find next not-preallocated vertex
      fill(newSpot, newStack);
      delete newStack;
    }
  }
}

Stack* init()
{
  last_progress = clock();
  grp_size = k + 1;
  grp_count = k - 1;
  orbit_count = k*k +1  + grp_count*(grp_count-1)*grp_size/2;
  adjacencies = new int[k*k*k*k];
  for (int i_K = 0 ; i_K < k*k; i_K++)
  {
    for (int j_K = i_K+1 ; j_K < k*k; j_K++)
    {
      int grp_i = group(i_K);
      int grp_j = group(j_K);
      int ind_i = index(i_K);
      int ind_j = index(j_K);
      int orbit;
      if (i_K == 0)
      {
        orbit = grp_j * grp_size;;
      } 
      else if (grp_i == grp_j) 
      {
        int diff = ind_j - ind_i;
        if ((diff > grp_size/2))
        {
          diff = grp_size - diff;
        }
        orbit = (grp_i-1) * grp_size + diff;
      } 
      else
      {
        int index_diff = (ind_j - ind_i)%grp_size;
        if (index_diff < 0) index_diff += grp_size;
        orbit = k*k+1   // Skip orbits including infinity or inside a group
                + grp_size*(grp_count*(grp_i-1) - ((grp_i-1)*(grp_i-1) + (grp_i-1))/2)
                      // Skip orbits from groups before grp_i
                + (grp_j-grp_i-1)*grp_size  // Skip orbits from grp_i to groups before grp_j
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
    vertices_used[i] = false;
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
    std::cout << "Usage: decomp k [-a] [spot 1 alloc1 [spot 2 alloc2 [ ...] ]]" << std::endl;
    std::cout << " where vertex alloc1 will be put in spot1 and so-on" << std::endl;
    std::cout << " [-a] indicates that the program should attempt to find " << std::endl;
    std::cout << " all such decompositions, and count them." << std::endl;
    return 0;
  }
  k = atoi(argv[1]);
  if ((argc >= 3) && ( strcmp("-a", argv[2])==0))
  {
    count = 0;
  }
  else
    count = -1; // -1 indicates that we should exit on finding one.

#ifndef QUIET
  first = true;
  progress = new int[k*k];
#endif
  prealloc = new bool[k*k];
  for (int i = 0; i < k*k; i++)
    prealloc[i] = false;

  Stack *stack = init();
  int spot_R = 0;
  for (int i = 3 + count; i < argc; i+=2) // Note that count is 0 if "-a" is
    // passed in, and -1 otherwise, so this starts at the right spot
  {
    int spot = atoi(argv[i]);
    int vert = atoi(argv[i+1]);
    Stack *newStack = use(spot, vert, stack);
    if ( newStack == NULL)
    {
      int grp = group(vert);
      int ind = index(vert);
      std::cout << "Your allocation failed when putting vert_K " <<
        grp << "," << ind << " into position_R " << spot << std::endl;
      //print_orbits(stack);
      print_forced(stack);
      return -1;
    }
    prealloc[spot] = true;
    stack = newStack;
  }
  std::cout << "Forced allocations" << std::endl;
  print_forced(stack);

  while (prealloc[++spot_R]) ; // Empty while loop to find next
                                // not-preallocated vertex
  //std::cout << "Starting at " << spot_R << std::endl;
  fill(spot_R, stack);

#ifndef QUIET
  print_progress();
#endif
  if (count > 1)
    std::cout << "Found " << count << " decompositions." << std::endl;
  return -1; // Did not find decomp.
}
