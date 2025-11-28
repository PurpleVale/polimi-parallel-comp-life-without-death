%%writefile life_wout_death.cpp

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <omp.h>

#define MAX_STEPS 1000000 //* simulation runs for that generations a t max
#define NUM_GROUPS 7
#define COLOR_VARIATION 0.08f
#define CLUSTER_FILL_DENSITY 0.55f
#define CLUSTER_RADIUS_FACTOR 0.08f
#define CLUSTER_MIN_RADIUS 3

// === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
// =================================================

// if you need any additional standard library up to C11 or C++17, add it here

// if you need any additional define, add it here

// toggle this to save time when profiling
// be wary: this also disables the results check
#define DISABLE_SEQUENTIAL false

// if you happen to need extra information in the grid, you can modify this object, however,
// only add, do not remove anything, as to preserve compatibility with the sequential version!
//
// notes about the use of C++ here:
// - Grid is a C++ object owning two heap-allocated arrays
// - copying or assigning a Grid (Grid B = A; or B = A;) performs a deep copy of both vectors
// - passing Grid by value also triggers this deep copy, thus always pass by reference
// => this is different from C: the struct looks cheap, but copies are expensive
struct Grid {
  int W, H;
  //* value 1 for liveness and 0 for dead
  //! disposed in row major order
  std::vector<unsigned char> alive; // size: W*H
  std::vector<float> hue; // size: W*H
  // constructor: allocates and owns two dynamic arrays via std::vector
  Grid(int w, int h) : W(w), H(h), alive(w*h,0), hue(w*h,0.0f) {}
};

// =================================================
// === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

// game rules B368 / S012345678:
// - a dead cell is born if alive neighbor count is 3, 6, or 8
// - a live cell survives with ANY count 0..8
inline bool birth_rule(int n) {
  return (n == 3 || n == 6 || n == 8); //* n represents the number of alive neighbors
}
//! For us it is not needed
inline bool survive_rule(int n) {
  return (n >= 0 && n <= 8); //* n represents the number of dead neighbors
}

// HELPER FUNCTION [YOU CAN USE IT AS-IS]
// hue averaging using circular mean around the unit circle
//* When a cell is born its color is the average of its parents
float hue_average(const float *h, int n) {
  float sumx = 0.0f, sumy = 0.0f;
  for (int i = 0; i < n; i++) { //* iterate for all the neighboring points
    float angle = h[i] * 2.0f * M_PI; //* convert the hue into radiants
    sumx += cosf(angle);
    sumy += sinf(angle);
  }
  float angle = atan2f(sumy, sumx); //* transform the point obtained into an angle
  if (angle < 0) angle += 2.0f * M_PI;
  return angle / (2.0f * M_PI);
}

// HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
// random hue selection
float pick_group_hue(int group_id) {
  float base = float(group_id) / NUM_GROUPS;
  float h = base + ((float)rand() / RAND_MAX)*2*COLOR_VARIATION - COLOR_VARIATION;
  if (h < 0) h += 1.0f;
  if (h >= 1) h -= 1.0f;
  return h;
}

// HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
// grid initialization with a few clusters
void initialize_grid(Grid &g) {
  int W = g.W, H = g.H;

  //* calculate the initial cluster size (same for each one)
  int min_dim = (W < H ? W : H);
  int r = int(min_dim * CLUSTER_RADIUS_FACTOR);
  if (r < CLUSTER_MIN_RADIUS) r = CLUSTER_MIN_RADIUS;

  //* create an array of hues
  std::vector<float> group_hue(NUM_GROUPS);
  for (int i = 0; i < NUM_GROUPS; i++)
    group_hue[i] = pick_group_hue(i);

  //* Define the center point of each cluster
  std::vector<int> cx(NUM_GROUPS), cy(NUM_GROUPS);
  for (int i = 0; i < NUM_GROUPS; i++) {
    cx[i] = rand() % W;
    cy[i] = rand() % H;
  }


  for (int g_id = 0; g_id < NUM_GROUPS; g_id++) {
    float base_h = group_hue[g_id]; //* take the color
    int gx = cx[g_id], gy = cy[g_id]; //* take the central point

    //* populate the cluster with the points
    for (int dx = -r; dx <= r; dx++) {
      for (int dy = -r; dy <= r; dy++) {
        if (dy*dy + dx*dx <= r*r) { //* So we obtain a cirle shape
          int x = (gx + dx + 2*W) % W, y = (gy + dy + 2*H) % H; //* handle the edge points (add the *2W  because enforces a positive number without changing the relative position)
          if (((float)rand()/RAND_MAX) < CLUSTER_FILL_DENSITY) { //* 55% of chances of coloring a point inside the cluster perimeter
            int idx = y*W + x;
            g.alive[idx] = 1; //* Add an active point on the grid
            g.hue[idx] = base_h; //* Add the color
          }
        }
      }
    }
  }
}

// HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
// compares two grids
bool compare_grids(const Grid &a, const Grid &b) {
  int N = a.W * a.H;
  for (int i = 0; i < N; i++) {
    if (a.alive[i] != b.alive[i]) return false;
    if (a.alive[i]) {
      float ha = a.hue[i];
      float hb = b.hue[i];
      // compare hue with a tolerance on floats (due to associativity)
      if (fabs(ha - hb) > 1e-4f) return false;
    }
  }
  return true;
}

// HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
// write grid to file for later visualization :)
void write_grid_to_file(const Grid &g, const char *filename) {
  FILE *f = fopen(filename, "wb");
  if (!f) {
    fprintf(stderr, "ERROR: cannot open %s for writing\n", filename);
    return;
  }
  uint32_t W = g.W;
  uint32_t H = g.H;
  fwrite(&W, sizeof(uint32_t), 1, f);
  fwrite(&H, sizeof(uint32_t), 1, f);
  for (size_t i = 0; i < W * H; i++)
    fwrite(&g.alive[i], sizeof(uint8_t), 1, f);
  for (size_t i = 0; i < W * H; i++)
    fwrite(&g.hue[i], sizeof(float), 1, f);
  fclose(f);
}

// advance the simulation by one step; return the count of changed cells
int evolve_step(const Grid &cur, Grid &next) { //* need to pass the current and the next grid
  int W = cur.W, H = cur.H;
  int changes = 0;

  // iterate over all cells
  for (int x = 0; x < W; x++) {
    for (int y = 0; y < H; y++) {
    int idx = y*W + x;
    unsigned char alive = cur.alive[idx];

    int alive_neighbors = 0;
    float parent_hues[8]; //* neighboring colors

    //! Possiamo mettere questo ciclo dopo e farlo andare solo se la cella per il momento non è alive?
    //* count alive moore neighbors and collect their hues
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        if (dx == 0 && dy == 0) continue; //* skip the current cell because it is itself
        // wrap around the grid (torus)
        int xx = (x + dx + W) % W;
        int yy = (y + dy + H) % H;
        int nidx = yy*W + xx; //* row major global id
        if (cur.alive[nidx]) {
          parent_hues[alive_neighbors] = cur.hue[nidx];
          alive_neighbors++;
        }
      }
    }

    unsigned char new_alive = alive;
    if (!alive) { //* Look if the current cell is alive or not
      if (birth_rule(alive_neighbors)) {
        new_alive = 1;
        next.hue[idx] = hue_average(parent_hues, alive_neighbors); //* pass the hue to the new grid
      }
    } else {
      // you are free to skip this check, here it was kept for sake of completeness
      if (survive_rule(alive_neighbors)) {
        new_alive = 1;
        next.hue[idx] = cur.hue[idx];
      } else {
        new_alive = 0;
      }
    }

    next.alive[idx] = new_alive;
    if (new_alive != alive)
      changes++;
    }
  }
  return changes;
}

// sequential simulation entry point
void simulate_sequential(Grid &g) {
  // it's hard to perform updates in place, we ping-pong between grid copies
  Grid tmp(g.W, g.H);

  for (long step = 0; step < MAX_STEPS; step++) {
    // one step at a time
    // note: fully overwrites tmp with the next state
    int changes = evolve_step(g, tmp);

    // swap g <-> tmp
    g.alive.swap(tmp.alive);
    g.hue.swap(tmp.hue);

    if (changes == 0) break;
  }
}

// === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
// =================================================

// ===> IMPLEMENT YOUR VERSION OF THIS FUNCTION <===

int evolve_step_parallel(const Grid &cur, Grid &next) { //* need to pass the current and the next grid
  int W = cur.W, H = cur.H;
  int changes = 0;

  // iterate over all cells
  # pragma omp parallel for collapse(2) reduction(+:changes)
  for (int x = 0; x < W; x++) {
    for (int y = 0; y < H; y++) {
    int idx = y*W + x;
    unsigned char alive = cur.alive[idx];

    if (alive) {
      next.alive[idx] = 1;
      next.hue[idx] = cur.hue[idx];
    }
    else {
      int alive_neighbors = 0;
      float parent_hues[8];

      for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        if (dx == 0 && dy == 0) continue;

        int xx = (x + dx + W) % W;
        int yy = (y + dy + H) % H;
        int nidx = yy*W + xx;

        if (cur.alive[nidx]) {
          parent_hues[alive_neighbors] = cur.hue[nidx];
          alive_neighbors++;
        }
      }
    }

    if (birth_rule(alive_neighbors)){
        next.alive[idx] = 1;
        next.hue[idx] = hue_average(parent_hues, alive_neighbors);
        changes++;
      } else {
        next.alive[idx] = 0;
      }
    }
  }
}
  return changes;
}

void simulate_parallel(Grid &g) {
  Grid tmp(g.W, g.H);

  for (long step = 0; step < MAX_STEPS; step++) {
    // one step at a time
    // note: fully overwrites tmp with the next state
    int changes = evolve_step_parallel(g, tmp);

    // swap g <-> tmp
    g.alive.swap(tmp.alive);
    g.hue.swap(tmp.hue);

    if (changes == 0) break;
  }
}

// =================================================
// === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

int main(int argc, char **argv) {
    omp_set_nested(true); // just in case, enable nested parallelism

    if (argc < 4 || argc > 5) {
        fprintf(stderr, "Usage: %s <grid-width> <grid-height> <seed> <[opt]-output-filename>\n", argv[0]);
        return 1;
    }

    int W = atoi(argv[1]);
    int H = atoi(argv[2]);
    int seed = atoi(argv[3]);
    srand(seed);

    Grid gs(W,H), gp(W,H);
    initialize_grid(gs);
    // copy the initial state
    gp = gs;

    // === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
    // =================================================

    // if you need to initialize additional things, do that here!

    // =================================================
    // === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

    // sequential
    #if !DISABLE_SEQUENTIAL
    double t1 = omp_get_wtime();
    simulate_sequential(gs);
    double t2 = omp_get_wtime();
    #endif

    // parallel
    double t3 = omp_get_wtime();
    simulate_parallel(gp);
    double t4 = omp_get_wtime();

    #if !DISABLE_SEQUENTIAL
    printf("Sequent. time: %.6f s\n", t2 - t1);
    #endif
    printf("Parallel time: %.6f s\n", t4 - t3);

    #if !DISABLE_SEQUENTIAL
    bool equal = compare_grids(gs, gp);
    printf("Results check: %s\n", equal ? "PASS" : "FAIL");
    #endif

    // === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
    // =================================================

    // if you need more logging, put it here!

    // =================================================
    // === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

    if (argc == 5) {
      write_grid_to_file(gp, (char *)argv[4]);
      printf("Results written to %s\n", (char *)argv[4]);
    }

    return 0;
}