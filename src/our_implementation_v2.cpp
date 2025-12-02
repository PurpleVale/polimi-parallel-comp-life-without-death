#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <ranges>

#define MAX_STEPS 1000000
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

#define INDEX(x,y) (y*W + x)
#define WRAP(x,y) INDEX((x+W)%W,(y+H)%H)
#define INLINE_KER_SUM(vec) vec[WRAP(x-1,y-1)]+vec[WRAP(x,y-1)]+vec[WRAP(x+1,y-1)]+\
                            vec[WRAP(x-1,y)]+                   vec[WRAP(x+1,y)]+\
                            vec[WRAP(x-1,y+1)]+vec[WRAP(x,y+1)]+vec[WRAP(x+1,y+1)]

inline bool birth_rule(int n);
float hue_average(const float *h, int n);

size_t part(std::vector<size_t> vec,  size_t high, size_t low) {
    const size_t pivot = vec[high];
    size_t moving_pointer = (low-1);

    for (auto j = low; j<= high-1; ++j ) {
        if (vec[j] <= pivot) {
            moving_pointer++;
            std::swap(vec.at(moving_pointer), vec.at(j));
        }
    }
    std::swap(vec.at(moving_pointer+1), vec.at(high));
    return moving_pointer+1;
}

void qsort(std::vector<size_t> vec, size_t high, size_t low, size_t depth = 0) {
    if (high>=low) return;
    #pragma omp single
    {
        size_t pivot  = part(vec, high, low);
        #pragma omp task shared(vec) if(depth<3)
        qsort(vec,pivot-1,low,depth+1);
        #pragma omp task shared(vec) if(depth<3)
        qsort(vec,high,pivot+1,depth+1);
        #pragma omp taskwait
    }
}

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
    std::vector<unsigned char> alive; // size: W*H
    std::vector<float> hue; // size: W*H
    std::vector<float> new_hue;
    std::vector<unsigned char> new_alive;
    std::vector<size_t> new_borders;
    std::vector<size_t> borders;
    // constructor: allocates and owns two dynamic arrays via std::vector
    Grid(int w, int h) : W(w), H(h),
                alive(w * h, 0), hue(w * h, 0.0f),new_hue(h*w,0.0f),new_alive(w*h,0)
    {};

    void step() {
        std::vector<size_t> new_borders_local;
        std::vector<size_t> borders_local(W*H);
        #pragma omp parallel private(new_borders_local) shared(borders_local)
        {
            new_borders_local.reserve((W*H/omp_get_num_teams()));
            #pragma omp for schedule(dynamic)
            for(const auto& index : borders) {
                new_alive[index] = 0;
                if (!alive[index]) {
                    const size_t x = index % W;
                    const size_t y = index / W;
                    const size_t NW_N   = WRAP(x-1,y-1);
                    const size_t N_N    = WRAP(x+0,y-1);
                    const size_t NE_N   = WRAP(x+1,y-1);
                    const size_t W_N    = WRAP(x-1,y+0);
                    const size_t E_N    = WRAP(x+1,y+0);
                    const size_t SW_N   = WRAP(x-1,y+1);
                    const size_t S_N    = WRAP(x+0,y+1);
                    const size_t SE_N   = WRAP(x+1,y+1);

                    int neighbours = 0;
                    float par_hues[8];

                    if (alive[NW_N])  par_hues[neighbours++] = hue[NW_N];
                    if (alive[N_N])   par_hues[neighbours++] = hue[N_N];
                    if (alive[NE_N])   par_hues[neighbours++] = hue[NE_N];

                    //lateral neighbours
                    if (alive[W_N])   par_hues[neighbours++] = hue[W_N];
                    if (alive[E_N])   par_hues[neighbours++] = hue[E_N];

                    //south neighbours
                    if (alive[SW_N])   par_hues[neighbours++] = hue[SW_N];
                    if (alive[S_N])   par_hues[neighbours++] = hue[S_N];
                    if (alive[SE_N])   par_hues[neighbours++] = hue[SE_N];

                    if (birth_rule(neighbours)) {
                        new_alive[index] = 1;
                        new_hue[index] = hue_average(par_hues,neighbours);

                        if (!alive[NW_N])  new_borders_local.push_back(NW_N);
                        if (!alive[N_N])   new_borders_local.push_back(N_N);
                        if (!alive[NE_N])   new_borders_local.push_back(NE_N);

                        //lateral neighbours
                        if (!alive[W_N])   new_borders_local.push_back(W_N);
                        if (!alive[E_N])   new_borders_local.push_back(E_N);

                        //south neighbours
                        if (!alive[SW_N])   new_borders_local.push_back(SW_N);
                        if (!alive[S_N])   new_borders_local.push_back(S_N);
                        if (!alive[SE_N])   new_borders_local.push_back(SE_N);
                    }
                }
            }

            #pragma omp parallel for schedule(dynamic)
            for(const auto& index : borders) {
                if (new_alive[index] == 1) {
                    alive[index] = 1;
                    hue[index] = new_hue[index];
                }
            }

            #pragma omp single
            {
                borders.clear();
            }
            #pragma omp critical
            {
                borders.insert(borders.begin(), new_borders_local.begin(), new_borders_local.end());
            }

            #pragma omp barrier

            if (!borders.empty()) {
                #pragma omp single
                std::sort(borders.begin(), borders.end());

                #pragma omp for schedule(dynamic)
                for(size_t i = 0; i < borders.size()-1; ++i) {
                    if (borders[i]!=borders[i+1]) {
                        #pragma omp critical
                        borders_local.push_back(borders[i]);
                    }
                }
                #pragma omp single
                if (borders_local[borders_local.size()-1]!=borders[borders.size()-1]) {
                    borders_local.push_back(borders[borders.size()-1]);
                }
                #pragma omp single
                std::swap(borders, borders_local);
            }
        }
    }

    int evolve() {
        int k = 0;
        while (k++ < MAX_STEPS && !borders.empty()) {
            step();
        }
        return k;
    }

    void initialize() {
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                if (!alive[INDEX(x,y)] && birth_rule(INLINE_KER_SUM(alive))) {
                    borders.push_back(INDEX(x,y));
                }
            }
        }
    }

};

// =================================================
// === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

// game rules B368 / S012345678:
// - a dead cell is born if alive neighbor count is 3, 6, or 8
// - a live cell survives with ANY count 0..8
inline bool birth_rule(int n) {
  return (n == 3 || n == 6 || n == 8);
}

inline bool survive_rule(int n) {
  return (n >= 0 && n <= 8);
}

// HELPER FUNCTION [YOU CAN USE IT AS-IS]
// hue averaging using circular mean around the unit circle
float hue_average(const float *h, int n) {
  float sumx = 0.0f, sumy = 0.0f;
  for (int i = 0; i < n; i++) {
    float angle = h[i] * 2.0f * M_PI;
    sumx += cosf(angle);
    sumy += sinf(angle);
  }
  float angle = atan2f(sumy, sumx);
  if (angle < 0) angle += 2.0f * M_PI;
  return angle / (2.0f * M_PI);
}

// HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
// random hue selection
float pick_group_hue(int group_id) {
  float base = float(group_id) / NUM_GROUPS;
  float h = base + ((float) rand() / RAND_MAX) * 2 * COLOR_VARIATION - COLOR_VARIATION;
  if (h < 0) h += 1.0f;
  if (h >= 1) h -= 1.0f;
  return h;
}

// HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
// grid initialization with a few clusters
void initialize_grid(Grid &g) {
    int W = g.W, H = g.H;

    int min_dim = (W < H ? W : H);
    int r = int(min_dim * CLUSTER_RADIUS_FACTOR);
    if (r < CLUSTER_MIN_RADIUS) r = CLUSTER_MIN_RADIUS;

    std::vector<float> group_hue(NUM_GROUPS);
    for (int i = 0; i < NUM_GROUPS; i++)
        group_hue[i] = pick_group_hue(i);

    std::vector<int> cx(NUM_GROUPS), cy(NUM_GROUPS);
    for (int i = 0; i < NUM_GROUPS; i++) {
        cx[i] = rand() % W;
        cy[i] = rand() % H;
    }

    for (int g_id = 0; g_id < NUM_GROUPS; g_id++) {
        float base_h = group_hue[g_id];
        int gx = cx[g_id], gy = cy[g_id];

        for (int dx = -r; dx <= r; dx++) {
            for (int dy = -r; dy <= r; dy++) {
                if (dy*dy + dx*dx <= r*r) {
                    int x = (gx + dx + 2*W) % W, y = (gy + dy + 2*H) % H;
                    if (((float)rand()/RAND_MAX) < CLUSTER_FILL_DENSITY) {
                        int idx = y*W + x;
                        g.alive[idx] = 1;
                        g.hue[idx] = base_h;
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
    bool ret = true;
    for (int i = 0; i < N; i++) {
        if (a.alive[i] != b.alive[i]) {ret = false; std::cerr << "Error alive at" << i << std::endl;}
        if (a.alive[i]) {
            float ha = a.hue[i];
            float hb = b.hue[i];
            // compare hue with a tolerance on floats (due to associativity)
            if (fabs(ha - hb) > 1e-4f) {ret = false; std::cerr << "Error hue at" << i << std::endl;}
        }
    }
    return ret;
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
int evolve_step(const Grid &cur, Grid &next) {
  int W = cur.W, H = cur.H;
  int changes = 0;

  // iterate over all cells
  for (int x = 0; x < W; x++) {
    for (int y = 0; y < H; y++) {
      int idx = y * W + x;
      unsigned char alive = cur.alive[idx];

      int alive_neighbors = 0;
      float parent_hues[8];

      // count alive moore neighbors and collect their hues
      for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
          if (dx == 0 && dy == 0) continue;
          // wrap around the grid (torus)
          int xx = (x + dx + W) % W;
          int yy = (y + dy + H) % H;
          int nidx = yy * W + xx;
          if (cur.alive[nidx]) {
            parent_hues[alive_neighbors] = cur.hue[nidx];
            alive_neighbors++;
          }
        }
      }

      unsigned char new_alive = alive;

      if (!alive) {
        if (birth_rule(alive_neighbors)) {
          new_alive = 1;
          next.hue[idx] = hue_average(parent_hues, alive_neighbors);
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
    long step;
  for (step = 0; step < MAX_STEPS; step++) {
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
// as of now this is a stub: it just calls the sequential version...
int simulate_parallel(Grid &g) {
    return g.evolve();
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

  Grid gs(W, H), gp(W, H);
  initialize_grid(gs);

  // copy the initial state
  gp = gs;

  // === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
  // =================================================

  // if you need to initialize additional things, do that here!

    gp.initialize();

  // =================================================
  // === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

  // sequential

  // parallel
    std::cout << "start par" << std::endl;
  double t3 = omp_get_wtime();
  simulate_parallel(gp);
  double t4 = omp_get_wtime();

    std::cout << "start seq" << std::endl;
#if !DISABLE_SEQUENTIAL
    double t1 = omp_get_wtime();
    simulate_sequential(gs);
    double t2 = omp_get_wtime();
#endif

#if !DISABLE_SEQUENTIAL
    printf("Sequent. time: %.6f s\n", t2 - t1);
#endif
    printf("Parallel time: %.6f s\n", t4 - t3);

    std::cout << "checking" << std::endl;
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
        write_grid_to_file(gs, "base.lwd");
        printf("Results written to %s\n", (char *)argv[4]);
    }

  return 0;
}
