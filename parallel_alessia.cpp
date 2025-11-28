#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <omp.h>

#define MAX_STEPS 1000000 //* limite massimo di iterazioni
#define NUM_GROUPS 7
#define COLOR_VARIATION 0.08f
#define CLUSTER_FILL_DENSITY 0.55f
#define CLUSTER_RADIUS_FACTOR 0.08f
#define CLUSTER_MIN_RADIUS 3

//? === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
//? =================================================

//? if you need any additional standard library up to C11 or C++17, add it here

//? if you need any additional define, add it here

//? toggle this to save time when profiling
//? be wary: this also disables the results check
#define DISABLE_SEQUENTIAL false //! da modificare quando si farà profiling

//? if you happen to need extra information in the grid, you can modify this object, however,
//? only add, do not remove anything, as to preserve compatibility with the sequential version!
//
//? notes about the use of C++ here:
//? - Grid is a C++ object owning two heap-allocated arrays
//? - copying or assigning a Grid (Grid B = A; or B = A;) performs a deep copy of both vectors
//? - passing Grid by value also triggers this deep copy, thus always pass by reference
//? => this is different from C: the struct looks cheap, but copies are expensive

//TODO struttura dati fondamentali
//TODO invece di una matrice bidimensionale [x][y] usa un singolo vettore 1D (alive e hue) di dimensione W*H

//! NB l'accesso alla cella (x,y) si fa con l'indice idx = y * W + x

struct Grid {
  int W, H;
  //* contiene 0 (morta) o 1 (viva)
  std::vector<unsigned char> alive; //? size: W*H
  //* contiene il colore (float)
  std::vector<float> hue; //? size: W*H
  //? constructor: allocates and owns two dynamic arrays via std::vector
  Grid(int w, int h) : W(w), H(h), alive(w*h,0), hue(w*h,0.0f) {}
};

//? =================================================
//? === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

//? game rules B368 / S012345678:
//? - a dead cell is born if alive neighbor count is 3, 6, or 8
//? - a live cell survives with ANY count 0..8
inline bool birth_rule(int n) {
  return (n == 3 || n == 6 || n == 8);
}
inline bool survive_rule(int n) {
  return (n >= 0 && n <= 8);
} //* una volta viva, una cella non muore mai

//? HELPER FUNCTION [YOU CAN USE IT AS-IS]
//? hue averaging using circular mean around the unit circle

//TODO calcola il nuovo colore
//TODO è computazionalmente costosa perché usa funzioni trigonometriche
//TODO più celle nascono, più il carico computazionale aumenta

float hue_average(const float *h, int n) {
  float sumx = 0.0f, sumy = 0.0f; //* accumulatori per le coordinate cartesiane
  //TODO ciclo su tutti i colori dei vicini
  for (int i = 0; i < n; i++) {
    float angle = h[i] * 2.0f * M_PI; //* converte il valore [0,1] in radianti [0,2pi]
    sumx += cosf(angle); //* somma la componente X (coseno)
    sumy += sinf(angle); //* somma la componente Y (seno)
  }
  float angle = atan2f(sumy, sumx); //* converte il vettore risultante (X,Y) di nuovo in angolo
  if (angle < 0) angle += 2.0f * M_PI; //* normalizza l'angolo se è negativo
  return angle / (2.0f * M_PI); //* riporta l'angolo nel range [0,1]
}

//^ NOTE PER PARALLELIZZAZIONE
//^ poiché viene chiamata ogni volta che una cella nasce, le zone della griglia dove nascono molte celle richiederanno
//^ molto più tempo di calcolo rispetto alle zone stabili -> sbilanciamento del carico (Load Imbalance)

//? HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
//? random hue selection

//TODO genera un colore per un gruppo di celle iniziale

float pick_group_hue(int group_id) {
  float base = float(group_id) / NUM_GROUPS; //* dividiamo lo spettro dei colori in base al numero dei gruppi
  //TODO aggiunge un "rumore" casuale al colore base per non averli tutti identici
  float h = base + ((float)rand() / RAND_MAX)*2*COLOR_VARIATION - COLOR_VARIATION; //* crea un valore casuale tra [-COLOR_VARIATION, +COLOR_VARIATION]
  //* gestisce il "wrap-around" -> se esce da [0,1], lo riporta dentro
  if (h < 0) h += 1.0f; //* 1.1 diventa 0.1
  if (h >= 1) h -= 1.0f; //* -0.1 diventa 0.9
  return h;
}

//? HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
//? grid initialization with a few clusters

//TODO rimpie la griglia vuota con dei "blob" (cluster) circolari di celle vive

void initialize_grid(Grid &g) {
  int W = g.W, H = g.H;

  //* calcola il raggio dei cerchi in base alla dimensione minima della griglia
  int min_dim = (W < H ? W : H);
  int r = int(min_dim * CLUSTER_RADIUS_FACTOR);
  if (r < CLUSTER_MIN_RADIUS) r = CLUSTER_MIN_RADIUS; //* assicura dimensione minima

  //* pre-calcola i colori per ogni gruppo
  std::vector<float> group_hue(NUM_GROUPS);
  for (int i = 0; i < NUM_GROUPS; i++)
    group_hue[i] = pick_group_hue(i);

  //* sceglie posizioni casuali (centri) per ogni gruppo
  std::vector<int> cx(NUM_GROUPS), cy(NUM_GROUPS);
  for (int i = 0; i < NUM_GROUPS; i++) {
    cx[i] = rand() % W;
    cy[i] = rand() % H;
  }

  //* ciclo principale per creare i gruppi
  for (int g_id = 0; g_id < NUM_GROUPS; g_id++) {
    float base_h = group_hue[g_id];
    int gx = cx[g_id], gy = cy[g_id];

    //* itera in un quadrato attorno al centro del gruppo
    for (int dx = -r; dx <= r; dx++) {
      for (int dy = -r; dy <= r; dy++) {
        //* controlla se il punto è dentro al cerchio (tramite distanza euclidea)
        if (dy*dy + dx*dx <= r*r) {
            //* calcola coordinate gestendo il toroide
          int x = (gx + dx + 2*W) % W, y = (gy + dy + 2*H) % H; //* +2*W garantisce che il numero sia positivo prima del modulo
          
          //* riempie la cella in modo probabilistico
          if (((float)rand()/RAND_MAX) < CLUSTER_FILL_DENSITY) {
            int idx = y*W + x; //* indice lineare
            g.alive[idx] = 1; //* rende viva la cella
            g.hue[idx] = base_h; //* assegna il colore
          }
        }
      }
    }
  }
}

//? HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
//? compares two grids

//TODO confronta la griglia calcolata in sequenziale con quella calcolata in parallelo

bool compare_grids(const Grid &a, const Grid &b) {
  int N = a.W * a.H;
  for (int i = 0; i < N; i++) { //* scansione tutte le celle
    //! se lo stato (viva/morta) è diverso, errore immediato
    if (a.alive[i] != b.alive[i]) return false;
    //! se sono entrambe vive, controlla il colore
    if (a.alive[i]) {
      float ha = a.hue[i];
      float hb = b.hue[i];
      //? compare hue with a tolerance on floats (due to associativity)
      if (fabs(ha - hb) > 1e-4f) return false; //* confrontiamo i float con una tolleranza
    }
  }
  return true;
}

//? HELPER FUNCTION [YOU DON'T NEED TO USE THIS]
//? write grid to file for later visualization :)

//TODO salva lo stato della griglia in un file binario per poterlo visualizzare

void write_grid_to_file(const Grid &g, const char *filename) {
  FILE *f = fopen(filename, "wb"); //* apre il file in scrittura binaria
  if (!f) {
    fprintf(stderr, "ERROR: cannot open %s for writing\n", filename);
    return;
  }
  uint32_t W = g.W;
  uint32_t H = g.H;
  //* scrive le dimensioni (Width e Height) come interi a 32 bit
  fwrite(&W, sizeof(uint32_t), 1, f);
  fwrite(&H, sizeof(uint32_t), 1, f);

  //* scrive tutto l'array di byte "alive" in blocco
  for (size_t i = 0; i < W * H; i++)
    fwrite(&g.alive[i], sizeof(uint8_t), 1, f);

  //* scrive tutto l'array di float "hue" in blocco
  for (size_t i = 0; i < W * H; i++)
    fwrite(&g.hue[i], sizeof(float), 1, f);
  fclose(f);
}

//TODO questa funzione calcola lo stato successivo (next) basandosi su quello corrente (cur)

//? advance the simulation by one step; return the count of changed cells
int evolve_step(const Grid &cur, Grid &next) {
  int W = cur.W, H = cur.H;
  int changes = 0;

  //? iterate over all cells
  for (int x = 0; x < W; x++) {
    for (int y = 0; y < H; y++) {
      //! il fatto che abbiamo x esterno e y interno e l'indice calcolato come y*W+x,
      //! abbiamo che l'indice idx fa salti di grandezza -> inefficiente per la memoria 
      int idx = y*W + x; //* calcolo indice
      unsigned char alive = cur.alive[idx];

    int alive_neighbors = 0;
    float parent_hues[8];

    //? count alive moore neighbors and collect their hues
    //* questi cicli scansionano le 9 celle attorno
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        if (dx == 0 && dy == 0) continue;
        //? wrap around the grid (torus)
        //TODO questo modulo gestisce i bordi -> se esci a destra, rientri a sinistra
        //* questo viene gestito tramite il modulo in quanto, se siamo sul bordo destrp (x=W-1) e guardiamo a sinistra (dx=+1)
        //* la somma è W -> facendo W % W otteniamo 0 (il bordo sinistro)
        int xx = (x + dx + W) % W;
        int yy = (y + dy + H) % H;
        int nidx = yy*W + xx;
        //TODO se il vicino è vivo, lo contiamo e salviamo il suo colore
        if (cur.alive[nidx]) {
          parent_hues[alive_neighbors] = cur.hue[nidx];
          alive_neighbors++;
        }
      }
    }

    unsigned char new_alive = alive;

    if (!alive) { //* se la cella è morta
      //TODO deve nascere se ha 3, 6, 8 vicini
      if (birth_rule(alive_neighbors)) {
        new_alive = 1;
        //TODO media dei colori trigonometrica
        next.hue[idx] = hue_average(parent_hues, alive_neighbors); //* chiamata solo quando una cella nasce
      }
    } else { //* se la cella è viva
      //? you are free to skip this check, here it was kept for sake of completeness
      //TODO sopravvive sempre (0,...,8 vicini)
      if (survive_rule(alive_neighbors)) {
        new_alive = 1;
        next.hue[idx] = cur.hue[idx]; //* mantiene il vecchio colore
      } else {
        new_alive = 0;
      }
    }

    //TODO scrive il risultato nella griglia next
    next.alive[idx] = new_alive;
    //TODO se lo stato è cambiato (da morta a viva), incrementa il contatore
    if (new_alive != alive)
      changes++;
    }
  }
  return changes; //* restituisce changes, il numero di celle che hanno cambiato stato
  //* se ritorna 0 vuole dire che la simulazione è finita
}

//? sequential simulation entry point
void simulate_sequential(Grid &g) {
  //? it's hard to perform updates in place, we ping-pong between grid copies
  //TODO creiamo una seconda griglia "tmp" vuota con
  //TODO g -> stato corrente (sola lettura per il calcolo)
  //TODO tmp -> stato prossimo (sola scrittura)
  Grid tmp(g.W, g.H);

  for (long step = 0; step < MAX_STEPS; step++) {
    //? one step at a time
    //? note: fully overwrites tmp with the next state
    //TODO questa funzione legge da "g" e scrive su "tmp"
    //TODO restituisce quanti cambi ci sono stati
    int changes = evolve_step(g, tmp);

    //? swap g <-> tmp
    //TODO cambiamo i puntatori interni dei vettori
    //TODO in questo modo "g" punta ai dati nuovi appena calcolati e "tmp" ai dati vecchi
    g.alive.swap(tmp.alive);
    g.hue.swap(tmp.hue);

    //! se nessuno è nato o morto la convergenza è stabile
    if (changes == 0) break;
  }
}

//? === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
//? =================================================

int evolve_step_par(const Grid &cur, Grid &next) {
  int W = cur.W, H = cur.H;
  int changes = 0;

  //? iterate over all cells
  #pragma omp for collapse(2) schedule(dynamic)
  for (int y = 0; y < H; y++) {
    for (int x = 0; x < W; x++) {
      int idx = y*W + x;
      unsigned char alive = cur.alive[idx];

    int alive_neighbors = 0;
    float parent_hues[8];

    //? count alive moore neighbors and collect their hues
    //* questi cicli scansionano le 9 celle attorno
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        if (dx == 0 && dy == 0) continue;
        //? wrap around the grid (torus)
        //TODO questo modulo gestisce i bordi -> se esci a destra, rientri a sinistra
        //* questo viene gestito tramite il modulo in quanto, se siamo sul bordo destrp (x=W-1) e guardiamo a sinistra (dx=+1)
        //* la somma è W -> facendo W % W otteniamo 0 (il bordo sinistro)
        int xx = (x + dx + W) % W;
        int yy = (y + dy + H) % H;
        int nidx = yy*W + xx;
        //TODO se il vicino è vivo, lo contiamo e salviamo il suo colore
        if (cur.alive[nidx]) {
          parent_hues[alive_neighbors] = cur.hue[nidx];
          alive_neighbors++;
        }
      }
    }

    unsigned char new_alive = alive;

    if (!alive) { //* se la cella è morta
      //TODO deve nascere se ha 3, 6, 8 vicini
      if (birth_rule(alive_neighbors)) {
        new_alive = 1;
        //TODO media dei colori trigonometrica
        next.hue[idx] = hue_average(parent_hues, alive_neighbors); //* chiamata solo quando una cella nasce
      }
    } else { //* se la cella è viva
      //? you are free to skip this check, here it was kept for sake of completeness
      //TODO sopravvive sempre (0,...,8 vicini)
      if (survive_rule(alive_neighbors)) {
        new_alive = 1;
        next.hue[idx] = cur.hue[idx]; //* mantiene il vecchio colore
      } else {
        new_alive = 0;
      }
    }

    //TODO scrive il risultato nella griglia next
    next.alive[idx] = new_alive;
    //TODO se lo stato è cambiato (da morta a viva), incrementa il contatore
    if (new_alive != alive)
      changes++;
    }
  }
  return changes; //* restituisce changes, il numero di celle che hanno cambiato stato
  //* se ritorna 0 vuole dire che la simulazione è finita
}

//? as of now this is a stub: it just calls the sequential version...
void simulate_parallel(Grid &g) {
  Grid tmp(g.W, g.H);
  int total_changes = 0;
  bool converged = false;

  #pragma omp parallel
  {
    for (long step = 0; step < MAX_STEPS; step++) {
      if (converged) break;

      #pragma omp single
      total_changes = 0; //* bisogna resettarla all'inizio di ogni step
      
      int step_changes = evolve_step_par(g, tmp);

      #pragma omp atomic
      total_changes += step_changes;

      #pragma omp barrier
      
      #pragma omp single
      {
        g.alive.swap(tmp.alive);
        g.hue.swap(tmp.hue);

        if (total_changes == 0) converged = true;
      }
    }
  }
}

//? =================================================
//? === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

int main(int argc, char **argv) {
    omp_set_nested(true); //? just in case, enable nested parallelism

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
    //? copy the initial state
    gp = gs;

    //? === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
    //? =================================================

    //? if you need to initialize additional things, do that here!

    //? =================================================
    //? === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

    //? sequential
    #if !DISABLE_SEQUENTIAL
    double t1 = omp_get_wtime();
    simulate_sequential(gs);
    double t2 = omp_get_wtime();
    #endif

    //? parallel
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

    //? === DO NOT CHANGE ANYTHING ABOVE THIS COMMENT ===
    //? =================================================

    //? if you need more logging, put it here!

    //? =================================================
    //? === DO NOT CHANGE ANYTHING BELOW THIS COMMENT ===

    if (argc == 5) {
      write_grid_to_file(gp, (char *)argv[4]);
      printf("Results written to %s\n", (char *)argv[4]);
    }

    return 0;
}