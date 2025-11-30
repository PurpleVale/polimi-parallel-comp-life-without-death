# Project

Life without Death challenge for Parallel Computing Course a.y. 25/26 Polimi

# Authors


| Name        | Surname    | Person Code |
| ----------- | ---------- | ----------- |
| Vale        | Turco      | 10809855    |
| Alessia     | Rigoni     | 10859832    |
| Giulio Enzo | Donninelli | 10823453    |

# Code Profiling

The profiling of the code was done with a custom script `compile_and_log.sh`
The results can be found here: [Profiling on GoogleSheet](https://docs.google.com/spreadsheets/d/1q1YAxiXo2tiVPbItVfo0__0ttwRlqZgJAdyFnSGF6EM/edit?usp=sharing)

# Project Description

<img src="https://upload.wikimedia.org/wikipedia/commons/2/23/Without_death.gif" align="right" width="200px">

Realize a parallel implementation of **[Life without Death](https://en.wikipedia.org/wiki/Life_without_Death)** in OpenMP.
<br>
This is a spin on the rules of the well-known [Game of Life](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life) automata by John Conway!

<br>

For the uninitiated about classic automata:
<br>
You have a grid of cells of known $width$ and $height$.
<br>
Each cell can be either dead or alive.
<br>
The simulation updates in global time steps.
<br>
<img src="https://upload.wikimedia.org/wikipedia/commons/4/4d/Moore_neighborhood_with_cardinal_directions.svg" align="right" width="200px">
Each cell has a neighborhood (called "Moore neighborhood" for us, see figure on the right) constituted by the eight cells up-down-right-left of itself, plus those diagonally adjacent.
<br>
At every step, a cell:

- if it is dead, be born if it is surrounded by a number $b \in B$ of alive cells, otherwise, remain dead
- if it is alive, stay alive if it is surrounded by a number $s \in S$ of alive cells, otherwise die

The sets $B, S \subseteq \{0, \dots 8\}$ define the rules of the game.
<br>
The classic game of life has rules $B = \{3\}, S = \{2, 3\}$, typically written as B3/S23 in the family of [life-like automata](https://en.wikipedia.org/wiki/Life-like_cellular_automaton).

<br>

Our automata variation:
<br>
Rules are much closer to life without death than to game of life.
<br>
In particular $B = \{3, 6, 8\}, S = \{0, 1, 2, 3, 4, 5, 6, 7, 8\}$, or B368/S012345678. Thus cells can never die.
<br>
The grid "wraps-around" as in it behaves like a toroid. The neighbors of the cells on the edges of the grid are those on the other side w.r.t. their current border. In practice, an arbitrary neighbor with coordinates $(x, y)$ always wraps back to $(x \text{ mod } width, y \text{ mod } height)$.

<br>

Additional rules:
<br>
Each cell is assigned to a color encoded as a float in the range [0, 1] representing "hue" in the [HSV](https://it.wikipedia.org/wiki/Hue_Saturation_Brightness) format.
<br>
Whenever a cell is born, its color is permanently set for the length of its life to the average of all its "parent cells" colors (see the definition of "hue average" in the code).
<br><!--If a cell ever dies (the poor guy had zero neighbors, I feel sorry for him), its color is reset. Not like this could ever happen.
<br>-->
The simulation starts with randomly generated groups of alive cells, each assigned to a random color (see initialization in the code).
<br>
The simulation ends once two consecutive steps results in the same identical grid state, in other words, when a step does not cause any cell to change status.

<br>

Your task is to parallelize in OpenMP, to the best of your abilities, the simulation.
<br>
You are provided, through CLI arguments, with the grid's $width$ and $height$, and a seed that is used to determine a random starting grid configuration.
<br>
You must return the simulation's final, converged, grid state, complete for each cell with its status and color.
<br>
How you perform everything in between is up to you (and OpenMP).

<br>

*Note: you have [here](https://colab.research.google.com/drive/1c1nw8XrW7diQfiQgge58TxzOpSkNerDK?usp=sharing) a python script animating our particular version of the game. Before starting, I suggest you take a look at it, just to grasp on where and when meaningful computation actually occurs.*

*Note: there are notorious ways to speedup the simulation of life-like automata, most notably the [Hashlife](https://en.wikipedia.org/wiki/Hashlife) algorithm and its implementation in [Golly](https://en.wikipedia.org/wiki/Golly_(program)). These are beyond the scope of the course and you are not expected to implement similar techniques. However, their ideas very good readig material if you want to deepen your repertoir of algorithms, and may serve as inspiration for your solution to this challenge. To further discourage any wasted efforths, the specific life-without-death-like set of rules you are hereby tasked with, is not necessary the best target for optimizations like "memoization", as used in the above algorithms. Still, I am reporting here these tools and ideas to pick on your curiosity ^-^ !*

---

A reference sequential implementation is given. Implement your parallel version alongside it in `simulate_parallel`. Rely on the provided sequential version to check the correctness of results.

To get a general performance metric rely on `omp_get_wtime`. The most significant measurement is how much time your code takes to reach the final state of the simulation.

Note that your submissions will be profiled on 16 logical cores (8 physical, with hyperthreading) of an AMD EPYC 7453 @ 2.75GHz with 32GB of RAM.
You are not expected to own a machine with similar specifications, but it is recommended that you test your code, especially its scalability, on your local machine. Doing so will provide you with far more realistic metrics than Colab's two little cores.

To ensure the scalability of your code on the machine it will be tested on, you can assume that the number of threads to use will always be determined through the `OMP_NUM_THREADS` environment variable, which will be set to 16. However, if you want to, you are still entitled to override such behaviour with the `num_threads` clause or other by calling `omp_set_num_threads` after acquiring the available number of CPUs at runtime.

Furthermore, take note that when profiling your code, widely arbitrary (and large) grid sizes and seeds will be used, but none of the initialization logic as you see it now will be altered.

Step one is beating the reference implementation, that should be easy, then you can use all tricks in the book to push it further.
Anything goes, but if you use "exotic" tricks we want an explanation.
In fact, before submitting your work, be sure to fill out the [report](#report) with brief insights of what you did.
