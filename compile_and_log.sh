#!/usr/bin/zsh

g++ src/our_implementation.cpp -O2 -pg -fopenmp -o bin/LwD
GRID=(100 180 500 700 1040 720 1920 1080)
for key val in "${(@kv)GRID}"; do
  for i in {1..10}; do
    echo $key $val
    bin/LwD $key $val 67 >> ./output/${key}x${val}.txt
  done
done