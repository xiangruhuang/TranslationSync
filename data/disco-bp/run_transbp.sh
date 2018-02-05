#!/bin/sh
export LD_LIBRARY_PATH=./corona-1.0.2/src/.libs/:$LD_LIBRARY_PATH

# Run transbp on the Acropolis dataset.  The energy and labels are printed each iteration.
./bin/bp-trans --iters 15 --threads 10 --anticollapse 0.00001 1.0 --labelspace 101 0.2 --unweight 2 --dumpmessages --translate  data/acropolis/global.translations.txt data/acropolis/global.poses.txt 

