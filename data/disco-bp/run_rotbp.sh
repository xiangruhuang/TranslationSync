#!/bin/sh
export LD_LIBRARY_PATH=./corona-1.0.2/src/.libs/:$LD_LIBRARY_PATH

ROTBP_LABELSPACE='11 5' # MUST be odd

# Run bp-rot on the Acropolis dataset.  The energy and labels are printed each iteration.
bin/bp-rot --nopreprop --threads 10 --labelspace ${ROTBP_LABELSPACE} --iters 100 \
    --noconf --twoconf --vanish data/acropolis/tilt_twist.txt --geoplanar data/acropolis/geoplanar.geotags.txt \
    data/acropolis/pairs.txt
