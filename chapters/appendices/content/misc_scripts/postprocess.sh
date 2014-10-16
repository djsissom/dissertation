#!/usr/bin/env bash

startdir=`pwd`

for snapdir in {2lpt,za}/*; do
  echo Working on $snapdir...
  cd $startdir/$snapdir

  ./postprocess

done

# - end of script
