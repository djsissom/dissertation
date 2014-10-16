#!/usr/bin/env bash

minsnap=0
maxsnap=61

minbox=1
maxbox=3

for ((i=$minbox; i<=$maxbox; i++)); do
  if [ ! -e ../box$i/crossmatch ]; then
    mkdir -v ../box$i/crossmatch
  fi

  cp -v run_crossmatch.pbs ../box$i/.

  for ((snap=$minsnap; snap<=$maxsnap; snap++)); do
    if [ $snap -lt 10 ]; then
      j=00$snap
    elif [ $snap -lt 100 ]; then
      j=0$snap
    fi

    if [ ! -e ../box$i/crossmatch/snap$j ]; then
      mkdir -v ../box$i/crossmatch/snap$j
    fi

    cp -v -r crossmatch_proto/* ../box$i/crossmatch/snap$j/.

    echo "OUTPUT_DIR        /home/sissomdj/projects/simulations/rockstar/box$i/crossmatch/snap$j" >> ../box$i/crossmatch/snap$j/rockstar_2lpt.param
    echo "FIRST_GROUPDIR    /home/sissomdj/projects/simulations/rockstar/box$i/2lpt/snap$j/halos" >> ../box$i/crossmatch/snap$j/rockstar_2lpt.param
    echo "SECOND_GROUPDIR   /home/sissomdj/projects/simulations/rockstar/box$i/za/snap$j/halos" >> ../box$i/crossmatch/snap$j/rockstar_2lpt.param

    echo "OUTPUT_DIR        /home/sissomdj/projects/simulations/rockstar/box$i/crossmatch/snap$j" >> ../box$i/crossmatch/snap$j/rockstar_za.param
    echo "FIRST_GROUPDIR    /home/sissomdj/projects/simulations/rockstar/box$i/za/snap$j/halos" >> ../box$i/crossmatch/snap$j/rockstar_za.param
    echo "SECOND_GROUPDIR   /home/sissomdj/projects/simulations/rockstar/box$i/2lpt/snap$j/halos" >> ../box$i/crossmatch/snap$j/rockstar_za.param

  done

done


