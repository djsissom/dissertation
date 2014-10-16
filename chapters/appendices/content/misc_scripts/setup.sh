#!/usr/bin/env bash

minsnap=0
maxsnap=61

minbox=1
maxbox=3

for ((i=$minbox; i<=$maxbox; i++)); do
  if [ ! -e ../box$i ]; then
    mkdir -v ../box$i
  fi
  if [ ! -e ../box$i/2lpt ]; then
    mkdir -v ../box$i/2lpt
  fi
  if [ ! -e ../box$i/za ]; then
    mkdir -v ../box$i/za
  fi
  if [ ! -e ../box$i/crossmatch ]; then
    mkdir -v ../box$i/crossmatch
  fi

  cp -v run_*.pbs ../box$i/.
  cp -v postprocess.sh ../box$i/.

  for ((snap=$minsnap; snap<=$maxsnap; snap++)); do
    if [ $snap -lt 10 ]; then
      j=00$snap
    elif [ $snap -lt 100 ]; then
      j=0$snap
    fi

    if [ ! -e ../box$i/2lpt/snap$j ]; then
      mkdir -v ../box$i/2lpt/snap$j
    fi
    if [ ! -e ../box$i/za/snap$j ]; then
      mkdir -v ../box$i/za/snap$j
    fi

    cp -v -r proto/* ../box$i/2lpt/snap$j/.
    cp -v -r proto/* ../box$i/za/snap$j/.

    ln -v -s ~/projects/data/2lpt/box$i/2lpt_512_z300_PM_$j ../box$i/2lpt/snap$j/particles/2lpt_512_z300_PM_$j
    ln -v -s ~/projects/data/za/box$i/za_512_z300_PM_$j ../box$i/za/snap$j/particles/za_512_z300_PM_$j

    echo /home/sissomdj/projects/simulations/rockstar/box$i/2lpt/snap$j/particles/2lpt_512_z300_PM_$j > ../box$i/2lpt/snap$j/particles/snapnames.lst
    echo /home/sissomdj/projects/simulations/rockstar/box$i/za/snap$j/particles/za_512_z300_PM_$j > ../box$i/za/snap$j/particles/snapnames.lst

    echo "BGC2_SNAPNAMES = \"/home/sissomdj/projects/simulations/rockstar/box$i/2lpt/snap$j/particles/snapnames.lst\"">> ../box$i/2lpt/snap$j/onenode.cfg
    echo "BGC2_SNAPNAMES = \"/home/sissomdj/projects/simulations/rockstar/box$i/za/snap$j/particles/snapnames.lst\"">> ../box$i/za/snap$j/onenode.cfg

    echo "FILENAME = \"2lpt_512_z300_PM_$j\"" >> ../box$i/2lpt/snap$j/onenode.cfg
    echo "FILENAME = \"za_512_z300_PM_$j\"" >> ../box$i/za/snap$j/onenode.cfg
  done

done


