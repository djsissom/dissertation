#!/bin/bash

echo 'running finish_bgc2...'
~/projects/programs/nbody/rockstar/Rockstar-0.99.9/util/finish_bgc2 -c onenode.cfg -s 0

echo 'running bgc2_to_ascii...'
~/projects/programs/nbody/rockstar/Rockstar-0.99.9/util/bgc2_to_ascii -c onenode.cfg -s 0 > halos/all_halos.bgc2.ascii

echo 'running find_parents...'
~/projects/programs/nbody/rockstar/Rockstar-0.99.9/util/find_parents halos/out_0.list 10.0 > halos/out_0.list.parents

echo 'finished'
