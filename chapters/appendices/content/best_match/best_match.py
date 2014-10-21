#!/usr/bin/env python

import sys
import getopt
import numpy as np


def main():
	# read in files
	print 'reading files...'
	with open(sys.argv[1]) as f:
		matches1 = f.readlines()
	with open(sys.argv[2]) as f:
		matches2 = f.readlines()
	print 'done reading files'

	header = matches1[2:6]
	header.insert(0, '# Best matches for bi-directional crossmatch\n')
	header.insert(1, '#\n')

	matches1 = matches1[7:]
	matches2 = matches2[7:]

	# convert to numpy arrays
	print 'converting to numpy arrays...'
	match_array1 = np.asarray([line.split() for line in matches1], dtype=int)
	match_array2 = np.asarray([line.split() for line in matches2], dtype=int)
	print 'done converting'

	# find matches that exist in both lists
	print 'finding matches...'
	mask = np.zeros(len(match_array1), dtype=bool)
	for i, line in enumerate(match_array1):
		id1 = line[id1_col]
		id2 = line[id2_col]
		tmp = (match_array2[:,id1_col] == id2)
		tmp = (match_array2[tmp,id2_col] == id1)
		mask[i] = tmp.any()
		if i % 1000 == 0:
			print "Finished line ", i

	print 'done matching'
	
	out_array = match_array1[mask]

	# write results
	print 'writing results...'
	with open(sys.argv[3], 'w') as f:
		f.writelines(("%s" % line for line in header))
		np.savetxt(f, out_array, fmt='%10d')

	print 'Finished.'


id1_col             =  4
npart1_col          =  5
id2_col             =  1
npart2_col          =  2
ncommon_col         =  6
hnum1_col           =  3
hnum2_col           =  0


if __name__ == '__main__':
  main()
