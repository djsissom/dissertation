#!/usr/bin/env python

import sys
import numpy as np
from ipdb import set_trace

def main():
	# Read in particle files
	header, halos = read_files(sys.argv[1:], header_line = 3)

	if remove_nonfit_halos:
		print 'Removing NaNs...'
		halos = halos[np.isfinite(halos[:,c_lpt_col])]
		halos = halos[np.isfinite(halos[:,c_za_col])]

	if global_filter_halos:
		print 'Filtering data...'
		for col, val in zip(glob_lt_cols, glob_lt_vals):
			halos = halos[halos[:, col] <= val]
		for col, val in zip(glob_gt_cols, glob_gt_vals):
			halos = halos[halos[:, col] >= val]
		for col, val in zip(glob_eq_cols, glob_eq_vals):
			halos = halos[halos[:, col] == val]
		for col, val in zip(glob_ne_cols, glob_ne_vals):
			halos = halos[halos[:, col] != val]
	

	if sort_col != None:
		halos = sort_by_column(halos, sort_col)
	if (nhalos != None) or (nhalos != 0):
		halos = halos[:nhalos]
	#if (nhalos == 'perc25'):
	#	halos = halos[:len(halos)/10]
	if bad_halo_pairs != None:
		mask = np.arange(len(halos))
		mask = np.in1d(mask, bad_halo_pairs)
		mask = np.invert(mask)
		halos = halos[mask]

	c_rockstar_2lpt = halos[:, Rv1_col] / halos[:, Rs1_col]
	c_rockstar_za   = halos[:, Rv2_col] / halos[:, Rs2_col]
	if use_klypin:
		mask = (halos[:,4] < 100)
		print "changed %d halos" % (mask.sum())
		print "c_2lpt before ", c_rockstar_2lpt[mask][0]
		c_rockstar_2lpt[mask] = halos[mask, Rv1_col] / halos[mask, 79]
		print "c_2lpt klypin ", c_rockstar_2lpt[mask][0]
		mask = (halos[:,5] < 100)
		print "changed %d halos" % (mask.sum())
		print "c_za before ", c_rockstar_za[mask][0]
		c_rockstar_za[mask] = halos[mask, Rv2_col] / halos[mask, 80]
		print "c_za klypin ", c_rockstar_za[mask][0]
	c_diff_2lpt = 2.0 * (c_rockstar_2lpt - halos[:, c_lpt_col]) / (c_rockstar_2lpt + halos[:, c_lpt_col])
	c_diff_za   = 2.0 * (c_rockstar_za   - halos[:, c_za_col])  / (c_rockstar_za   + halos[:, c_za_col]) 
	#halos = np.column_stack((halos, c_rockstar_2lpt, c_rockstar_za, c_diff_2lpt, c_diff_za))
	#header.append('c_rockstar')
	#header.append('c_rockstar')
	#header.append('c_diff')
	#header.append('c_diff')

	c_diff_2lpt = c_diff_2lpt[np.isfinite(c_diff_2lpt)]
	c_diff_za = c_diff_za[np.isfinite(c_diff_za)]
	c_diff_tot = np.append(c_diff_2lpt, c_diff_za)

	c_diff_2lpt_frac = (np.abs(c_diff_za) <= cutoff_diff_frac).sum() / float(len(c_diff_2lpt))
	c_diff_za_frac = (np.abs(c_diff_za) <= cutoff_diff_frac).sum() / float(len(c_diff_za))
	c_diff_tot_frac = (np.abs(c_diff_tot) <= cutoff_diff_frac).sum() / float(len(c_diff_tot))

	with open(c_diff_file, 'w') as fd:
		fd.write("%g  %g  %g\n" % (c_diff_tot_frac, c_diff_za_frac, c_diff_2lpt_frac))

	print 'Finished snapshot.'


def read_files(files, header_line = None, comment_char = '#'):
	header = None
	data = None
	if type(files) == str:
		files = [files]

	if header_line != None:
		with open(files[0], 'r') as fd:
			for line in range(header_line):
				fd.readline()
			header = fd.readline()
		if header[0] != comment_char:
			print "Header must start with a '%s'" % comment_char
			sys.exit(4)
		header = header[1:]
		header = header.split()

	for file in files:
		print 'Reading file %s...' % (file)
		if data == None:
			data = np.genfromtxt(file, comments=comment_char)
		else:
			data = np.append(data, np.genfromtxt(file, comments=comment_char), axis=0)

	print 'Finished reading files.'
	if header_line == None:
		return data
	else:
		return header, data


def sort_by_column(halos, col):
	print 'Sorting halos...'
	mask = np.argsort(halos[:, col])
	mask = mask[::-1]
	halos = halos[mask]
	return halos



remove_nonfit_halos = False
global_filter_halos = True
use_klypin = False

nhalos = 100
#nhalos = 'perc25'
#sort_col = None
sort_col = 9

cutoff_diff_frac = 0.2


Rv1_col = 53
Rv2_col = 54
Rs1_col = 55
Rs2_col = 56

c_lpt_col = 17
c_za_col  = 18


lt_cols = [17, 18]
lt_vals = [100.0, 100.0]

gt_cols = [17, 18, 31, 32]
gt_vals = [1.0, 1.0, 0.0, 0.0]

eq_cols = []
eq_vals = []

ne_cols = []
ne_vals = []


# global filters
glob_lt_cols = []
glob_lt_vals = []

glob_gt_cols = [4, 5]
glob_gt_vals = [100, 100]

glob_eq_cols = [109, 110]
glob_eq_vals = [-1, -1]

glob_ne_cols = []
glob_ne_vals = []

bad_halo_pairs = None

c_diff_file = 'stats/c_diff.dat'



if __name__ == '__main__':
	main()

