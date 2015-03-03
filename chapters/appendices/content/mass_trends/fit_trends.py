#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import gamma as Gamma
from scipy.special import psi as digamma
from ipdb import set_trace


def main():
	#for filenum, file in enumerate(sys.argv[1:]):
	if (len(sys.argv[1:]) == 3):
		data1 = read_files(sys.argv[1], header_line = None)
		data2 = read_files(sys.argv[2], header_line = None)
		rsnap_data = read_files(sys.argv[3], header_line = None)
	else:
		print 'need 3 files'
		sys.exit(15)
	
	if fit_trend:
		with open(statsfile, 'w') as fd:
			fd.write("#plot slope slope_err intercept intercept_err\n")

	if minsnap > 0:
		#for data in data1, data2, data3:
		#	data = data[data[:,0] >= minsnap]
		data1 = data1[data1[:,0] >= minsnap]
		data2 = data2[data2[:,0] >= minsnap]

	z = 1.0 / rsnap_data[:,1] - 1.0
	if (len(data1) == len(data2)):
		z = z[-len(data1):]
	else:
		sys.exit(16)
	
	data1 = np.column_stack((data1, z))
	data2 = np.column_stack((data2, z))

	for data in [data1, data2]:
		data[:,slope_err_col] = np.sqrt(data[:,slope_err_col])       # var to stdev
		data[:,intercept_err_col] = np.sqrt(data[:,intercept_err_col])       # var to stdev

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# plots                                                                   #
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	for (data, ylabel, color, label, name) in zip([data1, data2], ylabels1, colors, labels1, names):
		print "Making %s plot..." % (name)
		fig = plt.figure(figsize=(9.0, 6.0))
		ax = fig.add_subplot(111)

		ax = make_plot(ax, data[:,z_col] - offset, data[:,slope_col], err = data[:,slope_err_col], color = 'blue', marker='o', label=label)

		if plot_intercept:
			if separate_axes:
				ax = ax.twinx()
			ax = make_plot(ax, data[:,z_col] + offset, data[:,intercept_col], err = data[:,intercept_err_col], color = 'red', marker='o', label=label)

		if fit_trend:
			ax, slope, slope_err, intercept, intercept_err = add_fit(ax, data[:,z_col], data[:,slope_col], err=data[:,slope_err_col], color='red')
			save_fits(statsfile, name, slope, np.sqrt(slope_err), intercept, np.sqrt(intercept_err))

		#ax.legend(loc='lower right')
		ax.set_xlim(z[0] + 1.0, z[-1] - 1.0)
		#ax.invert_xaxis()

		ax.set_xlabel(xlabel, fontsize='x-large')
		ax.set_ylabel(ylabel, fontsize='x-large')

		fig.tight_layout()
		fig.savefig(plot_base + 'mean_stdev_' + name + plot_ext, bbox_inches='tight')

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# make skew and kurtosis plots                                            #
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	'''
	for (data, ylabel_kurt, ylabel_skew, color, name, ylim_low1, ylim_high1, ylim_low2, ylim_high2) in zip([data1, data2, data3], ylabels2_kurt, ylabels2_skew, colors, names, [-10.0, -10.0, -1.0], [20.0, 20.0, 1.5], [-0.2, -1.5, -0.4], [0.5, 3.5, 0.1]):
		print "Making %s plot..." % (name)
		fig = plt.figure(figsize=(9.0, 6.0))
		ax = fig.add_subplot(111)

		ax = make_plot(ax, data[:,z_col] - offset, data[:,kurt_col], err = data[:,kurt_err_col], color = 'red', marker='o', linestyle=':', label='Kurtosis')
		legend_lines1, legend_labels1 = ax.get_legend_handles_labels()

		ax.set_xlabel(xlabel, fontsize='x-large')
		ax.set_ylabel(ylabel_kurt, fontsize='x-large')
		ax.set_ylim(ylim_low1, ylim_high1)

		if separate_axes:
			ax = ax.twinx()
		ax = make_plot(ax, data[:,z_col] + offset, data[:,skew_col], err = data[:,skew_err_col], color = 'blue', marker='o', linestyle=':', label='Skew')
		legend_lines2, legend_labels2 = ax.get_legend_handles_labels()

		ax.set_ylabel(ylabel_skew, fontsize='x-large')
		ax.legend(legend_lines1 + legend_lines2, legend_labels1 + legend_labels2, loc='lower right')
		ax.set_xlim(z[0] + 1.0, z[-1] - 1.0)
		ax.set_ylim(ylim_low2, ylim_high2)

		fig.tight_layout()
		fig.savefig(plot_base + 'skew_kurtosis_' + name + plot_ext, bbox_inches='tight')

	'''
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	print 'Finished all plots.'


def make_plot(ax, x, y, err=None, color='black', marker='None', linestyle='None', label=None):
	if err == None:
		if label == None:
			ax.plot(x, y, color=color, marker=marker, linestyle=linestyle)
		else:
			ax.plot(x, y, color=color, marker=marker, linestyle=linestyle, label=label)
	else:
		if label == None:
			ax.errorbar(x, y, yerr=err, color=color, marker=marker, linestyle=linestyle)
		else:
			ax.errorbar(x, y, yerr=err, color=color, marker=marker, linestyle=linestyle, label=label)
	return ax



def add_fit(ax, x, y, err=None, color='red'):
	from scipy.optimize import curve_fit
	p0 = [0.0, 0.0]
	try:
		coeffs, pcov = curve_fit(linear, x, y, sigma=err, p0=p0)
	except RuntimeError:
		print '********* Curve fit failed *********'
		return np.nan, np.nan
	xmin, xmax = ax.get_xlim()
	x_fit = np.linspace(xmin, xmax, 20)
	y_fit = linear(x_fit, coeffs[0], coeffs[1])
	ax.plot(x_fit, y_fit, color=color, linestyle='--')
	return ax, coeffs[0], pcov[0,0], coeffs[1], pcov[1,1]


def linear(x, slope, intercept):
	return slope * x + intercept


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


def save_fits(file, name, slope, slope_err, intercept, intercept_err):
	with open(file, 'a') as fd:
		fd.write("%s %g %g %g %g\n" % (name, slope, slope_err, intercept, intercept_err))


plot_dest_type = 'paper'
if plot_dest_type == 'paper':
	mpl.rcParams['font.family'] = 'serif'
	mpl.rcParams['font.size'] = 16
	mpl.rcParams['axes.linewidth'] = 3
	mpl.rcParams['lines.linewidth'] = 4
	mpl.rcParams['patch.linewidth'] = 4
	mpl.rcParams['xtick.major.width'] = 3
	mpl.rcParams['ytick.major.width'] = 3
	mpl.rcParams['xtick.major.size'] = 8
	mpl.rcParams['ytick.major.size'] = 8

#colors = ['red', 'green', 'blue']
colors = ['black', 'black']
labels1 = [r'$c$', r'$M_{\mathrm{vir}}$']
names = ['c_rockstar', 'Mvir']
xlabel = 'Redshift'
ylabels1 = [r'$\Delta c$ Slope $((\log(M_{\odot}))^{-1})$', r'$\Delta M_{\mathrm{vir}}$ Slope $((\log(M_{\odot}))^{-1})$']
ylabels2_kurt = [r'Kurtosis for $\Delta c$', r'Kurtosis for $\Delta M_{\mathrm{vir}}$']
ylabels2_skew = [r'Skew for $\Delta c$', r'Skew for $\Delta M_{\mathrm{vir}}$']
plot_base = 'plots/'
plot_ext  = '.eps'

statsfile = 'plots/stats.dat'

z_col        = -1
snap_col     = 0
slope_col     = 1
slope_err_col = 2
intercept_col      = 3
intercept_err_col  = 4

data_mean_col = 1
data_rms_col  = 15

#z_col        = -1
#snap_col     = 0
#mean_col     = 1
#mean_err_col = -2
#var_col      = 2
#var_err_col  = -2
#skew_col     = 3
#skew_err_col = -2
#kurt_col     = 4
#kurt_err_col = -2

offset = 0.06
#offset = 0.0

minsnap = 39
#minsnap = None

fit_trend     = True
separate_axes = True
plot_intercept = False


if __name__ == '__main__':
	main()

