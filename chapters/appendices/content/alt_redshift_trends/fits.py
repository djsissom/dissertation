#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import gamma as Gamma
from scipy.special import psi as digamma
from ipdb import set_trace


def main():
	if (len(sys.argv[1:]) == 4):
		data1 = read_files(sys.argv[1], header_line = None)
		data2 = read_files(sys.argv[2], header_line = None)
		data3 = read_files(sys.argv[3], header_line = None)
		rsnap_data = read_files(sys.argv[4], header_line = None)
	else:
		print 'need 4 files'
		sys.exit(15)
	
	if fit_mean_trend:
		with open(statsfile, 'w') as fd:
			fd.write("#plot slope slope_err intercept intercept_err\n")

	if minsnap > 0:
		data1 = data1[data1[:,0] >= minsnap]
		data2 = data2[data2[:,0] >= minsnap]
		data3 = data3[data3[:,0] >= minsnap]

	z = 1.0 / rsnap_data[:,1] - 1.0
	if (len(data1) == len(data2)) and (len(data1) == len(data3)):
		z = z[-len(data1):]
	else:
		sys.exit(16)
	
	data1 = np.column_stack((data1, z))
	data2 = np.column_stack((data2, z))
	data3 = np.column_stack((data3, z))


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# make mean and stdv plots                                                #
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	for (data, ylabel, label, name) in zip([data1, data2, data3], ylabels1, labels1, names):
		print "Making %s plot..." % (name + ' xvals')
		fig = plt.figure(figsize=(9.0, 6.0))
		ax = fig.add_subplot(111)

		ax = make_plot(ax, data[:,z_col], data[:,peak_col], err = None, color = 'black', marker='o', linestyle='-', label=None)

		for (x_val_col, color) in zip(x_val_cols, colors1):
			ax = make_plot(ax, data[:,z_col], data[:,x_val_col], err = None, color = color, marker='o', linestyle='--', label=None)

		#if add_rms_line:
		#	ax = make_plot(ax, data[:,z_col], data[:,data_rms_col], color = 'green', linestyle=':')

		#if fit_mean_trend:
		#	ax, slope, slope_err, intercept, intercept_err = add_fit(ax, data[:,z_col], data[:,mean_col], err=data[:,mean_err_col], color='red')
		#	save_fits(statsfile, name, slope, np.sqrt(slope_err), intercept, np.sqrt(intercept_err))

		#ax.legend(loc='lower right')
		ax.set_xlim(z[0] + 1.0, z[-1] - 1.0)
		#ax.invert_xaxis()

		ax.set_xlabel(xlabel, fontsize='x-large')
		ax.set_ylabel(ylabel, fontsize='x-large')

		fig.tight_layout()
		fig.savefig(plot_base + name + '_xvals' + plot_ext, bbox_inches='tight')

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	for (data, ylabel, label, name) in zip([data1, data2, data3], ylabels2, labels1, names):
		print "Making %s plot..." % (name + ' sumfrac')
		fig = plt.figure(figsize=(9.0, 6.0))
		ax = fig.add_subplot(111)

		for (sum_frac_col, color) in zip(sum_frac_cols, colors2):
			ax = make_plot(ax, data[:,z_col], data[:,sum_frac_col], err = None, color = color, marker='o', linestyle='-', label=None)
		for (doublesum_frac_col, color) in zip(doublesum_frac_cols, colors2):
			ax = make_plot(ax, data[:,z_col], data[:,doublesum_frac_col], err = None, color = color, marker='o', linestyle='--', label=None)

		ax.set_xlabel(xlabel, fontsize='x-large')
		ax.set_ylabel(ylabel, fontsize='x-large')
		ax.set_xlim(z[0] + 1.0, z[-1] - 1.0)
		ax.set_yscale('log')

		fig.tight_layout()
		fig.savefig(plot_base + name + '_sumfrac' + plot_ext, bbox_inches='tight')


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
colors = ['black', 'black', 'black']
labels1 = [r'$c$', r'$M_{\mathrm{vir}}$', r'$X_{\mathrm{off}}$']
names = ['c_rockstar', 'Mvir', 'Xoff']
xlabel = 'Redshift'
ylabels1 = [r"$\Delta' c(f_{h},z)$ and $\Delta' c_{\mathrm{peak}}$", r"$\Delta' M_{\mathrm{vir}}(f_{h},z)$ and $\Delta' M_{\mathrm{vir, peak}}$", r"$\Delta' X_{\mathrm{off}}(f_{h},z)$ and $\Delta' X_{\mathrm{off, peak}}$"]
ylabels2 = [r"$f_{h}(\Delta' c,z)$", r"$f_{h}(\Delta' M_{\mathrm{vir}},z)$", r"$f_{h}(\Delta' X_{\mathrm{off}},z)$"]
plot_base = 'plots/'
plot_ext  = '.eps'

statsfile = 'plots/stats.dat'

z_col        = -1
snap_col     = 0
mean_col     = 7
mean_err_col = 8
var_col      = 9
var_err_col  = 10
skew_col     = 3
skew_err_col = -2
#skew_col     = 7
#skew_err_col = 8
#kurt_col     = 4
#kurt_err_col = -2
kurt_col     = 13
kurt_err_col = 14
beta_col     = 13
beta_err_col = 14

data_mean_col = 1
data_rms_col  = 15



peak_col  = 1
x_val_cols = np.array([4, 6, 8]) + 2
sum_frac_cols = np.array([2, 4, 6, 8]) + 2 + 9
doublesum_frac_cols = sum_frac_cols + 9

colors1 = ['red', 'green', 'blue']
colors2 = ['blue', 'green', 'red', 'black']

offset = 0.06
#offset = 0.0

minsnap = 39
#minsnap = None

fit_mean_trend     = False
add_rms_line       = False


if __name__ == '__main__':
	main()

