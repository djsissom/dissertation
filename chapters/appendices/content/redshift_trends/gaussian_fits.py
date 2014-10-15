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

	if skew_err_boxes:
		skew_err1 = get_skew_err(sys.argv[1])
		skew_err2 = get_skew_err(sys.argv[2])
		skew_err3 = get_skew_err(sys.argv[3])

	if minsnap > 0:
		#for data in data1, data2, data3:
		#	data = data[data[:,0] >= minsnap]
		data1 = data1[data1[:,0] >= minsnap]
		data2 = data2[data2[:,0] >= minsnap]
		data3 = data3[data3[:,0] >= minsnap]
		if skew_err_boxes:
			skew_err1 = skew_err1[-len(data1):]
			skew_err2 = skew_err2[-len(data2):]
			skew_err3 = skew_err3[-len(data3):]

	if skew_err_col == -2:
		data1 = np.column_stack((data1, skew_err1))
		data2 = np.column_stack((data2, skew_err2))
		data3 = np.column_stack((data3, skew_err3))

	#if (mean_err_col == -2) or (var_err_col == -2) or (skew_err_col == -2) or (kurt_err_col == -2):
	#	fake_err = np.zeros(len(data1))
	#	data1 = np.column_stack((data1, fake_err))
	#	data2 = np.column_stack((data2, fake_err))
	#	data3 = np.column_stack((data3, fake_err))

	z = 1.0 / rsnap_data[:,1] - 1.0
	if (len(data1) == len(data2)) and (len(data1) == len(data3)):
		z = z[-len(data1):]
	else:
		sys.exit(16)
	
	data1 = np.column_stack((data1, z))
	data2 = np.column_stack((data2, z))
	data3 = np.column_stack((data3, z))

	#data1[:,-1] = data1[:,-1] - 0.12
	#data2[:,-1] = data2[:,-1] + 0.12

	for data in [data1, data2, data3]:
		if expand_error:
			mask = (np.abs(data[:,data_mean_col] - data[:,mean_col]) > data[:,mean_err_col])
			data[mask,mean_err_col] = np.abs(data[mask,data_mean_col] - data[mask,mean_col])
		if transform_variance:
			data[:,var_col] = data[:,var_col]**2 * Gamma(3.0 / data[:,beta_col]) / Gamma(1.0 / data[:,beta_col])
			data[:,var_err_col] = data[:,var_err_col]**2 * Gamma(3.0 / data[:,beta_col]) / Gamma(1.0 / data[:,beta_col])
		if transform_kurtosis:
			#data[:,kurt_col] = ( Gamma(5.0 / data[:,kurt_col]) * Gamma(1.0 / data[:,kurt_col]) / Gamma(3.0 / data[:,kurt_col]) ) - 3.0
			beta = data[:,beta_col]
			beta_err = data[:,beta_err_col]
			kurtosis = ( Gamma(5.0 / beta) * Gamma(1.0 / beta) / Gamma(3.0 / beta) ) - 3.0
			kurtosis_err = beta_err * (1.0 / beta**2) * (kurtosis + 3) * (6.0 * digamma(3.0/beta) - 5.0 * digamma(5.0/beta) - digamma(1.0/beta))

			data[:,kurt_col] = kurtosis[:]
			data[:,kurt_err_col] = kurtosis_err[:]

		data[:,var_col] = np.sqrt(data[:,var_col])       # var to stdev
		data[:,var_err_col] = np.sqrt(data[:,var_err_col])       # var to stdev


	if save_transformed_data:
		for data, path in zip([data1, data2, data3], sys.argv[1:4]):
			fname = transform_file_base + os.path.basename(path)
			with open(fname, 'w') as fd:
				fd.write(transformed_data_header)
				np.savetxt(fd, np.column_stack((z, data)), fmt='%g')

	

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# make mean and stdv plots                                                #
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	for (data, ylabel, color, label, name) in zip([data1, data2, data3], ylabels1, colors, labels1, names):
		print "Making %s plot..." % (name)
		fig = plt.figure(figsize=(9.0, 6.0))
		ax = fig.add_subplot(111)

		ax = make_plot(ax, data[:,z_col], data[:,mean_col], err = data[:,mean_err_col], color = 'blue', marker='o', label=label)
		ax = make_plot(ax, data[:,z_col], data[:,mean_col] + data[:,var_col], color = 'black', linestyle='--')
		ax = make_plot(ax, data[:,z_col], data[:,mean_col] - data[:,var_col], color = 'black', linestyle='--')

		if add_rms_line:
			ax = make_plot(ax, data[:,z_col], data[:,data_rms_col], color = 'green', linestyle=':')

		if fit_mean_trend:
			ax, slope, slope_err, intercept, intercept_err = add_fit(ax, data[:,z_col], data[:,mean_col], err=data[:,mean_err_col], color='red')
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

	#data1[:,-1] = data1[:,-1] + 0.12
	#data2[:,-1] = data2[:,-1] - 0.12


	for (data, ylabel_kurt, ylabel_skew, color, name, ylim_low1, ylim_high1, ylim_low2, ylim_high2) in zip([data1, data2, data3], ylabels2_kurt, ylabels2_skew, colors, names, [-10.0, -10.0, -1.0], [20.0, 20.0, 1.5], [-0.2, -1.5, -0.4], [0.5, 3.5, 0.1]):
		print "Making %s plot..." % (name)
		fig = plt.figure(figsize=(9.0, 6.0))
		ax = fig.add_subplot(111)

		#ax = make_plot(ax, data[:,z_col] - offset, data[:,kurt_col], err = data[:,kurt_err_col], color = 'red', marker='o', linestyle='-', label='Kurtosis')
		#ax = make_plot(ax, data[:,z_col] + offset, data[:,skew_col], err = data[:,skew_err_col], color = 'blue', marker='o', linestyle='-', label='Skew')
		ax = make_plot(ax, data[:,z_col] - offset, data[:,kurt_col], err = data[:,kurt_err_col], color = 'red', marker='o', linestyle=':', label='Kurtosis')
		legend_lines1, legend_labels1 = ax.get_legend_handles_labels()

		ax.set_xlabel(xlabel, fontsize='x-large')
		ax.set_ylabel(ylabel_kurt, fontsize='x-large')
		ax.set_ylim(ylim_low1, ylim_high1)

		if separate_skew_axes:
			ax = ax.twinx()
		ax = make_plot(ax, data[:,z_col] + offset, data[:,skew_col], err = data[:,skew_err_col], color = 'blue', marker='o', linestyle=':', label='Skew')
		legend_lines2, legend_labels2 = ax.get_legend_handles_labels()

		ax.set_ylabel(ylabel_skew, fontsize='x-large')
		ax.legend(legend_lines1 + legend_lines2, legend_labels1 + legend_labels2, loc='lower right')
		ax.set_xlim(z[0] + 1.0, z[-1] - 1.0)
		ax.set_ylim(ylim_low2, ylim_high2)
		#ax.invert_xaxis()

		fig.tight_layout()
		fig.savefig(plot_base + 'skew_kurtosis_' + name + plot_ext, bbox_inches='tight')

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


def get_skew_err(filebase):
	z = None
	skew = None
	for i in range(3):
		filename = filebase.replace('plots_all_snaps', 'plots_all_snaps_box'+str(i+1))
		data = read_files(filename, header_line = None)

		if i == 0:
			min_length = len(data)
		elif len(data) < min_length:
			min_length = len(data)

		if z == None:
			z = data[-min_length:,snap_col]
		else:
			z = np.column_stack((z[-min_length:], data[-min_length:,snap_col]))

		if skew == None:
			skew = data[-min_length:,skew_col]
		else:
			skew = np.column_stack((skew[-min_length:], data[-min_length:,skew_col]))
	
	if (z[:,0] != z[:,1]).all() or (z[:,0] != z[:,2]).all():
		print 'Need matching snapshots for skew error from individual boxes.'
		print z
		sys.exit(-1)
	
	skew_err = np.std(skew, axis=1) / np.sqrt(3.0)
	return skew_err


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
ylabels1 = [r'$\mu$ and $\sigma$ for $\Delta c$', r'$\mu$ and $\sigma$ for $\Delta M_{\mathrm{vir}}$', r'$\mu$ and $\sigma$ for $\Delta X_{\mathrm{off}}$']
ylabels2_kurt = [r'Kurtosis for $\Delta c$', r'Kurtosis for $\Delta M_{\mathrm{vir}}$', r'Kurtosis for $\Delta X_{\mathrm{off}}$']
ylabels2_skew = [r'Skew for $\Delta c$', r'Skew for $\Delta M_{\mathrm{vir}}$', r'Skew for $\Delta X_{\mathrm{off}}$']
plot_base = 'plots/'
plot_ext  = '.eps'

statsfile = 'plots/stats.dat'
transform_file_base = 'plots/'
transformed_data_header = '#z  snap  data_mean  data_stdev  data_skew  data_kurt  fit_height  +/-err  fit_mean  +/-err  fit_stdev  +/-err  fit_skew  +/-err  fit_kurt  +/-err  data_rms  data_gt_epsilon  chi2  pval  skew_err  z\n'

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

transform_variance = True
transform_kurtosis = True
expand_error       = True
fit_mean_trend     = True
separate_skew_axes = True
skew_err_boxes     = True
add_rms_line       = True
save_transformed_data = True


if __name__ == '__main__':
	main()

