#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
#from ipdb import set_trace



def main():
	# Read in particle files
	header, halos = read_files(sys.argv[1:], header_line = 3)

	if c_source == 'density_profile':
		print 'len(halos) = ', len(halos)
		halos = halos[np.isfinite(halos[:,c_2lpt_col])]
		halos = halos[np.isfinite(halos[:,c_za_col])]
		print 'len(halos) = ', len(halos)
	
	print 'Filtering data...'
	for col, val in zip(lt_cols, lt_vals):
		halos = halos[halos[:, col] <= val]
	for col, val in zip(gt_cols, gt_vals):
		halos = halos[halos[:, col] >= val]
	for col, val in zip(eq_cols, eq_vals):
		halos = halos[halos[:, col] == val]
	for col, val in zip(ne_cols, ne_vals):
		halos = halos[halos[:, col] != val]
		
	m_avg = (halos[:,47] + halos[:,48])/2.0
	halos = np.column_stack((halos, m_avg))
	header = np.append(header, 'M_avg')

	if x_min_lim > 0:
		print 'nhalos =', len(halos)
		mask = (m_avg >= x_min_lim)
		halos = halos[mask]
		print 'nhalos =', len(halos)

	if c_source == 'rockstar':
		c1 = halos[:, Rv1_col] / halos[:, Rs1_col]
		c2 = halos[:, Rv2_col] / halos[:, Rs2_col]
		if use_klypin:
			mask = (halos[:,4] < 100)
			c1[mask] = halos[mask, Rv1_col] / halos[mask, 79]
			mask = (halos[:,5] < 100)
			c1[mask] = halos[mask, Rv2_col] / halos[mask, 80]
	if c_source == 'density_profile':
		c1 = halos[:, c_2lpt_col]
		c2 = halos[:, c_za_col]

	dc = 2.0 * (c1 - c2) / (c1 + c2)
	#dc = c1 - c2

	m1 = halos[:,47]
	m2 = halos[:,48]
	dm = 2.0 * (m1 - m2) / (m1 + m2)

	for x_col, xlabel in zip(x_cols, xlabels):
		make_plot(halos[:, x_col], dm, x_col, header[x_col], xlabel, ylabel_m, plot_base_m, stats_file_m, y_lim_m, use_log=False)
		make_plot(halos[:, x_col], dc, x_col, header[x_col], xlabel, ylabel_c, plot_base_c, stats_file_c, y_lim_c, use_log=False)
	for x_col, xlabel in zip(x_log_cols, xlabels_log):
		make_plot(halos[:, x_col], dm, x_col, header[x_col], xlabel, ylabel_m, plot_base_m, stats_file_m, y_lim_m, use_log=True)
		make_plot(halos[:, x_col], dc, x_col, header[x_col], xlabel, ylabel_c, plot_base_c, stats_file_c, y_lim_c, use_log=True)
	
	print 'Finished all plots.'


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


def make_plot(x, y, x_col, header, xlabel, ylabel, plot_base, stats_file, y_lim, use_log):
	print 'generating plot...'
	fig = plt.figure(figsize=(9.0,6.0))
	ax = fig.add_subplot(1,1,1)
	ax = draw_hist2d(ax, x, y, y_lim)
	if fit_to_data:
		ax = draw_data_fit(ax, x, y, x.min(), x.max(), use_log=use_log)
	if fit_to_binned_data:
		mid_bins, mean, stdev, n = get_bin_avgs(x, y, use_log=use_log)
		ax = draw_bin_fit(ax, mid_bins, mean, stdev/np.sqrt(n), x.min(), x.max(), stats_file, use_log=use_log)
		ax = draw_bin_avgs(ax, mid_bins, mean, stdev, n, use_log=use_log)

	ax.set_xlim([x.min(), x.max()])
	#ax.set_yscale("log")
	ax.set_xlabel(xlabel, fontsize="x-large")
	ax.set_ylabel(ylabel, fontsize="x-large")
	
	fig.tight_layout()
	header = header.replace("/", "over")
	plot_name = "%s%s%0.3d%s%s%s" % (plot_base, '(', x_col, ')_', header, plot_ext)
	plt.savefig(plot_name, bbox_inches='tight')
	print 'finished plot ' + plot_name


def draw_hist2d(ax, x, y, y_lim):
	if use_log:
		xbins = np.logspace(np.log10(x.min()), np.log10(x.max()), num=nbins+1)
	else:
		xbins = np.linspace(x.min(), x.max(), num=nbins+1)

	ybins = np.linspace(y.min(), y.max(), num=nbins+1)

	if use_log:
		ax.set_xscale("log")
		im = my_hist2d(ax, x, y, bins=[xbins, ybins], zorder=-50)
	else:
		im = ax.hist2d(x, y, bins=[xbins, ybins], cmap=colormap, zorder=-50)

	if y_lim > 0.0:
		ax.set_ylim([-y_lim, y_lim])

	line = ax.plot([x.min(), x.max()], [0.0, 0.0], color='0.65', linestyle='--', linewidth=1, zorder=-20)
	return ax


def my_hist2d(ax, x, y, bins=10, range=None, normed=False, weights=None,
						cmin=None, cmax=None, **kwargs):
	import matplotlib as mpl

	bin_range = range
	range = mpl.axes.__builtins__["range"]
	h, xedges, yedges = np.histogram2d(x, y, bins=bins, range=bin_range,
	                                   normed=normed, weights=weights)

	if cmin is not None:
		h[h < cmin] = None
	if cmax is not None:
		h[h > cmax] = None

	if z_log:
		h[h<1.0] = 0.5
		h = np.log10(h)
	
	h = gaussian_filter(h, len(h) / 75.0)
	
	pc = ax.imshow(h[:,::-1].T, cmap=colormap, extent=[x.min(), x.max(), y.min(), y.max()], interpolation='gaussian', **kwargs)
	ax.set_xlim(xedges[0], xedges[-1])
	ax.set_ylim(yedges[0], yedges[-1])
	return h, xedges, yedges, pc


def get_bin_avgs(x, y, use_log):
	if use_log:
		fit_bins = np.logspace(np.log10(x.min()), np.log10(x.max()), num=nfit_bins+1)
	else:
		fit_bins = np.linspace(x.min(), x.max(), num=nfit_bins+1)
	
	mid_bins = (fit_bins[:-1] + fit_bins[1:]) / 2.0

	mean = np.array([])
	stdev = np.array([])
	n = np.array([])
	for xmin, xmax in zip(fit_bins[:-1], fit_bins[1:]):
		mask = np.logical_and(x > xmin, x <= xmax)
		if mask.sum() > 0:
			mean_el = y[mask].mean()
			#stdev_el = y[mask].std() / np.sqrt(len(y))
			stdev_el = y[mask].std()
			#stdev_el = stdev / np.sqrt(len(y[mask]))
			n_el = len(y[mask])
		else:
			mean_el = 0.0
			stdev_el = -1.0
			n_el = 0
		mean = np.append(mean, mean_el)
		stdev = np.append(stdev, stdev_el)
		n = np.append(n, n_el)
	
	mask = (n > 0)
	mean = mean[mask]
	stdev = stdev[mask]
	n = n[mask]
	mid_bins = mid_bins[mask]

	return mid_bins, mean, stdev, n


def draw_bin_avgs(ax, mid_bins, mean, stdev, n, use_log):
	ax.errorbar(mid_bins, mean, yerr=stdev/np.sqrt(n), fmt='o', color='black', linewidth=2)
	
	if draw_stdev_lines:
		ax.plot(mid_bins, mean + stdev, color='black', linestyle=':', linewidth=3, zorder=-15)
		ax.plot(mid_bins, mean - stdev, color='black', linestyle=':', linewidth=3, zorder=-15)
	return ax


def draw_bin_fit(ax, mid_bins, mean, stdev, x_min, x_max, stats_file, use_log):
	stdev[stdev == 0.0] = 0.1
	#fit data
	if use_log:
		#coefs, res, rank, singvals, rcond = np.polyfit(np.log10(mid_bins), mean, 1, full=True)
		coefs, pcov = curve_fit(linear, np.log10(mid_bins), mean, sigma=stdev, p0=[0.0, 0.0])
	else:
		#coefs, stats = np.polynomial.polynomial.polyfit(mid_bins, mean, 1, full=True)
		coefs, pcov = curve_fit(linear, mid_bins, mean, sigma=stdev, p0=[0.0, 0.0])
	print 'coefs = ', coefs

	
	m = coefs[0]
	b = coefs[1]
	m_err = pcov[0,0]
	b_err = pcov[1,1]

	if use_log:
		x = np.logspace(np.log10(x_min), np.log10(x_max), 100)
		y = m * np.log10(x) + b
	else:
		x = np.linspace(x_min, x_max, 100)
		y = m * x + b
	#y = x**m + b
	#line = ax.plot(x, y, color='white', linewidth=8)  # to avoid blending with colormap background
	line = ax.plot(x, y, color='magenta', zorder=-10)

	if print_fit_params:
		if use_log:
			textstr = '$y = m \log x + b$\n$m = %g$\n$b = %g$' % (m, b)
		else:
			textstr = '$y = m x + b$\n$m = %g$\n$b = %g$' % (m, b)
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=14,
				verticalalignment='top', bbox=props)

	if save_fit_params:
		with open(stats_file, "a") as fd:
			fd.write("%g  %g  %g  %g\n" % (m, m_err, b, b_err))

	return ax


def linear(x, slope, intercept):
	return slope * x + intercept


def draw_data_fit(ax, x, y, x_min, x_max, use_log):
	if remove_zero_strip:
		mask = (np.abs(y) >= y_epsilon)
		x = x[mask]
		y = y[mask]

	#fit data
	if use_log:
		coefs, residual, rank, singular_values, rcond = np.polyfit(np.log10(x), y, 1, full=True)
	#	coefs, stats = np.polynomial.polynomial.polyfit(np.log10(mid_bins), mean, 1, w=1.0/stdev, full=True)
	#	coefs, res, rank, singvals, rcond = np.polyfit(np.log10(mid_bins), mean, 1, full=True)
	else:
		coefs, residual, rank, singular_values, rcond = np.polyfit(x, y, 1, full=True)
	#	coefs, stats = np.polynomial.polynomial.polyfit(mid_bins, mean, 1, w=1.0/stdev, full=True)
	#	coefs, stats = np.polynomial.polynomial.polyfit(mid_bins, mean, 1, full=True)
	print 'coefs =', coefs, '+/-', residual 

	
	m = coefs[0]
	b = coefs[1]
	if use_log:
		x = np.logspace(np.log10(x_min), np.log10(x_max), 100)
		y = m * np.log10(x) + b
	else:
		x = np.linspace(x_min, x_max, 100)
		y = m * x + b
	#y = x**m + b
	line = ax.plot(x, y, color='red')

	if print_fit_params:
		if use_log:
			textstr = '$y = m \log x + b$\n$m = %g$\n$b = %g$' % (m, b)
		else:
			textstr = '$y = m x + b$\n$m = %g$\n$b = %g$' % (m, b)
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=14,
				verticalalignment='top', bbox=props)

	if save_fit_params:
		with open("fits_to_data.dat", "a") as fd:
			fd.write("%g %g %g\n" % (m, b, residual))

	return ax


use_log = True
#use_log = False
z_log = True

#fit_bins = True
#fit_data = True

print_fit_params = False
save_fit_params = True

use_klypin = True

remove_zero_strip = False
y_epsilon = 0.01

y_lim_m = 0.5
y_lim_c = 1.0
x_min_lim = 5.33e5 * 100

#if use_log:
#  x_cols = [4, 5, 6, 9, 10, 23, 24, 31, 32, 47, 48, 51, 52, 57, 58] # log10 columns
#else:
#  x_cols = [17, 18, 77, 78, 91, 92, 93, 94, 97, 98, 99, 100, 107, 108, 111, 112] # nolog columns

x_cols     = []
x_log_cols = [-1]
#x_log_cols = [47, 48, -1]

xlabels = []
xlabels_log = [r"$M_{\mathrm{vir,avg}} \, \mathrm{(M_{\odot})}$"]
#xlabels_log = [r"$\mathrm{M_{2LPT} (M_{\odot})}$",
#               r"$\mathrm{M_{ZA} (M_{\odot})}$",
#			   r"$\mathrm{M_{avg} (M_{\odot})}$"]

#ylabel = r"$\mathrm{(M_{2LPT} - M_{ZA}) / M_{avg}}$"
ylabel_m = r"$(M_{\mathrm{vir,2LPT}} \, - \, M_{\mathrm{vir,ZA}}) \, / \, M_{\mathrm{vir,avg}}$"
ylabel_c = r"$(c_{\mathrm{2LPT}} \, - \, c_{\mathrm{ZA}}) \, / \, c_{\mathrm{avg}}$"

#c_source = 'density_profile'
c_source = 'rockstar'

plot_base_m = 'plots/diff_M_-_vs_-_'
plot_base_c = 'plots/diff_c_-_vs_-_'
plot_ext  = '.eps'

stats_file_m = 'fits_to_bins_m.dat'
stats_file_c = 'fits_to_bins_c.dat'

#plot_name = 'test.eps'
#plot_name = 'c_v_M200c_2lpt.eps'
fit_to_binned_data = True
fit_to_data = False
draw_stdev_lines = True

Rv1_col = 53
Rv2_col = 54
Rs1_col = 55
Rs2_col = 56

c_2lpt_col = 17
c_za_col   = 18

nbins = 100
nfit_bins = 10

## c_2lpt, c_za, chi2_2lpt, chi2_za
#lt_cols = [17, 18, 37, 38]
#lt_vals = [100.0, 100.0, 10.0, 10.0]
#
## c_2lpt, c_za, chi2_2lpt, chi2_za
#gt_cols = [17, 18, 37, 38]
#gt_vals = [1.0, 1.0, 0.0, 0.0]

lt_cols = []
lt_vals = []

gt_cols = [4, 5]
gt_vals = [100, 100]

eq_cols = [109, 110]
eq_vals = [-1, -1]

ne_cols = []
ne_vals = []

#colormap = cm.PuBuGn
#colormap = cm.cubehelix_r
#colormap = cm.ocean_r
#colormap = cm.rainbow
#colormap = cm.gnuplot2_r
#colormap = cm.CMRmap_r

def add_white(orig_map, num):
	temp_cmap = cm.get_cmap(orig_map, num)
	vals = temp_cmap(np.arange(num))
	nfade = num / 7
	vals[:nfade,0] = np.linspace(1., vals[nfade-1,0], nfade)
	vals[:nfade,1] = np.linspace(1., vals[nfade-1,1], nfade)
	vals[:nfade,2] = np.linspace(1., vals[nfade-1,2], nfade)
	#vals[:nfade,3] = np.linspace(0., vals[nfade-1,3], nfade)
	#vals[0] = [1.0, 1.0, 1.0, 1.0]
	#vals[1] = (vals[1] + [1.0, 1.0, 1.0, 1.0]) / 2.0
	newcmap = mpl.colors.LinearSegmentedColormap.from_list("custom_1", vals)
	return newcmap

colormap = add_white('rainbow', 30)

plot_dest_type = 'paper'
if plot_dest_type == 'paper':
	mpl.rcParams['font.family'] = 'serif'
	mpl.rcParams['font.size'] = 16
	mpl.rcParams['axes.linewidth'] = 3
	mpl.rcParams['lines.linewidth'] = 4
	#mpl.rcParams['lines.linewidth'] = 3
	mpl.rcParams['patch.linewidth'] = 4
	#mpl.rcParams['patch.linewidth'] = 3
	mpl.rcParams['xtick.major.width'] = 3
	mpl.rcParams['ytick.major.width'] = 3
	mpl.rcParams['xtick.major.size'] = 8
	mpl.rcParams['ytick.major.size'] = 8
	mpl.rcParams['xtick.minor.width'] = 2
	mpl.rcParams['ytick.minor.width'] = 2
	mpl.rcParams['xtick.minor.size'] = 4
	mpl.rcParams['ytick.minor.size'] = 4
	#mpl.rcParams['lines.antialiased'] = True


if __name__ == '__main__':
	main()

