#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.special import gamma as gamma_func
from scipy.optimize import curve_fit
import statsmodels.sandbox.distributions.extras as extrastats
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
	halos = np.column_stack((halos, c_rockstar_2lpt, c_rockstar_za, c_diff_2lpt, c_diff_za))
	header.append('c_rockstar')
	header.append('c_rockstar')
	header.append('c_diff')
	header.append('c_diff')

	if mass_quartiles and len(halos) > 50:
		start_fracs = [0.0, 0.25, 0.50, 0.75, 0.0]
		end_fracs   = [0.25, 0.50, 0.75, 1.0, 1.0]
	else:
		start_fracs = [0.0]
		end_fracs   = [1.0]
		

	for start_frac, end_frac in zip(start_fracs, end_fracs):
		halos_to_pass = halos[start_frac * len(halos) : end_frac * len(halos)]
		if use_alt_frac and (start_frac == 0.0) and (end_frac == 1.0):
			alt_halos_to_pass = halos[alt_start_frac * len(halos) : alt_end_frac * len(halos)]
		else:
			alt_halos_to_pass = None
		if len(halos_to_pass) > 0:
			for (lpt_col, za_col, fancy_x_label) in zip(lpt_log_cols, za_log_cols, fancy_log_x_labels):
				make_plot(halos_to_pass, alt_halos_to_pass, lpt_col, za_col, start_frac, end_frac, fancy_x_label, header, use_log=True)
			for (lpt_col, za_col, fancy_x_label) in zip(lpt_cols, za_cols, fancy_x_labels):
				make_plot(halos_to_pass, alt_halos_to_pass, lpt_col, za_col, start_frac, end_frac, fancy_x_label, header, use_log=False)

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


def sort_by_column(halos, col):
	print 'Sorting halos...'
	mask = np.argsort(halos[:, col])
	mask = mask[::-1]
	halos = halos[mask]
	return halos


def make_plot(halos, alt_halos, lpt_col, za_col, start_frac, end_frac, fancy_x_label, header=None, use_log=False):
	print 'start =', start_frac
	print 'end =', end_frac
	x_lpt = halos[:, lpt_col]
	x_za  = halos[:, za_col]
	x_lpt, x_za = filter(x_lpt, x_za, lpt_col, za_col)

	if alt_halos != None:
		alt_x_lpt = alt_halos[:, lpt_col]
		alt_x_za  = alt_halos[:, za_col]
		alt_x_lpt, alt_x_za = filter(alt_x_lpt, alt_x_za, lpt_col, za_col)

	if header != None:
		header_lpt = header[lpt_col]
		header_za   = header[za_col]
		if header_lpt == header_za:
			xlabel = header_lpt
			xlabel = xlabel.replace('/', '_over_')
		else:
			print 'column mismatch... exiting'
			set_trace()
			sys.exit(123)

	if len(x_lpt) == 0 or len(x_za) == 0:
		print "Skipping range %f - %f for %s plot.  No halos found." % (start_frac, end_frac, xlabel)
		return
		#set_trace()
	
	if perc_diff:
		print 'Finding percent difference stats...'
		x_perc_diff = (x_lpt - x_za) / x_za
		perc_diff_file = "%s%s%0.3d%s%0.3d%s%s_(%s-%s)%s" % \
	                     (perc_diff_base, '(', lpt_col, ',', za_col, ')_', xlabel, start_frac, end_frac, stats_ext)
		perc_diff_stats(x_perc_diff, perc_diff_file, use_log=use_log)
		print 'done.'

	x = 2.0 * (x_lpt - x_za) / (x_lpt + x_za)
	x[np.logical_and(x_lpt == 0, x_za == 0)] = 0

	if alt_halos != None:
		alt_x = 2.0 * (alt_x_lpt - alt_x_za) / (alt_x_lpt + alt_x_za)
		alt_x[np.logical_and(alt_x_lpt == 0, alt_x_za == 0)] = 0

#	set_trace()
	
	if x_lim == None:
		#x_max = max(abs(x.max()), abs(x.min()))
		if lpt_col == 47:
			x_max = x.mean() + x.std() * 1.5
			x_min = x.mean() - x.std() * 1.5
		else:
			x_max = np.std(x) * 3.0
			x_min = -x_max
	else:
		x_max = x_lim
		x_min = -x_lim

	# get stats
	data_mean = x.mean()
	data_stdev = x.std()**2
	data_skew = stats.skew(x)
	data_kurt = stats.kurtosis(x)
	data_rms = np.sqrt(np.mean(x**2))
	data_gt_epsilon = float(len(x[np.abs(x) >= 0.1])) / float(len(x))

	# Generate plot
	print 'generating', xlabel, 'plot...'
	fig = plt.figure(figsize=(9.0, 6.0))
	if add_residuals_panel:
		grid = gridspec.GridSpec(2, 1, height_ratios=[1,4])
		ax = fig.add_subplot(grid[1])
	else:
		ax = fig.add_subplot(111)
	ax, n, bins, patches = draw_hist(ax, x, x_min=x_min, x_max=x_max, \
	                                 use_log=use_log, color='blue', fill=None)

	p0 = [1.0, data_mean, data_stdev, 2.0]
	ax, fit_height, fit_mean, fit_stdev, fit_skew, fit_kurt, fit_height_err, fit_mean_err, fit_stdev_err, fit_skew_err, fit_kurt_err, chi2, pval = draw_fit(ax, n, bins, p0)

	if draw_data_fit:
		ax = draw_data_gaussian(ax, x, n, bins)

	if alt_halos != None:
		ax, n_alt, bins_alt, patches_alt = draw_hist(ax, alt_x, x_min=x_min, x_max=x_max, \
												 use_log=use_log, color='green', fill="0.75")
		#ax = draw_fit(ax, n, bins)

	#ax.grid(color='gray', linestyle='dashed')
	ax.set_xlim([x_min, x_max])
	#ax.set_xlabel('(' + xlabel + '_2lpt - ' + xlabel + '_za) / ' + xlabel + '_avg')
	#ax.set_ylabel(ylabel)
	if label_axes:
		ax.set_xlabel(fancy_x_label, fontsize="xx-large")
		ax.set_ylabel(fancy_y_label, fontsize="xx-large")
	#ax.legend()
	
	if add_residuals_panel:
		ax = fig.add_subplot(grid[0])
		ax = draw_residuals(ax, n, bins, fit_height, fit_mean, fit_stdev, fit_kurt)
		ax.tick_params(axis='x', labelbottom='off')

	fig.tight_layout()
	plot_name = "%s%s%0.3d%s%0.3d%s%s_(%s-%s)%s" % \
	            (plot_base, '(', lpt_col, ',', za_col, ')_', xlabel, start_frac, end_frac, plot_ext)
	fig.savefig(plot_name, bbox_inches='tight')

	if save_stats:
		statsfile = "%s%s%0.3d%s%0.3d%s%s_(%s-%s)%s" % \
	                (stats_base, '(', lpt_col, ',', za_col, ')_', xlabel, start_frac, end_frac, stats_ext)
		with open(statsfile, 'w') as fd:
			if bin_test:
				for ntestbins in range(nbins_min, nbins_max+1, 5):
					fit_mean, fit_stdev = rebin_stats(ntestbins, x, x_min=x_min, x_max=x_max, use_log=use_log)
					fd.write("%d %g %g %g %g\n" % (ntestbins, data_mean, data_stdev, fit_mean, fit_stdev))
			else:
				fd.write("%d   %g %g %g %g   %g %g %g %g %g %g %g %g %g %g   %g %g   %g %g\n" % \
				         (nbins, data_mean, data_stdev, data_skew, data_kurt, \
						  fit_height, fit_height_err, fit_mean, fit_mean_err, fit_stdev, fit_stdev_err, fit_skew, fit_skew_err, fit_kurt, fit_kurt_err, \
						  data_rms, data_gt_epsilon, chi2, pval))

	print 'finished plot ' + plot_name
	return


def perc_diff_stats(x, filename, use_log=False):
	data_mean = x.mean()
	data_stdev = x.std()**2
	data_skew = stats.skew(x)
	data_kurt = stats.kurtosis(x)
	data_rms = np.sqrt(np.mean(x**2))
	data_gt_epsilon = float(len(x[np.abs(x) >= 0.1])) / float(len(x))

	if x_lim == None:
		x_max = min((x.mean() + x.std() * 3.0), x.max())
		x_min = max((x.mean() - x.std() * 3.0), x.min())
	else:
		x_max = x_lim
		x_min = -x_max

	global nbins
	if nbins <= 0:
		nbins = np.sqrt(len(x))
		if nbins % 2 == 0:
			nbins = nbins - 1
	if nbins < nbins_min:
		nbins = nbins_min
	elif nbins > nbins_max:
		nbins = nbins_max

	if use_log:
		xbins = np.logspace(np.log10(x_min), np.log10(x_max), num=nbins+1)
		mid_bins = 10.0**(0.5 * (np.log10(xbins[1:]) + np.log10(xbins[:-1])))
	else:
		xbins = np.linspace(x_min, x_max, num=nbins+1)
		mid_bins = 0.5 * (xbins[1:] + xbins[:-1])

	hist, bin_edges = np.histogram(x, bins=xbins)
	x_peak = mid_bins[hist == hist.max()][0]

	x_sorted = np.sort(x)
	n_halos = len(x_sorted)

	x_vals = []
	for frac in fractions:
		x_vals.append(x_sorted[len(x_sorted)*frac])
	x_vals = np.array(x_vals)

	sum_frac_halos = []
	for diff_val in diff_vals:
		n_gt_val = (x_sorted >= diff_val).sum()
		sum_frac_halos.append(float(n_gt_val) / float(n_halos))
	sum_frac_halos = np.array(sum_frac_halos)

	doublesum_frac_halos = []
	for right_diff_val in diff_vals:
		left_diff_val = (1.0 / (right_diff_val + 1.0)) - 1.0
		n_gt_val = (x_sorted >= right_diff_val).sum() + (x_sorted <= left_diff_val).sum()
		doublesum_frac_halos.append(float(n_gt_val) / float(n_halos))
	doublesum_frac_halos = np.array(doublesum_frac_halos)

	with open(filename, 'w') as fd:
		fd.write("%d   %g   %s   %s   %s   %g %g %g %g   %g %g\n" % \
		         (nbins, x_peak, \
				  ' '.join("%g" % x for x in x_vals), \
				  ' '.join("%g" % x for x in sum_frac_halos), \
				  ' '.join("%g" % x for x in doublesum_frac_halos), \
		          data_mean, data_stdev, data_skew, data_kurt, \
		          data_rms, data_gt_epsilon))

	return


def find_frac_bounds(hist, start_bin, frac):
	n_tot = hist.sum()
	n_sum = hist[start_bin]

	left_tot = hist[:start_bin].sum() + hist[start_bin]/2.0
	right_tot = hist[start_bin+1:].sum() + hist[start_bin]/2.0

	if float(left_tot) / float(n_tot) <= frac / 2.0:
		right_only = True
	if float(right_tot) / float(n_tot) <= frac / 2.0:
		left_only = True

	left_bound = start_bin
	right_bound = start_bin
	while(float(n_sum) / float(n_tot) < frac):
		
		pass

	return left_bound, right_bound


def filter(x_lpt, x_za, lpt_col, za_col):
	mask = np.isfinite(x_lpt)
	x_lpt = x_lpt[mask]
	x_za  = x_za[mask]
	mask = np.isfinite(x_za)
	x_lpt = x_lpt[mask]
	x_za  = x_za[mask]

	if column_filter_halos:
		x_lpt, x_za = filter_columns(lpt_col, x_lpt, x_za)
		x_za, x_lpt = filter_columns(za_col, x_za, x_lpt)
	
	return x_lpt, x_za


def filter_columns(x_col, x1, x2):
	print 'Filtering data...'

	mask = np.isfinite(x1)
	x1 = x1[mask]
	x2 = x2[mask]

	mask = (x1 != -9999)
	x1 = x1[mask]
	x2 = x2[mask]

	if x_col in lt_cols:
		val = lt_vals[lt_cols.index(x_col)]
		mask = (x1 <= val)
		x1 = x1[mask]
		x2 = x2[mask]
	if x_col in gt_cols:
		val = gt_vals[gt_cols.index(x_col)]
		mask = (x1 >= val)
		x1 = x1[mask]
		x2 = x2[mask]
	if x_col in eq_cols:
		val = eq_vals[eq_cols.index(x_col)]
		mask = (x1 == val)
		x1 = x1[mask]
		x2 = x2[mask]
	if x_col in ne_cols:
		val = ne_vals[ne_cols.index(x_col)]
		mask = (x1 != val)
		x1 = x1[mask]
		x2 = x2[mask]
	return x1, x2


def draw_hist(ax, x, x_min=None, x_max=None, use_log=False, color=None, fill=None, label=None):
	global nbins
	if nbins <= 0:
		nbins = np.sqrt(len(x))
		if nbins % 2 == 0:
			nbins = nbins - 1
	if nbins < nbins_min:
		nbins = nbins_min
	elif nbins > nbins_max:
		nbins = nbins_max

	if use_log:
		xbins = np.logspace(np.log10(x_min), np.log10(x_max), num=nbins+1)
		ax.set_xscale('log')
	else:
		xbins = np.linspace(x_min, x_max, num=nbins+1)

	if fill == None:
		type='step'
	else:
		type='stepfilled'

	n, bins, patches = ax.hist(x, bins=xbins, histtype=type, facecolor=fill, normed=hist_normed, cumulative=hist_cumulative, log=ylog, edgecolor=color, label=label)
	return ax, n, bins, patches


def draw_fit(ax, hist, bin_edges, p0):
	bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

	if ignore_central_bin:
		mask = (np.abs(bin_centers) > 0.000001)
		bin_centers = bin_centers[mask]
		hist = hist[mask]
	
	hist[hist==0] = 1  #fix devide by zero error

	try:
		if poisson_weight:
			sigma=np.sqrt(hist)/hist
			sigma = sigma / float(hist.max())
		else:
			sigma=None

		if fit_in_log:
			#if sigma != None:
			#	sigma = np.log10(sigma)

			coeffs, var_matrix = curve_fit(log_generalized_normal, bin_centers, np.log10(hist/float(hist.max())), p0=p0, sigma=sigma)

			coeffs[0] = coeffs[0]**2
			var_matrix[0,0] = var_matrix[0,0]**2
		else:
			coeffs, var_matrix = curve_fit(generalized_normal, bin_centers, hist/float(hist.max()), p0=p0, sigma=sigma)

		if prevent_small_shape_param and coeffs[3] < 1.0:
			coeffs[3] = 1.0 / coeffs[3]
		print 'coeffs =', coeffs

	except RuntimeError:
		print '*******curve_fit failed!'
		return ax, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

	height, mean, stdv, skew, kurt = coeffs[0] * hist.max(), coeffs[1], coeffs[2], 0.0, coeffs[3]
	height_err, mean_err, stdv_err, skew_err, kurt_err = np.sqrt(var_matrix[0,0]*hist.max()), np.sqrt(var_matrix[1,1]), np.sqrt(var_matrix[2,2]), 0.0, np.sqrt(var_matrix[3,3])

	fit_x = np.linspace(bin_edges[0], bin_edges[-1], nfitpoints+1)
	hist_fit = generalized_normal(fit_x, height, mean, stdv, kurt)
	ax.plot(fit_x, hist_fit, color='red', linestyle='--')

	chi2_fit = generalized_normal(bin_centers, height, mean, stdv, kurt)
	chi2, pval = stats.chisquare(hist / hist.max(), chi2_fit / hist.max())

	return ax, height, mean, stdv, skew, kurt, height_err, mean_err, stdv_err, skew_err, kurt_err, chi2, pval


def draw_residuals(ax, hist, bin_edges, fit_height, fit_mean, fit_stdev, fit_kurt):
	bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
	fit = generalized_normal(bin_centers, fit_height, fit_mean, fit_stdev, fit_kurt)
	ratio = (hist - fit) / hist.max()
	#ax.plot(bin_centers, ratio, linestyle='steps-mid-')
	ax.plot(bin_centers, ratio, linestyle='steps-mid-')
	return ax


def draw_data_gaussian(ax, x, hist, bins):
	bin_centers = (bins[:-1] + bins[1:]) / 2.0
	x_min = bins[0]
	x_max = bins[-1]

	mean = np.mean(x)
	stdv = np.std(x)**2
	skew = stats.skew(x)
	kurt = stats.kurtosis(x)

	print "data stats:  mean = %g  stdv = %g  skew = %g  kurt = %g" % (mean, stdv, skew, kurt)

	coeffs, var_matrix = curve_fit(gaussian_height(mean, stdv, skew, kurt), bin_centers, hist, p0=[hist.max()])
	height = coeffs[0]

	fit_x = np.linspace(x_min, x_max, nfitpoints+1)
	hist_fit = gaussian(fit_x, height, mean, stdv, skew, kurt)
	ax.plot(fit_x, hist_fit, color='0.25', linestyle='-.')
	return ax


#def gaussian(x, A, mu, sigma, skew, kurtosis):
#	pdf_function = extrastats.pdf_mvsk([mu, sigma, skew, kurtosis])
#	return A * pdf_function(x)


def double_gaussian(x, A, mu, sigma, skew, kurtosis, A2, mu2, sigma2, skew2, kurtosis2):
	return gaussian(x, A, mu, sigma, skew, kurtosis) + gaussian(x, A2, mu2, sigma2, skew2, kurtosis2)


def gaussian_height(mu, sigma, skew, kurtosis):
	def func(x, A):
		pdf_function = extrastats.pdf_mvsk([mu, sigma, skew, kurtosis])
		return A * pdf_function(x)
	return func


#def log_gaussian(x, A, mu, sigma, skew=0.0, kurtosis=0.0):
def log_gaussian(x, A, mu, sigma):
	A = A**2  # remember to also square fit value for A
	y = gaussian(x, A, mu, sigma)
	#y = gaussian(x, A, mu, sigma, skew, kurtosis)
	if (y <= 0.0).any():
		#y[y<=0] = -y[y<=0] + 1
		y[y<=0] = (y[y<=0] + 0.0001)**2
	return np.log10(y)


#def log_double_gaussian(x, A1, mu1, sigma1, skew1, kurtosis1, A2, sigma2, skew2, kurtosis2):    # for common mean
#def log_double_gaussian(x, A1, mu1, sigma1, skew1, kurtosis1, A2, mu2, sigma2, skew2, kurtosis2):
def log_double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
	#mu2 = mu1  # for common mean
	A1 = A1**2  # remember to also square fit value for A
	A2 = A2**2
	skew1 = 0.0
	skew2 = 0.0
	kurtosis1 = 0.0
	kurtosis2 = 0.0
	y = double_gaussian(x, A1, mu1, sigma1, skew1, kurtosis1, A2, mu2, sigma2, skew2, kurtosis2)
	if (y <= 0.0).any():
		#y[y<=0] = -y[y<=0] + 1
		y[y<=0] = (y[y<=0] + 0.0001)**2
	return np.log10(y)


def gaussian(x, A, mu, sigma):
	return A * np.exp(-(x - mu)**2 / (2.0 * sigma**2))


def generalized_normal(x, A, mu, alpha, beta):
	if prevent_small_shape_param and beta < 1.0:
		beta = 1.0 / beta
	return A * ( beta / (2.0 * alpha * gamma_func(1.0 / beta)) ) * np.exp(-(np.abs(x - mu)/alpha)**beta)


def log_generalized_normal(x, A, mu, alpha, beta):
	A = A**2
	y = generalized_normal(x, A, mu, alpha, beta)
	if (y <= 0.0).any():
		#y[y<=0] = -y[y<=0] + 1.0
		y[y<=0] = (y[y<=0] + 0.0001)**2
	return np.log10(y)


def add_text(fig, ax, textstr):
	#props = dict(boxstyle='round', facecolor='white', alpha=0.25)
	props = dict(edgecolor='none', facecolor='none')
	ax.text(0.02, 0.16, textstr, transform=ax.transAxes, fontsize=14, \
	        verticalalignment='top', bbox=props)
	return fig, ax


def rebin_stats(ntestbins, x, x_min=None, x_max=None, use_log=False):
	if use_log:
		xbins = np.logspace(np.log10(x_min), np.log10(x_max), num=ntestbins+1)
	else:
		xbins = np.linspace(x_min, x_max, num=ntestbins+1)

	hist, bin_edges = np.histogram(x, bins=xbins)

	bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
	if ignore_central_bin:
		mask = (np.abs(bin_centers) > 0.000001)
		bin_centers = bin_centers[mask]
		hist = hist[mask]
	#p0 = [hist.max(), 0.0, 0.2]
	p0 = [hist.max(), hist.mean(), hist.std(), stats.skew(hist), stats.kurtosis(hist)]
	hist[hist==0] = 1  #fix devide by zero error
	try:
		if poisson_weight:
			coeffs, var_matrix = curve_fit(gaussian, bin_centers, hist, p0=p0, sigma=(np.sqrt(hist)/hist))
		else:
			coeffs, var_matrix = curve_fit(gaussian, bin_centers, hist, p0=p0)
	except RuntimeError:
		print '*******curve_fit failed!'
		return np.nan, np.nan

	mean, stdev = coeffs[1], coeffs[2]
	return mean, stdev


nbins = 35
#nbins = 25
#nbins = -1
nbins_min = 15
nbins_max = 200
#nbins_max = 200
nfitpoints = 100
remove_nonfit_halos = False
global_filter_halos = True
column_filter_halos = True
use_klypin = False
label_axes = True
ignore_central_bin = False
save_stats = True
bin_test = False
poisson_weight = True
fit_in_log = True
draw_data_fit = False
mass_quartiles = False
prevent_small_shape_param = False
add_residuals_panel = False
perc_diff = True

hist_normed = False
hist_cumulative = False
ylog = False
ylabel= 'Number of Halos'

#                                     v           v           v
fractions = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
diff_vals = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00, 2.00, 4.00]
#                         ^           ^           ^           ^

#nhalos = 100
nhalos = None
sort_col = 9

#lpt_log_cols  = [ 9, 23, 31, 47,  51, 57]
#za_log_cols   = [10, 24, 32, 48,  52, 58]
#lpt_cols      = [17, 77, 91, 93, 97,  99, 107, 111, -4, -2]
#za_cols       = [18, 78, 92, 94, 98, 100, 108, 112, -3, -1]

lpt_log_cols  = []
za_log_cols   = []
#lpt_cols      = [-4, 47, 91, 107, 111]
#za_cols       = [-3, 48, 92, 108, 112]
#lpt_cols      = [-4, 31, 47, 91, 107, 111]
#za_cols       = [-3, 32, 48, 92, 108, 112]
#lpt_cols      = [-4, 31, 47, 91, 111]
#za_cols       = [-3, 32, 48, 92, 112]
lpt_cols      = [-4, 47, 91, 93, 107]
za_cols       = [-3, 48, 92, 94, 108]
# conentration, mass, x_off, v_off, T/|U|

fancy_log_x_labels = []
#fancy_x_labels = [r"$\mathrm{\frac{c_{2LPT} - c_{ZA}}{c_{avg}}}$",
#                  r"$\mathrm{\frac{\rho_{0, 2LPT} - \rho_{0, ZA}}{\rho_{0, avg}}}$",
#                  r"$\mathrm{\frac{M_{vir, 2LPT} - M_{vir, ZA}}{M_{vir, avg}}}$",
#                  r"$\mathrm{\frac{X_{off, 2LPT} - X_{off, ZA}}{X_{off, avg}}}$",
#                  r"$\mathrm{\frac{N_{subs, 2LPT} - N_{subs, ZA}}{N_{subs, avg}}}$"]

fancy_x_labels = [r"$\mathrm{\frac{c_{2LPT} - c_{ZA}}{c_{avg}}}$",
                  r"$\mathrm{\frac{M_{vir, 2LPT} - M_{vir, ZA}}{M_{vir, avg}}}$",
                  r"$\mathrm{\frac{X_{off, 2LPT} - X_{off, ZA}}{X_{off, avg}}}$",
                  r"$\mathrm{\frac{V_{off, 2LPT} - V_{off, ZA}}{V_{off, avg}}}$",
                  r"$\mathrm{\frac{(T/|U|)_{2LPT} - (T/|U|)_{ZA}}{(T/|U|)_{avg}}}$"]

fancy_y_label  = r"$\mathrm{N_{halos}}$"

Rv1_col = 53
Rv2_col = 54
Rs1_col = 55
Rs2_col = 56

c_lpt_col = 17
c_za_col  = 18


# c_2lpt, c_za, chi2_2lpt, chi2_za
#lt_cols = [17, 18, 37, 38]
#lt_vals = [100.0, 100.0, 10.0, 10.0]
lt_cols = [17, 18]
lt_vals = [100.0, 100.0]

# c_2lpt, c_za, rho_0_2lpt, rho_0_za, chi2_2lpt, chi2_za
#gt_cols = [17, 18, 31, 32, 37, 38]
#gt_vals = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
gt_cols = [17, 18, 31, 32]
gt_vals = [1.0, 1.0, 0.0, 0.0]

eq_cols = []
eq_vals = []

ne_cols = []
ne_vals = []


# global filters
glob_lt_cols = [] glob_lt_vals = []

glob_gt_cols = [4, 5]
glob_gt_vals = [100, 100]

glob_eq_cols = [109, 110]
glob_eq_vals = [-1, -1]

glob_ne_cols = []
glob_ne_vals = []



use_alt_frac = True
alt_start_frac = 0.75
alt_end_frac = 1.0

#x_lim = 0.5
#x_lim = 1.0
x_lim = None

bad_halo_pairs = None
#bad_halo_pairs = [9, 28, 39, 51, 59, 95]

perc_diff_base = 'plots/perc_diff_'
#statsfile = 'plots/stats.txt'
stats_base = 'plots/stats_'
stats_ext  = '.txt'
plot_base = 'plots/hist_'
plot_ext  = '.eps'

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


if __name__ == '__main__':
	main()

