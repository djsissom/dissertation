#!/usr/bin/env python

import sys
import bgc2
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import patheffects
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy.stats import ks_2samp
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from ipdb import set_trace


#### Note:  only run one box pair at a time.
#### ex:    ./compare.py /crossmatch_dir/halos.dat /2lpt_dir/halos_0.*.bgc2 /za_dir/halos_0.*.bgc2

def main():
	crossmatched_halo_file, bgc2_2lpt_files, bgc2_za_files = parse_args(sys.argv[1:])

	header, halos = read_files(crossmatched_halo_file, header_line = 3)

	bgc2_2lpt_header, bgc2_2lpt_halos, bgc2_2lpt_particles = get_bgc2_data(bgc2_2lpt_files)
	bgc2_za_header,   bgc2_za_halos,   bgc2_za_particles   = get_bgc2_data(bgc2_za_files)

	header = np.asarray(header)
	bgc2_2lpt_halos, bgc2_za_halos = map(np.asarray, (bgc2_2lpt_halos, bgc2_za_halos))

	if sort_col != None:
		halos = sort_by_column(halos, sort_col)
	if remove_nonfit_halos:
		halos = remove_nans(halos)
	if global_filter_halos:
		halos = filter_halos(halos)
	if (nhalos != None) or (nhalos != 0):
		#halos = halos[:nhalos]
		halos = halos[[0,70]]        ############################# hard coded for the moment
		#halos = halos[10000:10050]
	
	header, halos = add_c_columns(header, halos)
	header = reduce_header(header)

	for i, halo_pair in enumerate(halos):
		make_plot(i, header, halo_pair, bgc2_2lpt_halos, bgc2_za_halos, \
		          bgc2_2lpt_particles, bgc2_za_particles)

	print 'Finished all plots.'


def parse_args(args):
	crossmatched_halo_file = args[0]
	if len(args[1:]) % 2 != 0.0:
		print 'Must call with even number of bgc2 files...exiting.'
		sys.exit(-1)
	bgc2_files = args[1:]
	bgc2_2lpt_files = bgc2_files[:len(bgc2_files)/2]
	bgc2_za_files = bgc2_files[len(bgc2_files)/2:]
	return crossmatched_halo_file, bgc2_2lpt_files, bgc2_za_files


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


def get_bgc2_data(bgc2_files):
	header    = None
	halos     = None
	particles = None
	for bgc2_file in bgc2_files:
		print 'Reading file %s...' % (bgc2_file)
		tmp_header, tmp_halos, tmp_particles = bgc2.read_bgc2(bgc2_file)
		if header == None:
			header = tmp_header
			halos = tmp_halos
			particles = tmp_particles
		else:
			halos = np.append(halos, tmp_halos, axis=0)
			particles = np.append(particles, tmp_particles, axis=0)
	print 'Finished reading bgc2 files.'
	return header, halos, particles


def sort_by_column(halos, col):
	print 'Sorting halos...'
	mask = np.argsort(halos[:, col])
	mask = mask[::-1]
	halos = halos[mask]
	return halos


def remove_nans(halos):
	print 'Removing NaNs...'
	halos = halos[halos[:,c_2lpt_col] != -9999]
	halos = halos[np.isfinite(halos[:,c_2lpt_col])]
	halos = halos[np.isfinite(halos[:,c_za_col])]
	return halos


def filter_halos(halos):
	print 'Filtering data...'
	for col, val in zip(lt_cols, lt_vals):
		halos = halos[halos[:, col] <= val]
	for col, val in zip(gt_cols, gt_vals):
		halos = halos[halos[:, col] >= val]
	for col, val in zip(eq_cols, eq_vals):
		halos = halos[halos[:, col] == val]
	for col, val in zip(ne_cols, ne_vals):
		halos = halos[halos[:, col] != val]
	return halos


def add_c_columns(header, halos):
	c1_rockstar = halos[:, Rv1_col] / halos[:, Rs1_col]
	c2_rockstar = halos[:, Rv2_col] / halos[:, Rs2_col]
	halos = np.column_stack((halos, c1_rockstar, c2_rockstar))
	header = np.append(header, 'c_rockstar')
	header = np.append(header, 'c_rockstar')
	return header, halos


def reduce_header(header):
	header_2lpt = header[print_cols_2lpt]
	header_za   = header[print_cols_za]
	if (header_2lpt == header_za).all():
		header = header_2lpt
	else:
		print 'column mismatch... exiting'
		set_trace()
		sys.exit(123)
	return header


def make_plot(itteration, header, halo_pair, bgc2_halos_2lpt, bgc2_halos_za, \
              bgc2_particles_2lpt, bgc2_particles_za):
	id_2lpt = halo_pair[id_col_2lpt]
	id_za   = halo_pair[id_col_za]
	properties_2lpt = halo_pair[print_cols_2lpt]
	properties_za   = halo_pair[print_cols_za]

	# find 2lpt and za halo from id
	halo_index_2lpt = np.where(bgc2_halos_2lpt[:, halo_id_col] == id_2lpt)[0][0]
	halo_index_za   = np.where(bgc2_halos_za[:, halo_id_col]   == id_za)[0][0]

	bgc2_halos_2lpt = bgc2_halos_2lpt[halo_index_2lpt]
	bgc2_halos_za   = bgc2_halos_za[halo_index_za]

	# convert particles to numpy arrays
	bgc2_particles_2lpt = np.asarray(bgc2_particles_2lpt[halo_index_2lpt])
	bgc2_particles_za   = np.asarray(bgc2_particles_za[halo_index_za])

	# make density profiles
	r_2lpt, rho_2lpt, rho_err_2lpt, r_vir_2lpt = density_profile(bgc2_halos_2lpt, bgc2_particles_2lpt)
	r_za,   rho_za,   rho_err_za,   r_vir_za   = density_profile(bgc2_halos_za, bgc2_particles_za)

	# fit density profiles
	nfw_r_2lpt, nfw_rho_2lpt, r_s_2lpt = fit_profile( r_2lpt / r_vir_2lpt, rho_2lpt / rho_2lpt.max(), err = rho_err_2lpt / rho_2lpt.max() )
	nfw_r_za  , nfw_rho_za  , r_s_za   = fit_profile( r_za   / r_vir_za,   rho_za   / rho_za.max(),   err = rho_err_za   / rho_za.max()   )

	# de-normalize values
	nfw_r_2lpt = nfw_r_2lpt * r_vir_2lpt
	nfw_r_za   = nfw_r_za   * r_vir_za
	nfw_rho_2lpt = nfw_rho_2lpt * rho_2lpt.max()
	nfw_rho_za   = nfw_rho_za   * rho_za.max()
	r_s_2lpt = r_s_2lpt * r_vir_2lpt
	r_s_za   = r_s_za   * r_vir_za

	# find center of halos and plot limit
	halo_pos_2lpt = bgc2_halos_2lpt[:,halo_pos_cols] * dist_scale
	halo_pos_za   = bgc2_halos_za[:,halo_pos_cols] * dist_scale
	particle_pos_2lpt = bgc2_particles_2lpt[:,particle_pos_cols] * dist_scale
	particle_pos_za   = bgc2_particles_za[:,particle_pos_cols] * dist_scale

	if wrap_box:
		for i in range(3):
			if abs(halo_pos_2lpt[i] - halo_pos_za[i]) > box_size / 2.0:
				print "#################################### wrapping halo ####################################"
				if (halo_pos_2lpt[i] > halo_pos_za[i]):
					halo_pos_za[i] += box_size
					particle_pos_za[:,i] += box_size
				if (halo_pos_2lpt[i] < halo_pos_za[i]):
					halo_pos_2lpt[i] += box_size
					particle_pos_2lpt[:,i] += box_size
				else:
					print "error in wrapping"
					sys.exit()

	center = (halo_pos_2lpt + halo_pos_za) / 2.0
	halo_pos_2lpt = halo_pos_2lpt - center
	halo_pos_za   = halo_pos_za   - center
	particle_pos_2lpt = particle_pos_2lpt - center
	particle_pos_za   = particle_pos_za   - center

	if zoom_projections:
		plot_lim = zoom_scale
	else:
		plot_lim = np.append(particle_pos_2lpt, particle_pos_za).max()


	r_vir_2lpt = bgc2_halos_2lpt[halo_r_col] * dist_scale
	r_vir_za   = bgc2_halos_za[halo_r_col] * dist_scale

	if make_stats:
		print 'generating plot...'
		fig = plt.figure(figsize = (9.0, 6.0))
		fig = make_projections(fig, 221, halo_pos_2lpt, halo_pos_za, particle_pos_2lpt, particle_pos_za, \
							   r_vir_2lpt, r_vir_za, plot_lim)
		ax = fig.add_subplot(223)
		ax = draw_density_profile(ax, r_2lpt, rho_2lpt, err=rho_err_2lpt, color='blue', label='2lpt')
		ax = draw_density_profile(ax, r_za,   rho_za,   err=rho_err_za,   color='red',  label='za')

		ax = fig.add_subplot(122)
		ax = draw_parameters(ax, header, properties_2lpt, properties_za)

		fig.tight_layout()
		plot_name = "%s%0.3d_(%d,%d)%s" % (plot_base, itteration, id_2lpt, id_za, plot_ext)
		plt.savefig(plot_name, bbox_inches='tight')
		print 'finished plot ' + plot_name

	if make_projection:
		print 'generating density projection plot...'
		fig = plt.figure(figsize = (9.0, 6.0))

		if label_projection:
			ax = fig.add_subplot(111, aspect=2.0/3.2)
			ax = hide_axes(ax)
			ax.set_xlabel(proj_xlabel)
			ax.set_ylabel(proj_ylabel)

		fig = make_projections(fig, 111, halo_pos_2lpt, halo_pos_za, particle_pos_2lpt, particle_pos_za, \
							   r_vir_2lpt, r_vir_za, plot_lim)
		fig.tight_layout()
		plot_name = "%s%0.3d_(%d,%d)%s%s" % (plot_base, itteration, id_2lpt, id_za, proj_name, plot_ext)
		plt.savefig(plot_name, bbox_inches='tight')
		print 'finished density projection plot ' + plot_name

	if make_density_profile:
		print 'generating density profile plot...'
		fig = plt.figure(figsize = (9.0, 12.0))

		if label_projection:
			ax = fig.add_subplot(211, aspect=2.0/3.2)
			ax = hide_axes(ax)
			ax.set_xlabel(proj_xlabel)
			ax.set_ylabel(proj_ylabel)

		fig = make_projections(fig, 211, halo_pos_2lpt, halo_pos_za, particle_pos_2lpt, particle_pos_za, \
							   r_vir_2lpt, r_vir_za, plot_lim)

		ax = fig.add_subplot(212)
		ax = hide_axes(ax)
		ax.set_xlabel(prof_xlabel)
		ax.set_ylabel(prof_ylabel)

		#grid = ImageGrid(fig, 212, nrows_ncols=(2,1), axes_pad=0.24)

		ax1 = fig.add_subplot(413)
		ax1 = draw_density_profile(ax1, r_2lpt, rho_2lpt, err=rho_err_2lpt, color='blue')
		ax1 = draw_nfw_profile(ax1, nfw_r_2lpt, nfw_rho_2lpt, R_s=r_s_2lpt, color='red')

		ax2 = fig.add_subplot(414)
		ax2 = draw_density_profile(ax2, r_za, rho_za, err=rho_err_za, color='blue')
		ax2 = draw_nfw_profile(ax2, nfw_r_za, nfw_rho_za, R_s=r_s_za, color='red')

		if equal_profile_axes:
			ymin = min(ax1.get_ylim()[0], ax2.get_ylim()[0])
			ymax = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
			ax1.set_ylim(ymin, ymax)
			ax2.set_ylim(ymin, ymax)

			xmin = min(ax1.get_xlim()[0], ax2.get_xlim()[0])
			xmax = max(ax1.get_xlim()[1], ax2.get_xlim()[1])
			ax1.set_xlim(xmin, xmax)
			ax2.set_xlim(xmin, xmax)

		if print_text:
			ax1.text(0.95, 0.85, '2LPT', color='black', horizontalalignment='right', verticalalignment='center', transform=ax1.transAxes)
			ax2.text(0.95, 0.85, 'ZA',   color='black', horizontalalignment='right', verticalalignment='center', transform=ax2.transAxes)


		#fig.tight_layout()
		plot_name = "%s%0.3d_(%d,%d)%s%s" % (plot_base, itteration, id_2lpt, id_za, dens_name, plot_ext)
		plt.savefig(plot_name, bbox_inches='tight')
		print 'finished density profile plot ' + plot_name



def density_profile(halo, particles):
	r_vir = halo[halo_r_col] * dist_scale
	halo_pos = halo[halo_pos_cols] * dist_scale
	#mass = np.ones(particles.shape[0]) * common_mass * mass_scale
	mass = particles[:,particle_mass_col] * mass_scale
	pos = particles[:,particle_pos_cols] * dist_scale
	pos = pos - halo_pos
	
	r_bins, rho, rho_err = calc_density_profile(mass, pos)
	mid_bins = 10.0**(0.5 * (np.log10(r_bins[1:]) + np.log10(r_bins[:-1])))
	
	# Don't pass NaN's to fitting routine
	rho_err_nonan = np.copy(rho_err)
	nan_check = np.isnan(rho_err_nonan)
	for i in range(len(rho_err_nonan)):
		if (mid_bins[i] < res_limit) or (nan_check[i] == True):
			rho_err_nonan[i] = 1.0e10

	return mid_bins, rho, rho_err, r_vir


def calc_density_profile(mass, pos):
	r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
	max_r = r.max()
	min_r = res_limit
	log_range = np.log10(max_r) - np.log10(min_r)
	local_nbins = float(nbins + 1)
	while True:
		bins = np.arange(local_nbins)
		bins = max_r * 10.0**(log_range * bins / (local_nbins-1.0) - log_range)
		bin_mass, r_bins = np.histogram(r, bins, weights=mass)
		if (bin_mass == 0.0).any():
			local_nbins -= 1
			continue
		else:
			break
	rho = bin_mass / (sphere_vol(r_bins[1:]) - sphere_vol(r_bins[:-1]))
	N_bin, blah = np.histogram(r, bins)
	rho_err = poisson_error(N_bin) * rho
	return r_bins, rho, rho_err


def sphere_vol(r):
	volume = (4.0 / 3.0) * np.pi * r**3
	return volume


def poisson_error(N):
	err = np.sqrt(N) / N
	return err


def fit_profile(r, rho, err=None, R_vir=None):
	popt, pcov = curve_fit(nfw_profile, r, rho, sigma=err, p0=[0.1, 1.0])
	R_s, rho_0 = popt[0], popt[1]
	nfw_r = np.linspace(r[0], r[-1], nfit)
	nfw_rho = nfw_profile(nfw_r, R_s, rho_0)
	return nfw_r, nfw_rho, R_s


def nfw_profile(r, R_s, rho_0):
	if R_s >= 1.0:
		return (R_s - 1.0) * np.exp(r) + rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**2)
	return rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**2)


def filter_column(x, x_col):
	print 'Filtering data...'
	x = x[x != -9999]
	if x_col in lt_cols:
		val = lt_vals[lt_cols.index(x_col)]
		x = x[x <= val]
	if x_col in gt_cols:
		val = gt_vals[gt_cols.index(x_col)]
		x = x[x >= val]
	if x_col in eq_cols:
		val = eq_vals[eq_cols.index(x_col)]
		x = x[x == val]
	if x_col in ne_cols:
		val = ne_vals[ne_cols.index(x_col)]
		x = x[x != val]
	return x


def draw_hist(fig, ax, x, x_min=None, x_max=None, use_log=False, color=None, label=None):
	if use_log:
		xbins = np.logspace(np.log10(x_min), np.log10(x_max), num=nbins+1)
		ax.set_xscale('log')
	else:
		xbins = np.linspace(x_min, x_max, num=nbins+1)

	n, bins, patches = ax.hist(x, bins=xbins, histtype='step', log=ylog, color=color, label=label)
	return fig, ax, n, bins, patches


def add_text(fig, ax, textstr):
	props = dict(boxstyle='round', facecolor='white', alpha=0.7)
	ax.text(0.02, 0.08, textstr, transform=ax.transAxes, fontsize=14, \
	        verticalalignment='top', bbox=props)
	return fig, ax


def make_projections(fig, position, halo_pos1, halo_pos2, pos1, pos2, r_vir1, r_vir2, plot_lim):
	#grid = ImageGrid(fig, position, nrows_ncols=(2,3), axes_pad=0.05, cbar_mode='single')
	grid = ImageGrid(fig, position, nrows_ncols=(2,3), axes_pad=0.12, cbar_mode='single')
	for i, (x, y, hx, hy, r) in enumerate(zip( \
			(pos1[:,0], pos1[:,0], pos1[:,1], pos2[:,0], pos2[:,0], pos2[:,1]), \
			(pos1[:,1], pos1[:,2], pos1[:,2], pos2[:,1], pos2[:,2], pos2[:,2]), \
			(halo_pos1[0], halo_pos1[0], halo_pos1[1], halo_pos2[0], halo_pos2[0], halo_pos2[1]), \
			(halo_pos1[1], halo_pos1[2], halo_pos1[2], halo_pos2[1], halo_pos2[2], halo_pos2[2]), \
			(r_vir1, r_vir1, r_vir1, r_vir2, r_vir2, r_vir2))):
		ax = grid[i]
		draw_projection(ax, x, y, hx, hy, r, plot_lim)
		if print_text:
			if i == 0:
				ax.text(0.05, 0.12, '2LPT', color='white', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=3, foreground='black')])
			if i == 3:
				ax.text(0.05, 0.12, 'ZA',   color='white', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=3, foreground='black')])

			if i == 0:
				ax.text(0.95, 0.88, 'XY', color='white', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=3, foreground='black')])
			if i == 1:
				ax.text(0.95, 0.88, 'XZ', color='white', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=3, foreground='black')])
			if i == 2:
				ax.text(0.95, 0.88, 'YZ', color='white', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=3, foreground='black')])
	return fig


def draw_projection(ax, x, y, hx, hy, r, plot_lim):
	limits = [[-plot_lim, plot_lim], [-plot_lim, plot_lim]]
	z, xedges, yedges = np.histogram2d(x, y, bins=npixels, range=limits)
	if log_scale_projections:
		z[z<1.0] = 0.5
		#z = np.log10(z)
		#z = np.log10(z)
		#z[np.isinf(z)] = -0.1
		plot_norm = mpl.colors.LogNorm(vmin = 1, vmax = z.max(), clip=True)
		#plot_norm = None
	else:
		plot_norm = None
	if extra_smoothing:
		z = gaussian_filter(z, smoothing_radius)
	im = ax.imshow(z.T, extent=(-plot_lim, plot_lim, -plot_lim, plot_lim), \
			interpolation='gaussian', origin='lower', cmap=colormap, norm=plot_norm)
			#interpolation='gaussian', origin='lower', cmap=colormap)
	ax.locator_params(nbins=6)
	if draw_circle:
		ax.add_patch(Circle((hx, hy), r, fc="None", ec="black", lw=1))
	if draw_contours:
		x_midpoints = (xedges[:-1] + xedges[1:]) / 2.0
		y_midpoints = (yedges[:-1] + yedges[1:]) / 2.0
		X, Y = np.meshgrid(x_midpoints, y_midpoints)
		ax.contour(X, Y, z.T, 2, colors='black', linewidths=4)
		ax.contour(X, Y, z.T, 2, colors='white', linewidths=2)
	if label_colorbar:
		if log_scale_projections:
			log_format = mpl.ticker.LogFormatterMathtext(10, labelOnlyBase=False)
			ax.cax.colorbar(im, format=log_format)
		else:
			ax.cax.colorbar(im)
	else:
		bar = ax.cax.colorbar(im, ticks=[])
		bar.ax.set_yticklabels([])
		#plt.setp(bar.ax.get_yticklabels(), visible=False)


def draw_density_profile(ax, r, rho, err=None, color='black', label=None):
	im = ax.loglog(r, rho, linestyle='steps-mid-', color=color, label=label)
	line1 = ax.axvline(res_limit, color='black', linestyle=':')
	ax.set_xlim(r[0] - (r[1]-r[0]), r[-1] + (r[-1]-r[-2]))
	#ax.set_xlabel(xlabel_prof)
	#ax.set_ylabel(ylabel_prof)
	if err != None:
		err_bars = ax.errorbar(r, rho, yerr=err,linestyle='None', color=color)
	if label != None:
		ax.legend(fontsize='x-small')
	return ax


def draw_nfw_profile(ax, r, rho, R_s=None, color='black'):
	ax.loglog(r, rho, linestyle='-', color=color)
	if R_s != None:
		line = ax.axvline(R_s, color='purple', linestyle='-.')
	return ax


def draw_parameters(ax, header, params1, params2):
	strlen = 12
	header  = [str(item)[:strlen] for item in header]
	params1 = [str(item)[:strlen] for item in params1]
	params2 = [str(item)[:strlen] for item in params2]
	header.insert(0, 'simulation')
	params1.insert(0, '-- 2lpt --')
	params2.insert(0, '--- za ---')
	header = '\n'.join(header)
	params1 = '\n'.join(params1)
	params2 = '\n'.join(params2)
	ax.text(0.05, 0.5, header,  horizontalalignment="left", verticalalignment="center", transform=ax.transAxes)
	ax.text(0.40, 0.5, params1, horizontalalignment="left", verticalalignment="center", transform=ax.transAxes)
	ax.text(0.75, 0.5, params2, horizontalalignment="left", verticalalignment="center", transform=ax.transAxes)
	ax.axis('off')
	return ax


def hide_axes(ax):
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
	return ax




nhalos = 1
sort_col = 9  # density_profile 2lpt halo mass
#sort_col = 47  # rockstar 2lpt halo mass (M200c)

nbins = 40
nfit = 100
npixels = 30
#npixels = 100
smoothing_radius = 0.9
remove_nonfit_halos = True
global_filter_halos = True
column_filter_halos = True
log_scale_projections = True
wrap_box = False
label_colorbar = False
label_projection = True
zoom_projections = True
zoom_scale = 18.0 # kpc
draw_circle = False
draw_contours = True
extra_smoothing = True
label_proj = True
label_2lpt_za = True
equal_profile_axes = True
print_text = True

box_size = 10000.0 # kpc

id_col_2lpt = 0
id_col_za   = 1

print_cols_2lpt  = [43, 57, 6,  9, 17, 23, 31, 47, 51, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 91, 93, 97,  99, 101, 103, 105, 107, 111, 163, 201, -2]
print_cols_za    = [44, 58, 6, 10, 18, 24, 32, 48, 52, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 92, 94, 98, 100, 102, 104, 106, 108, 112, 164, 202, -1]

Rv1_col = 53
Rv2_col = 54
Rs1_col = 55
Rs2_col = 56

c_2lpt_col = 17
c_za_col  = 18

# c_2lpt, c_za, chi2_2lpt, chi2_za
lt_cols = [17, 18, 37, 38]
lt_vals = [100.0, 100.0, 10.0, 10.0]

# c_2lpt, c_za, rho_0_2lpt, rho_0_za, chi2_2lpt, chi2_za
gt_cols = [17, 18, 31, 32, 37, 38]
gt_vals = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

eq_cols = []
eq_vals = []

ne_cols = []
ne_vals = []

# bgc2 halo array columns
halo_id_col   = 0
halo_r_col    = 4
halo_mass_col = 5
halo_pos_cols = [6,7,8]

# bgc2 particle array columns
particle_mass_col = 0
particle_pos_cols = [1,2,3]
particle_vel_cols = [4,5,6]

mass_scale = 1.0
common_mass = 5.33423e5
dist_scale = 1.0e3
res_limit = 0.5  #changed from 4.0 to 0.5 to match density_profile.py <-- maybe check why it was 4.0?
nfit = 500

dist_units = 'kpc'
#xlabel_proj = [r'X Position (%s h$^{-1}$)' % (dist_units), r'X Position (%s h$^{-1}$)' % (dist_units), r'Y Position (%s h$^{-1}$)' % (dist_units)]
#ylabel_proj = [r'Y Position (%s h$^{-1}$)' % (dist_units), r'Z Position (%s h$^{-1}$)' % (dist_units), r'Z Position (%s h$^{-1}$)' % (dist_units)]
proj_xlabel = r'Position (kpc h$^{-1}$)'
proj_ylabel = r'Position (kpc h$^{-1}$)'
prof_xlabel = r'Radius (%s h$^{-1}$)' % (dist_units)
prof_ylabel = r'Density (M$_{\odot}$ %s$^{-3}$ h$^{2}$)' % (dist_units)

#colormap = 'ocean_r'
colormap = 'rainbow'
plot_base = 'plots/halo_pair_'
proj_name = '_proj'
dens_name = '_dens'
plot_ext  = '.eps'

make_stats = False
make_projection = False
make_density_profile = True

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

