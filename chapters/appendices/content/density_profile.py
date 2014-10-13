#!/usr/bin/env python

import sys
import bgc2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from scipy.stats import chisquare

#read_mode = 'ascii2'
read_mode = 'bgc2'

if read_mode == 'bgc2':
  use_bgc2 = True
  use_all = False
  individual_masses = False
  halo_id = 146289
  nbins = 50
  nfit = 500
  ooms = 3.0
  mass_scale = 1.0
  common_mass = 5.33423e5
  dist_scale = 1.0e3
  #res_limit = 0.488
  #res_limit = 4.0
  res_limit = 0.5
  #res_limit = 10.0
  draw_frac = 0.1
  tick_base_major = 100.0
  tick_base_minor = 10.0
  find_com = False
elif read_mode == 'ascii':
  use_bgc2 = False
  use_all = True
  individual_masses = True
  halo_id = 0
  nbins = 100
  nfit = 500
  ooms = 5.0
  mass_scale = 1.0e12
  dist_scale = 200.0
  res_limit = 1.0e-2
  draw_frac = 2.0e-4
  tick_base_major = 80.0
  tick_base_minor = 20.0
  find_com = True
elif read_mode == 'ascii2':
  use_bgc2 = False
  use_all = True
  individual_masses = True
  halo_id = 0
  nbins = 100
  nfit = 500
  ooms = 3.5
  mass_scale = 1.0e10
  dist_scale = 1.0
  #res_limit = 3.0e-1
  res_limit = 1.0
  draw_frac = 1.0e-2
  tick_base_major = 200.0
  tick_base_minor = 40.0
  find_com = True
else:
  sys.exit(98712)

#outfile = 'asciitest_halo_properties.txt'
outfile = 'density_profile_halos.dat'
comfile = 'center_of_mass.txt'

make_plot = False
#make_plot = True
draw_density = True
#plot_base = 'asciitest_density_profile.fig.'
plot_base = 'figure_'
plot_ext = '.eps'
dist_units = 'kpc'
xlabel_proj = [r'X Position (%s h$^{-1}$)' % (dist_units), r'X Position (%s h$^{-1}$)' % (dist_units), r'Y Position (%s h$^{-1}$)' % (dist_units)]
ylabel_proj = [r'Y Position (%s h$^{-1}$)' % (dist_units), r'Z Position (%s h$^{-1}$)' % (dist_units), r'Z Position (%s h$^{-1}$)' % (dist_units)]
xlabel_prof = r'Radius (%s h$^{-1}$)' % (dist_units)
ylabel_prof = r'Density (M$_{\odot}$ %s$^{-3}$ h$^{2}$)' % (dist_units)
npixels = 50

#common_mass = 1.0e-7
#common_mass = 1.0e5
mass_col = 0
pos_cols = (1,2,3)
vel_cols = (4,5,6)
halo_id_col = 0

grav_const = 4.3e-6 # kpc M_sol^-1 (km/s)^2


def read_files(files):
  data = 0
  for file in files:
    print 'Reading file %s...' % (file)
    if data == 0:
      data = np.genfromtxt(file, comments='#')
    else:
      data = np.append(data, np.genfromtxt(file, comments='#'), axis=0)
  print 'Finished reading files.'
  return data


def my_chisq(ydata,ymod,deg=2,sd=None):  
  """  
Returns the reduced chi-square error statistic for an arbitrary model,   
chisq/nu, where nu is the number of degrees of freedom. If individual   
standard deviations (array sd) are supplied, then the chi-square error   
statistic is computed as the sum of squared errors divided by the standard   
deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.  

ydata,ymod,sd assumed to be Numpy arrays. deg integer.  

Usage:  
>>> chisq=redchisqg(ydata,ymod,n,sd)  
where  
ydata : data  
ymod : model evaluated at the same x points as ydata  
n : number of free parameters in the model  
sd : uncertainties in ydata  

Rodrigo Nemmen  
http://goo.gl/8S1Oo  
   """  
  # Chi-square statistic  
  if sd==None:
    chisq=np.sum((ydata-ymod)**2)  
  else:
    chisq=np.sum( ((ydata-ymod)/sd)**2 )  

  # Number of degrees of freedom assuming 2 free parameters  
  nu=ydata.size-1-deg  
  return chisq/nu       


def calc_m_enclosed(mass, pos):
  r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
  r = np.sort(r)
  first_good_bin = 0
  for i in range(len(r)):
    if r[i] > res_limit:
      first_good_bin = i
      break
  print 'r1 =', r[first_good_bin-1]
  print 'r2 =', r[first_good_bin]
  print 'r3 =', r[first_good_bin+1]
  m_extra = mass[0] * first_good_bin
  r = r[first_good_bin:]
  #m_enclosed = np.zeros(len(r))
  #for i in range(len(r)):
  #  m_enclosed[i] = mass[0] * (i + 1.0)
  m_enclosed = (np.arange(len(r)) + 1.0) * mass[0] + m_extra
  return r, m_enclosed


def calc_density_profile(mass, pos):
  r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
  max_r = r.max()
  #min_r = max_r / 10**ooms
  min_r = res_limit
  log_range = np.log10(max_r) - np.log10(min_r)
  
  #global nbins
  local_nbins = float(nbins + 1)
  #nbins = len(r) / 1000
  while True:
    bins = np.arange(local_nbins)
    bins = max_r * 10.0**(log_range * bins / (local_nbins-1.0) - log_range)
    bin_mass, r_bins = np.histogram(r, bins, weights=mass)
    if (bin_mass == 0.0).any():
      local_nbins -= 1
      continue
    else:
      break

  #print 'Binning particles using bin edges of \n', r_bins

  rho = bin_mass / (sphere_vol(r_bins[1:]) - sphere_vol(r_bins[:-1]))
  
  N_bin, blah = np.histogram(r, bins)
  rho_err = poisson_error(N_bin) * rho

  return r_bins, rho, rho_err


def logbin(pos):
  r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
  max_r = r.max()
  min_r = max_r / 10**ooms
  log_range = np.log10(max_r) - np.log10(min_r)
  
  global nbins
  nbins = float(nbins + 1)
  bins = np.arange(nbins)
  bins = max_r * 10.0**(log_range * bins / (nbins-1.0) - log_range)

  hist, bin_edges = np.histogram(r, bins)
  #print 'Binning particles using bin edges of \n', bin_edges
  return hist, bin_edges


def poisson_error(N):
  err = np.sqrt(N) / N
  return err


def sphere_vol(r):
  volume = (4.0 / 3.0) * np.pi * r**3
  return volume


def get_rho_0(R_s, R_vir):
  H = 70.0e-3 # km s^-1 kpc^-1
  G = 4.3e-6 # kpc M_sol^-1 (km/s)^2
  rho_crit = 3.0 * H**2 / (8.0 * np.pi * G)

  v = 178
  c = R_vir / R_s
  g = 1.0 / (np.log(1.0+c) - c/(1.0+c))
  delta_char = v * c**3 * g / 3.0

  return rho_crit * delta_char


def nfw_fit_rho0(r, R_s, rho_0):
  if R_s >= 1.0:
    return (R_s - 1.0) * np.exp(r) + rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**2)
  return rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**2)


def nfw_fit_rho0_log(r, R_s, rho_0):
  r = 10.0**r
  R_s = 10.0**R_s
  rho_0 = 10.0**rho_0
  profile =  rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**2)
  return np.log10(profile)


def nfw_def_rho0(R_vir):
  def _nfw_def_rho0(r, R_s):
    rho_0 = get_rho_0(R_s, R_vir)
    return rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**2)
  return _nfw_def_rho0


def nfw_databin_rho0(rho_0):
  def _nfw_databin_rho0(r, R_s):
    return rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**2)
  return _nfw_databin_rho0


def dm_profile_fit_rho0_log(r, R_s, rho_0, alpha):
  r = 10.0**r
  R_s = 10.0**R_s
  rho_0 = 10.0**rho_0
  alpha = 10.0**alpha
  profile = rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**alpha)
  return np.log10(profile)


def dm_profile_fit_rho0(r, R_s, rho_0, alpha):
  return rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**alpha)


def dm_profile_def_rho0(R_vir):
  def _dm_profile_def_rho0(r, R_s, alpha):
    rho_0 = get_rho_0(R_s, R_vir)
    return rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**alpha)
  return _dm_profile_def_rho0


def dm_profile_databin_rho0(rho_0):
  def _dm_profile_databin_rho0(r, R_s, alpha):
    return rho_0 / (( r / R_s ) * ( 1.0 + r / R_s )**alpha)
  return _dm_profile_databin_rho0


def nfw_cdf(r, R_s, rho_0):
  r = 10.0**r
  R_s = 10.0**R_s
  rho_0 = 10.0**rho_0
  profile = rho_0 * R_s * (np.log(1.0 + r / R_s) - 1.0 / (1.0 + r / R_s))
  return np.log10(profile)


def nfw_cdf_nolog(r, R_s, rho_0):
  profile = rho_0 * R_s * (np.log(1.0 + r / R_s) - 1.0 / (1.0 + r / R_s))
  return profile


def mass_profile(s, c):
  g = 1.0 / (np.log(1.0 + c) - c / (1.0 + c))
  return g * (np.log(1.0 + c * s) - c * s / (1.0 + c * s))


def fit_mass_profile(s, m_enclosed, err=None, R_vir=None):
  #for i in range(len(s)):
  #  if s[i] > res_limit:
  #    first_good_bin = i
  #    break
  first_good_bin = 0

  #popt, pcov = curve_fit(nfw_cdf, np.log10(r), np.log10(m_outside), sigma=np.log10(err))
#  popt, pcov = curve_fit(nfw_cdf, np.log10(r), np.log10(m_outside))
#  popt = 10.0**popt
#  pcov = 10.0**pcov
  popt, pcov = curve_fit(mass_profile, s, m_enclosed)

  print 'fit_params =', popt
  print 'covariance =', pcov
  nfw_r = np.linspace(s[0], s[-1], nfit)
  nfw_fit = mass_profile(nfw_r, popt[0])
  chi2_fit = mass_profile(s, popt[0])

  chi2 = chisquare(np.log10(m_enclosed[first_good_bin:]), np.log10(chi2_fit[first_good_bin:]))
  chi2_nolog = chisquare(m_enclosed[first_good_bin:], chi2_fit[first_good_bin:])
  print 'chi_square =', chi2
  print 'chi_square_nolog =', chi2_nolog
  return nfw_r, nfw_fit, popt, pcov, chi2[0]


def fit_profile(r, rho, err=None, R_vir=None):
  first_good_bin = 0
#  for i in range(len(r)):
#    if r[i] > res_limit:
#      rho_0_databin = rho[i]
#      first_good_bin = i
#      break
#  print 'first_good_bin =', first_good_bin

  #--------- choose one fitting type ---------#
  #popt, pcov = curve_fit(nfw_fit_rho0, r, rho, sigma=err)
  #popt, pcov = curve_fit(nfw_def_rho0(R_vir), r, rho, p0=[10.0], sigma=err)
  #popt, pcov = curve_fit(nfw_databin_rho0(rho_0_databin), r, rho, sigma=err)
  blah = 3
  if blah == 0:
    for i in range(100):
      a = 2.0 * np.random.random() * 0.1 * r.max()
      b = 2.0 * np.random.random() * 10.0
      c = 2.0 * np.random.random() * 2.0
      try:
        popt, pcov = curve_fit(dm_profile_fit_rho0, r, rho, p0=[a,b,c], sigma=err)
      except RuntimeError:
        continue
      if (popt[0] < r.max()) and (popt[2] >= 0.0):
        break
      elif i >= 99:
        print 'no good fit found for this halo...'
#        return None, None, None, None, None
  elif blah == 1:
    #a = r.max() / 100.0
    a = 0.001
    b = rho[first_good_bin]
    c = 0.001
    #popt, pcov = curve_fit(dm_profile_fit_rho0, r, rho, sigma=err)
    print '-------------------------------------'
    print 'rho_0 before =', b
    #try:
      #popt, pcov = curve_fit(dm_profile_fit_rho0, r, rho, p0=[a,b,c], sigma=err, maxfev=1, xtol=100.0)
    popt, pcov = curve_fit(dm_profile_fit_rho0, r, rho, p0=[a,b,c], sigma=err, xtol=1.0e-1)
    #except RuntimeError:
    #  print 'just checking for now...'
    print 'rho_0 after =', popt[1]
    #sys.exit()
  elif blah == 2:
    #popt, pcov = curve_fit(dm_profile_fit_rho0_log, np.log10(r), np.log10(rho), sigma=np.log10(err))
    popt, pcov = curve_fit(nfw_fit_rho0_log, np.log10(r), np.log10(rho), sigma=np.log10(err))
    popt = 10.0**popt
    pcov = 10.0**pcov
  elif blah == 3:
    popt, pcov = curve_fit(nfw_fit_rho0, r, rho, sigma=err, p0 = [0.1, 1.0])
      
  #popt, pcov = curve_fit(dm_profile_def_rho0(R_vir), r, rho, sigma=err)
  #popt, pcov = curve_fit(dm_profile_databin_rho0(rho_0_databin), r, rho, sigma=err)
  #-------------------------------------------#

  print 'fit_params =', popt
  print 'covariance =', pcov

  nfw_r = np.linspace(r[0], r[-1], nfit)
  #--------- choose one fitting type ---------#
  nfw_fit = nfw_fit_rho0(nfw_r, popt[0], popt[1])
  #nfw_fit = nfw_def_rho0(R_vir)(nfw_r, popt[0])
  #nfw_fit = nfw_databin_rho0(rho_0_databin)(nfw_r, popt[0])
  #nfw_fit = dm_profile_fit_rho0(nfw_r, popt[0], popt[1], popt[2])
  #nfw_fit = dm_profile_def_rho0(R_vir)(nfw_r, popt[0], popt[1])
  #nfw_fit = dm_profile_databin_rho0(rho_0_databin)(nfw_r, popt[0], popt[1])
  #-------------------------------------------#
  #--------- choose one fitting type ---------#
  chi2_fit = nfw_fit_rho0(r, popt[0], popt[1])
  #chi2_fit = nfw_def_rho0(R_vir)(r, popt[0])
  #chi2_fit = nfw_databin_rho0(rho_0_databin)(r, popt[0])
  #chi2_fit = dm_profile_fit_rho0(r, popt[0], popt[1], popt[2])
  #chi2_fit = dm_profile_def_rho0(R_vir)(r, popt[0], popt[1])
  #chi2_fit = dm_profile_databin_rho0(rho_0_databin)(r, popt[0], popt[1])
  #-------------------------------------------#

  #chi2 = my_chisq(rho, chi2_fit, 2, err)
  chi2 = chisquare(rho, chi2_fit)
  print 'chi_square =', chi2
  chi2 = chi2[0]

  return nfw_r, nfw_fit, popt, pcov, chi2


def draw_projection(fig, place, plot_lim, x, y):
  ax = plt.subplot(2,3,place+1, aspect='equal')
  im = ax.plot(x, y, linestyle='', marker='.', markersize=1, markeredgecolor='blue')
  ax.set_xlabel(xlabel_proj[place])
  ax.set_ylabel(ylabel_proj[place])
  ax.set_xlim(-plot_lim, plot_lim)
  ax.set_ylim(-plot_lim, plot_lim)
  ax.xaxis.set_major_locator(MultipleLocator(tick_base_major))
  ax.xaxis.set_minor_locator(MultipleLocator(tick_base_minor))
  ax.yaxis.set_major_locator(MultipleLocator(tick_base_major))
  ax.yaxis.set_minor_locator(MultipleLocator(tick_base_minor))
  return fig


def draw_density_projection(fig, place, plot_lim, x, y):
  limits = [[-plot_lim, plot_lim], [-plot_lim, plot_lim]]
  ax = plt.subplot(2,3,place+1, aspect='equal')
  #ax.set_xlim(-plot_lim, plot_lim)
  #ax.set_ylim(-plot_lim, plot_lim)
  #im = ax.plot(x, y, linestyle='', marker='.', markersize=1, markeredgecolor='blue')
  z, xedges, yedges = np.histogram2d(x, y, bins = npixels, range = limits)
  #z = np.log10(z)
  im = ax.imshow(z.T, extent=(-plot_lim, plot_lim, -plot_lim, plot_lim), interpolation='gaussian', origin='lower')
  ax.locator_params(nbins=6)
  ax.set_xlabel(xlabel_proj[place])
  ax.set_ylabel(ylabel_proj[place])
#  ax.xaxis.set_major_locator(MultipleLocator(tick_base_major))
#  ax.xaxis.set_minor_locator(MultipleLocator(tick_base_minor))
#  ax.yaxis.set_major_locator(MultipleLocator(tick_base_major))
#  ax.yaxis.set_minor_locator(MultipleLocator(tick_base_minor))
  return fig


def draw_density_profile(fig, r, rho, err=None):
  ax = plt.subplot(2,1,2)
  im = ax.loglog(r, rho, linestyle='steps-mid-')
  line1 = ax.axvline(res_limit, color='black', linestyle=':')
  #ax.set_xlim(r_bins[0], r_bins[-1])
  ax.set_xlim(r[0] - (r[1]-r[0]), r[-1] + (r[-1]-r[-2]))
  ax.set_xlabel(xlabel_prof)
  ax.set_ylabel(ylabel_prof)
  if err != None:
    err_bars = ax.errorbar(r, rho, yerr=err,linestyle='None')
  return fig, ax


def draw_nfw_profile(fig, ax, r, rho, R_s=None):
  ax.loglog(r, rho, linestyle='-', color='red')
  if R_s != None:
    line = ax.axvline(R_s, color='purple', linestyle='-.')
  return fig


def calc_kinetic_energy(mass, vel):
  vsq = vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2
  energy = 0.5 * np.sum(mass*vsq)
  return energy


def calc_potential_energy(mass, pos):
  local_sqrt = np.sqrt
  partial_sum = 0.0
  for i in range(len(mass)):
    for j in range(len(mass)):
      if j != i:
        r_diff = local_sqrt((pos[i,0] - pos[j,0])**2 + (pos[i,1] - pos[j,1])**2 + (pos[i,2] - pos[j,2])**2)
        partial_sum = partial_sum - mass[i]*mass[j]/r_diff
  energy = partial_sum * grav_const / 2.0
  return energy


def calc_angular_momentum(mass, pos, vel):
  ang_mom_x = np.sum(mass * (pos[:,1] * vel[:,2] - pos[:,2] * vel[:,1]))
  ang_mom_y = np.sum(mass * (pos[:,2] * vel[:,0] - pos[:,0] * vel[:,2]))
  ang_mom_z = np.sum(mass * (pos[:,1] * vel[:,2] - pos[:,2] * vel[:,1]))
  ang_mom = np.sqrt(ang_mom_x**2 + ang_mom_y**2 + ang_mom_z**2)
  return ang_mom


def main():
  with open(outfile, 'w') as fd:
    #fd.write('#halo_mass  concentration  R_vir  R_s +- err  rho_0 +- err  alpha +- err  chi_square\n')
    fd.write('#halo_id        halo_mass           x_pos           y_pos           z_pos               c +-            err           R_vir             R_s +-            err           rho_0 +-            err      chi_square     nbins    N_part\n')
#  with open(comfile, 'w') as fd:
#    fd.write('#id mass dx dy dz\n')

#  if use_bgc2 == True:
#    header, halos, particles = bgc2.read_bgc2(sys.argv[1])
#    for i in range(len(halos)):
#      if halos[i][halo_id_col] == halo_id:
#        index = i
#    halo_particles = np.asarray(particles[index])
#    pos = halo_particles[:,pos_cols[0]:pos_cols[0]+3] * dist_scale
#    r_vir = halos[index][4] * dist_scale
#  else:
#    # Read in particle files
#    data = read_files(sys.argv[1:])
#    # Select particles with a given halo ID and convert positions from Mpc to kpc
#    if use_all == False:
#      halo_particles = data[np.where(data[:,halo_id_col] == halo_id)]
#    if use_all == True:
#      halo_particles = data
#    del data
#    pos = halo_particles[:,pos_cols[0]:pos_cols[0]+3] * dist_scale
#    r_vir = 241.48
#    #r_vir = pos.max()

  for input_file in sys.argv[1:]:
    if use_bgc2 == True:
      #header, halos, particles = bgc2.read_bgc2(sys.argv[1])
      header, halos, particles = bgc2.read_bgc2(input_file)
      halos = np.asarray(halos)
      indices = np.argsort(halos[:,2])      # sort by number of particles
      indices = indices[::-1]               # start with the biggest
    else:
      data = read_files([input_file])
      # Select particles with a given halo ID and convert positions from Mpc to kpc
      if use_all == False:
        particles = [data[np.where(data[:,halo_id_col] == halo_id)]]
      if use_all == True:
        particles = [data]
      del data


    itteration = 0
    #for index in range(len(halos)):
    #for index in range(1):
    #for index in indices[:10]:
    for index in indices:
      if ((len(particles[index]) >= 100) and (halos[index][1] == -1)):
        print '----------------------------------------------------------'

        halo_particles = np.asarray(particles[index])
        pos = halo_particles[:,pos_cols[0]:pos_cols[0]+3] * dist_scale
        vel = halo_particles[:,vel_cols[0]:vel_cols[0]+3]

        if use_bgc2 == True:
          halo_id = halos[index][0]
          r_vir = halos[index][4] * dist_scale
          halo_mass = halos[index][5]
          halo_pos = np.array([halos[index][6] * dist_scale, halos[index][7] * dist_scale, halos[index][8] * dist_scale])
          halo_vel = np.array([halos[index][9], halos[index][10], halos[index][11]])
        else:
          r_vir = 241.48
          halo_id = 0
          #halo_mass = mass[0] * len(halo_particles)
          halo_pos = np.array([0.0, 0.0, 0.0])
          halo_vel = np.array([0.0, 0.0, 0.0])


        if individual_masses == True:
          mass = halo_particles[:,mass_col] * mass_scale
        else:
          mass = np.ones(halo_particles.shape[0]) * common_mass * mass_scale

        if use_bgc2 == False:
          halo_mass = mass[0] * len(halo_particles)  #fix placement of this for ascii test

        print 'Using %d particles in halo %d.' % (halo_particles.shape[0], halo_id)

        # Find center of mass
        if find_com == True:
          mass_tot = mass.sum()
          m_pos = mass.reshape(mass.shape[0],1) * pos
          com = m_pos.sum(axis=0) / mass_tot
          pos = pos - com
          print 'Center of mass = (%g | %g | %g)' % (com[0], com[1], com[2])
        else:
          pos = pos - halo_pos
          vel = vel - halo_vel

        #with open(comfile, 'a') as fd:
        #  fd.write("%d %g %g %g %g\n" % (halo_id, halo_mass, halo_pos[0] - com[0], halo_pos[1] - com[1], halo_pos[2] - com[2]))

        # Bin halo particles into logarithmic shells and compute density
        r_bins, rho, rho_err = calc_density_profile(mass, pos)

        if len(r_bins) < 5:
          print 'Too few bins.  Skipping this halo.'
          with open(outfile, 'a') as fd:
            fd.write("%8d %16.12g  %14.10g  %14.10g  %14.10g  %14d +- %14d  %14d  %14d +- %14d  %14d +- %14d  %14d  %8d  %8d\n" % (halo_id, halo_mass, halo_pos[0], halo_pos[1], halo_pos[2], -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, len(halo_particles)))
          continue

      #  hist, r_bins = logbin(pos)
      #  err = poisson_error(hist)
      #  rho = mass * hist / (sphere_vol(r_bins[1:]) - sphere_vol(r_bins[:-1]))
      #  rho_err = err * rho
        mid_bins = 10.0**(0.5 * (np.log10(r_bins[1:]) + np.log10(r_bins[:-1])))
        print 'nbins = ', len(mid_bins)

        # Don't pass NaN's to fitting routine
        rho_err_nonan = np.copy(rho_err)
        nan_check = np.isnan(rho_err_nonan)
        for i in range(len(rho_err_nonan)):
          #if (nan_check[i] == True):
          #  rho[i] = 1.0e-10
          if (mid_bins[i] < res_limit) or (nan_check[i] == True):
            rho_err_nonan[i] = 1.0e10


#        r, m_enclosed = calc_m_enclosed(mass, pos)

        # Fit an NFW profile to the data
  #      try:
        nfw_r, nfw_fit, popt, pcov, chisq = fit_profile(mid_bins / r_vir, rho / rho.max(), err = rho_err_nonan / rho.max(), R_vir = 1.0)
        #nfw_r, nfw_fit, popt, pcov, chisq = fit_mass_profile(r / r_vir, m_enclosed / halo_mass)
        nfw_r = nfw_r * r_vir
        nfw_fit = nfw_fit * rho.max()
        #nfw_fit = nfw_fit * halo_mass
        scale_radius = popt[0] * r_vir
        scale_radius_err = pcov[0,0] * r_vir
        rho_0 = popt[1] * rho.max()
        rho_0_err = pcov[1,1] * rho.max()
        concentration = r_vir / scale_radius
        concentration_err = concentration * scale_radius_err / scale_radius

        # Print parameters
        print 'r_vir =', r_vir
        print "rho_0 = %g +/- %g" % (rho_0, rho_0_err)
        print "scale radius = %g +/- %g" % (scale_radius, scale_radius_err)
        print "concentration = %g +/- %g" % (concentration, concentration_err)

#put these back sometime##################################################################
#        kin_energy = calc_kinetic_energy(mass, vel)
#        pot_energy = calc_potential_energy(mass, pos)
#        ang_mom = calc_angular_momentum(mass, pos, vel)
#
#        ttow = 2.0 * abs(kin_energy / pot_energy)
#        lambda_spin = ang_mom * np.sqrt(abs(kin_energy + pot_energy)) / (grav_const * (np.sum(mass))**2.5)
        kin_energy = 0.0
        pot_energy = 0.0
        ang_mom = 0.0

        ttow = 0.0
        lambda_spin = 0.0
##########################################################################################


        if isinstance(pcov, float):
          print "inf covariance returned, skipping this halo..."
          with open(outfile, 'a') as fd:
            fd.write("%8d %16.12g  %14.10g  %14.10g  %14.10g  %14d +- %14d  %14d  %14d +- %14d  %14d +- %14d  %14d  %8d  %8d\n" % (halo_id, halo_mass, halo_pos[0], halo_pos[1], halo_pos[2], -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, len(halo_particles)))
          continue

        #Write parameters to file
        with open(outfile, 'a') as fd:
          #fd.write("%g  %g  %g  %g +- %g  %g +- %g  %g +- %g  %g\n" % (halo_mass, concentration, r_vir, scale_radius, pcov[0,0], rho_0, pcov[1,1], alpha, pcov[2,2], chisq))
          fd.write("%8d %16.12g  %14.10g  %14.10g  %14.10g  %14.10g +- %14.6g  %14.10g  %14.10g +- %14.6g  %14.10g +- %14.6g  %14.10g  %8d  %8d\n" % (halo_id, halo_mass, halo_pos[0], halo_pos[1], halo_pos[2], concentration, concentration_err, r_vir, scale_radius, scale_radius_err, rho_0, rho_0_err, chisq, len(r_bins), len(halo_particles)))
          

        ###################################################debug
        #blah_fit = nfw_fit_rho0(nfw_r, 20.0, 9.0e5)

        # Plot density profile histogram
        if (make_plot == True) and (itteration < 10):
          # Find the maximum of x, y, or z to be limit of projection plots
          plot_lim = pos.max()
          # Pick only a certain perentage of particles for projection plots
          if (draw_frac < 1.0):
            np.random.shuffle(pos)
            pos = pos[:(draw_frac*pos.shape[0])]

          fig = plt.figure()
          if draw_density == True:
            fig = draw_density_projection(fig, 0, plot_lim, pos[:,0], pos[:,1])
            fig = draw_density_projection(fig, 1, plot_lim, pos[:,0], pos[:,2])
            fig = draw_density_projection(fig, 2, plot_lim, pos[:,1], pos[:,2])
          else:
            fig = draw_projection(fig, 0, plot_lim, pos[:,0], pos[:,1])
            fig = draw_projection(fig, 1, plot_lim, pos[:,0], pos[:,2])
            fig = draw_projection(fig, 2, plot_lim, pos[:,1], pos[:,2])
          fig, ax = draw_density_profile(fig, mid_bins, rho, err=rho_err)  #put this back for binning
          #fig, ax = draw_density_profile(fig, r, m_enclosed)             #take this out for binning
          fig = draw_nfw_profile(fig, ax, nfw_r, nfw_fit, R_s=scale_radius)
          #fig = draw_nfw_profile(fig, ax, nfw_r, blah_fit, R_s=20.0)
          fig.tight_layout()
          plt.savefig(plot_base+str(itteration)+plot_ext)

          #sys.exit()

        itteration += 1


if __name__ == '__main__':
  main()

