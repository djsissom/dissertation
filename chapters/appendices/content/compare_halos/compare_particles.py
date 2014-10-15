#!/usr/bin/env python

import sys
import bgc2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from scipy.stats import chisquare

#id1, id2 = 727, 4420 # 2lpt first
#id1, id2 = 4416, 727 # za first

#id1, id2 = 4416, 4420 # both za
#id1, id2 = 4416, 4416 # both za

#id1, id2 = 653, 4355
#id1, id2 = 38, 3803
#id1, id2 = 155099, 80362
#id1, id2 = 98722, 14357
id1, id2 = 84289, 143514


#read_mode = 'ascii2'
read_mode = 'bgc2'

if read_mode == 'bgc2':
  use_bgc2 = True
  use_all = False
  multiple_halos = True
  individual_masses = False
  halo_id = 146289
  nbins = 50
  nfit = 500
  ooms = 3.0
  mass_scale = 1.0
  common_mass = 5.33423e5
  dist_scale = 1.0e3
  #res_limit = 0.488
  res_limit = 4.0
  #res_limit = 10.0
  #draw_frac = 1.0e-2
  draw_frac = 1.0
  tick_base_major = 10.0
  tick_base_minor = 1.0
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
  res_limit = 3.0e-1
  draw_frac = 1.0e-2
  tick_base_major = 200.0
  tick_base_minor = 40.0
else:
  sys.exit(98712)

outfile = 'halo_properties.txt'
comfile = 'center_of_mass.txt'

make_plot = True
plot_base = 'density_profile.fig.'
plot_ext = '.eps'
dist_units = 'kpc'
xlabel_proj = [r'X Position (%s h$^{-1}$)' % (dist_units), r'X Position (%s h$^{-1}$)' % (dist_units), r'Y Position (%s h$^{-1}$)' % (dist_units)]
ylabel_proj = [r'Y Position (%s h$^{-1}$)' % (dist_units), r'Z Position (%s h$^{-1}$)' % (dist_units), r'Z Position (%s h$^{-1}$)' % (dist_units)]
xlabel_prof = r'Radius (%s h$^{-1}$)' % (dist_units)
ylabel_prof = r'Density (M$_{\odot}$ %s$^{-3}$ h$^{2}$)' % (dist_units)

#common_mass = 1.0e-7
#common_mass = 1.0e5
mass_col = 0
pos_cols = (1,2,3)
vel_cols = (4,5,6)
halo_id_col = 0

grav_const = 4.3e-6 # kpc M_sol^-1 (km/s)^2

profile_type = 0  # 0 -> nfw, fit rho_0
                  # 1 -> nfw, calculate rho_0
                  # 2 -> nfw, rho_0 middle of leftmost bin above resolution
                  # 3 -> fit outer slope, fit rho_0
                  # 4 -> fit outer slope, calculate rho_0
                  # 5 -> fit outer slope, rho_0 middle of leftmost bin above resolution

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


def fit_profile(r, rho, err=None, R_vir=None):
  for i in range(len(r)):
    if r[i] > res_limit:
      rho_0_databin = rho[i]
      first_good_bin = i
      break
  #--------- choose one fitting type ---------#
  #popt, pcov = curve_fit(nfw_fit_rho0, r, rho, sigma=err)
  #popt, pcov = curve_fit(nfw_def_rho0(R_vir), r, rho, p0=[10.0], sigma=err)
  #popt, pcov = curve_fit(nfw_databin_rho0(rho_0_databin), r, rho, sigma=err)
  blah = 2
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
    popt, pcov = curve_fit(nfw_fit_rho0, r, rho, sigma=err)
      
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

  chi2 = chisquare(np.log10(rho[first_good_bin:]), np.log10(chi2_fit[first_good_bin:]))
  chi2_nolog = chisquare(rho[first_good_bin:], chi2_fit[first_good_bin:])
  print 'chi_square =', chi2
  print 'chi_square_nolog =', chi2_nolog
  return nfw_r, nfw_fit, popt, pcov, chi2[0]


def draw_projection(fig, place, plot_lim, x, y):
  ax = plt.subplot(1,3,place+1, aspect='equal')
  im = ax.plot(x, y, linestyle='', marker='.', markersize=1, markeredgecolor='blue')
  ax.set_xlabel(xlabel_proj[place])
  ax.set_ylabel(ylabel_proj[place])
  ax.set_xlim(-plot_lim, plot_lim)
  ax.set_ylim(-plot_lim, plot_lim)
#  ax.xaxis.set_major_locator(MultipleLocator(tick_base_major))
#  ax.xaxis.set_minor_locator(MultipleLocator(tick_base_minor))
#  ax.yaxis.set_major_locator(MultipleLocator(tick_base_major))
#  ax.yaxis.set_minor_locator(MultipleLocator(tick_base_minor))
  return fig, ax


def draw_projection_again(fig, ax, x, y):
  im = ax.plot(x, y, linestyle='', marker='.', markersize=1, markeredgecolor='red')
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
  #for input_file in sys.argv[1:]:
  #header1, halos1, particles1 = bgc2.read_bgc2(sys.argv[1])
  #header2, halos2, particles2 = bgc2.read_bgc2(sys.argv[2])

  nargs = len(sys.argv) - 1
  if (float(nargs) % 2.0) != 0.0:
    print 'number of arguments must be even'
    sys.exit()

  for i in range(nargs / 2):
    i += 1
    temp_header1, temp_halos1, temp_particles1 = bgc2.read_bgc2(sys.argv[i])
    temp_header2, temp_halos2, temp_particles2 = bgc2.read_bgc2(sys.argv[(nargs / 2) + i])
    if i == 1:
      halos1, particles1 = temp_halos1, temp_particles1
      halos2, particles2 = temp_halos2, temp_particles2
    else:
      halos1 = np.append(halos1, temp_halos1, axis=0)
      halos2 = np.append(halos2, temp_halos2, axis=0)
      particles1 = np.append(particles1, temp_particles1, axis=0)
      particles2 = np.append(particles2, temp_particles2, axis=0)

  halos1 = np.asarray(halos1)
  halos2 = np.asarray(halos2)
  #indices = np.argsort(halos[:,2])      # sort by number of particles
  #indices = indices[::-1]               # start with the biggest

  itteration = 0
  #for index in indices[:1000]:
  #for index in indices:
  for index in range(halos1.shape[0]):
    halo_id = halos1[index,0]
    if (halo_id == id1):
      print '----------------------------------------------------------'

      halo_particles1 = np.asarray(particles1[index])
      pos1 = halo_particles1[:,pos_cols[0]:pos_cols[0]+3] * dist_scale
      #vel1 = halo_particles1[:,vel_cols[0]:vel_cols[0]+3]

      r_vir1 = halos1[index][4] * dist_scale
      halo_mass1 = halos1[index][5]
      halo_pos1 = np.array([halos1[index][6] * dist_scale, halos1[index][7] * dist_scale, halos1[index][8] * dist_scale])
      #halo_vel1 = np.array([halos1[index][9], halos1[index][10], halos1[index][11]])

      print 'Using %d particles in halo %d.' % (halo_particles1.shape[0], halo_id)

      # Find center of mass
      #pos = pos - halo_pos
      #vel = vel - halo_vel

      # Pick only a certain perentage of particles for projection plots
      if (draw_frac < 1.0):
        np.random.shuffle(pos1)
        pos1 = pos1[:(draw_frac*pos1.shape[0])]

  for index in range(halos2.shape[0]):
    halo_id = halos2[index,0]
    if (halo_id == id2):
      print '----------------------------------------------------------'

      halo_particles2 = np.asarray(particles2[index])
      pos2 = halo_particles2[:,pos_cols[0]:pos_cols[0]+3] * dist_scale
      #vel2 = halo_particles2[:,vel_cols[0]:vel_cols[0]+3]

      r_vir2 = halos2[index][4] * dist_scale
      halo_mass2 = halos2[index][5]
      halo_pos2 = np.array([halos2[index][6] * dist_scale, halos2[index][7] * dist_scale, halos2[index][8] * dist_scale])
      #halo_vel2 = np.array([halos2[index][9], halos2[index][10], halos2[index][11]])

      print 'Using %d particles in halo %d.' % (halo_particles2.shape[0], halo_id)

      # Find center of mass
      #pos = pos - halo_pos
      #vel = vel - halo_vel

      # Pick only a certain perentage of particles for projection plots
      if (draw_frac < 1.0):
        np.random.shuffle(pos2)
        pos2 = pos2[:(draw_frac*pos2.shape[0])]

  # Find the maximum of x, y, or z to be limit of projection plots
  center = (halo_pos1 + halo_pos2) / 2.0
  pos1 = pos1 - center
  pos2 = pos2 - center
  halo_pos1 = halo_pos1 - center
  halo_pos2 = halo_pos2 - center
  plot_lim = np.append(pos1, pos2).max()

  # Plot density profile histogram
  if (make_plot == True):
    fig = plt.figure()

    fig, ax = draw_projection(fig, 0, plot_lim, pos1[:,0], pos1[:,1])
    fig = draw_projection_again(fig, ax, pos2[:,0], pos2[:,1])
    ax.add_patch(Circle((halo_pos1[0], halo_pos1[1]), r_vir1, fc="None", ec="black", lw=1))
    ax.add_patch(Circle((halo_pos2[0], halo_pos2[1]), r_vir2, fc="None", ec="black", lw=1))

    fig, ax = draw_projection(fig, 1, plot_lim, pos1[:,0], pos1[:,2])
    fig = draw_projection_again(fig, ax, pos2[:,0], pos2[:,2])
    ax.add_patch(Circle((halo_pos1[0], halo_pos1[2]), r_vir1, fc="None", ec="black", lw=1))
    ax.add_patch(Circle((halo_pos2[0], halo_pos2[2]), r_vir2, fc="None", ec="black", lw=1))

    fig, ax = draw_projection(fig, 2, plot_lim, pos1[:,1], pos1[:,2])
    fig = draw_projection_again(fig, ax, pos2[:,1], pos2[:,2])
    ax.add_patch(Circle((halo_pos1[1], halo_pos1[2]), r_vir1, fc="None", ec="black", lw=1))
    ax.add_patch(Circle((halo_pos2[1], halo_pos2[2]), r_vir2, fc="None", ec="black", lw=1))

    #fig, ax = draw_density_profile(fig, mid_bins, rho, err=rho_err)
    #fig = draw_nfw_profile(fig, ax, nfw_r, nfw_fit, R_s=scale_radius)
    fig.tight_layout()
    #plt.savefig(plot_base+str(itteration)+plot_ext)
    plt.savefig('test.eps')


if __name__ == '__main__':
  main()

