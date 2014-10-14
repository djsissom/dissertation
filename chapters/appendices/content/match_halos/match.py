#!/usr/bin/env python

import sys
import getopt
import numpy as np


def main():
  # read and parse command line arguments
  opts, args = get_args(sys.argv[1:])
  output_file, match_file, densprof_files, parents_files, ascii_files = parse_args(opts, args)

  # read in headers as lists and data as numpy arrays
  match_header,     match_data     = read_files(match_file, header_line = 3)
  densprof_header1, densprof_data1 = read_files(densprof_files[0], header_line = 0)
  densprof_header2, densprof_data2 = read_files(densprof_files[1], header_line = 0)
  parents_header1,  parents_data1  = read_files(parents_files[0], header_line = 0)
  parents_header2,  parents_data2  = read_files(parents_files[1], header_line = 0)
  ascii_header1,    ascii_data1    = read_files(ascii_files[:(len(ascii_files)/2)], header_line = 0)
  ascii_header2,    ascii_data2    = read_files(ascii_files[(len(ascii_files)/2):], header_line = 0)
  print 'Finished reading files.'

  # filter matches, remove duplicate halo matches, and reorder match columns
  print 'Flitering match data...'
  match_data = filter_matches(match_data)
  if filter_duplicate_matches:
    match_data = filter_dups(match_data, unique_col = match_id1_col)
    match_data = filter_dups(match_data, unique_col = match_id2_col)
  if reorder_match_columns:
    match_header, match_data = reorder_match_cols(match_header, match_data)
  
  # calculate number of subhalos and add column to parents data and headers
  print 'Finding number of subhalos...'
  parents_header1.append('N_subs')
  parents_header2.append('N_subs')
  parents_data1 = count_subs(parents_data1)
  parents_data2 = count_subs(parents_data2)

  # create header
  print 'Making header...'
  header = make_header(match_header, densprof_header1, densprof_header2, \
                       parents_header1, parents_header2, ascii_header1, ascii_header2)

  # match halos
  print 'Matching halos...'
  halos = match_halos(match_data, [densprof_data1, densprof_data2, \
                      parents_data1, parents_data2, ascii_data1, ascii_data2])

  # filter based on given criteria and sort
  print 'Filtering halo data...'
  if filter_halo_properties:
    halos = filter_halos(halos)
  if sort_col != None:
    sort_mask = halos[:,sort_col].argsort()
    sort_mask = sort_mask[::-1]
    halos = halos[sort_mask]

  # output matched table
  print 'Writing resluts...'
  write_results(output_file, header, halos)

  print 'Finished.'


def get_args(arglist):
  try:
    opts, args = getopt.gnu_getopt(arglist, shortopts, longopts)
  except getopt.GetoptError:
    print "Invalid option(s)."
    print help_string
    sys.exit(2)
  if opts == []:
    print 'No options given.'
    print help_string
    sys.exit(2)
  return opts, args


def parse_args(opts, args):
  densproffiles = None
  parentsfiles = None
  asciifiles = None
  use_ascii = False
  for opt in opts:
    if (opt[0] == '-h') or (opt[0] == '--help') or (opts == None):
      print help_string
      sys.exit(0)
    if (opt[0] == '-o') or (opt[0] == '--outfile'):
      outfile = opt[1]
    if (opt[0] == '-m') or (opt[0] == '--match'):
      matchfile = opt[1]
    if (opt[0] == '-d') or (opt[0] == '--density'):
      densproffiles = create_append(densproffiles, opt[1])
    if (opt[0] == '-p') or (opt[0] == '--parents'):
      parentsfiles = create_append(parentsfiles, opt[1])
    if (opt[0] == '-a'):
      use_ascii = True
  if use_ascii:
    if len(args) % 2 != 0:
      print 'Must have an even number of ascii files!'
      sys.exit(3)
    for arg in args:
      asciifiles = create_append(asciifiles, arg)
  return outfile, matchfile, densproffiles, parentsfiles, asciifiles


def create_append(lst, value):
  if lst == None:
    lst = [value]
  else:
    lst.append(value)
  return lst


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
      data = np.genfromtxt(file, comments='#')
    else:
      data = np.append(data, np.genfromtxt(file, comments='#'), axis=0)

  if header_line == None:
    return data
  else:
    return header, data


def filter_matches(halos):
  if filter_bad_matches:
    halos = halos[halos[:,match_id1_col] != -1]
    halos = halos[halos[:,match_id2_col] != -1]
  if (min_npart != 0) and (min_npart != None):
    halos = halos[halos[:, match_npart1_col] >= min_npart]
    halos = halos[halos[:, match_npart2_col] >= min_npart]
  if (minperc_ncommon != 0) and (minperc_ncommon != None):
    halos = halos[halos[:, match_ncommon_col] / halos[:, match_npart1_col] >= minperc_ncommon]
    halos = halos[halos[:, match_ncommon_col] / halos[:, match_npart2_col] >= minperc_ncommon]
  return halos


def filter_dups(halos, unique_col = 0):
  ncommon   = halos[:, match_ncommon_col]
  n1        = halos[:, match_npart1_col]
  n2        = halos[:, match_npart2_col]
  rank      = ncommon**2 / (n1 * n2)  -  np.abs(n1 - n2) / (n1 + n2)

  sort_mask = np.argsort(rank)
  halos     = halos[sort_mask]
  
  unique, mask = np.unique(halos[:, unique_col], return_index=True)
  halos = halos[mask]
  return halos


def reorder_match_cols(match_header, match_data):
  global match_id1_col
  global match_id2_col
  global match_hnum1_col
  global match_hnum2_col
  global match_npart1_col
  global match_npart2_col
  global match_ncommon_col

  order = [match_id1_col, match_id2_col, \
           match_hnum1_col, match_hnum2_col, \
           match_npart1_col, match_npart2_col, \
           match_ncommon_col]
  match_header = [match_header[index] for index in order]
  match_data   = match_data[:, order]

  match_id1_col     = 0
  match_id2_col     = 1
  match_hnum1_col   = 2
  match_hnum2_col   = 3
  match_npart1_col  = 4
  match_npart2_col  = 5
  match_ncommon_col = 6

  return match_header, match_data


def count_subs(halos):
  id      = halos[:, id_col]
  parents = halos[:, parents_col]
  parents = parents[parents != -1]
  nsubs   = (id[:, np.newaxis] == parents).sum(axis = 1)
  halos   = np.column_stack((halos, nsubs))
  return halos


def make_header(match, densprof1, densprof2, parents1, parents2, ascii1, ascii2):
  # zeroeth line just lists column number
  total_len = len(match + densprof1 + densprof2 + parents1 + parents2 + ascii1 + ascii2)
  header_line0 = [str(i) for i in range(total_len)]
  header_line0 = '  '.join(header_line0)
  header_line0 = '#' + header_line0

  # first line denotes which file columns are from
  match_repeat = len(match) - 4
  densprof_repeat = len(densprof1 + densprof2) - 4
  parents_repeat = len(parents1 + parents2) - 4
  ascii_repeat = len(ascii1 + ascii2) - 4

  match_part    = '  '.join(['|---', 'cross', 'match'] + ['----'] * match_repeat + ['---|'])
  densprof_part = '  '.join(['|---', 'density', 'profile'] + ['----'] * densprof_repeat + ['---|'])
  parents_part  = '  '.join(['|---', 'rockstar', 'parents'] + ['----'] * parents_repeat + ['---|'])
  ascii_part    = '  '.join(['|---', 'rockstar', 'ascii'] + ['----'] * ascii_repeat + ['---|'])

  header_line1 = '  '.join([match_part, densprof_part, parents_part, ascii_part])
  header_line1 = '#' + header_line1

  # second line labels 2lpt and za columns
  tot_len = len(match + densprof1 + densprof2 + parents1 + parents2 + ascii1 + ascii2)
  header_line2 = ['2lpt' if i % 2 == 0 else 'za' if i % 2 == 1 else 'blah' for i in range(tot_len - 1)]
  header_line2.insert(len(match) - 1, 'matched')
  header_line2 = '  '.join(header_line2)
  header_line2 = '#' + header_line2

  # third line pulls labels from original file headers
  match_part    = match
  densprof_part = interweave(densprof1, densprof2)
  parents_part  = interweave(parents1, parents2)
  ascii_part    = interweave(ascii1, ascii2)

  header_line3 = match_part + densprof_part + parents_part + ascii_part
  header_line3 = '  '.join(header_line3)
  header_line3 = '#' + header_line3

  header = [header_line0, header_line1, header_line2, header_line3]
  return header


def interweave(list1, list2):
  newlist = list1 + list2
  newlist[::2]  = list1
  newlist[1::2] = list2
  return newlist


def interweave_np_2d(array1, array2):
  newarray = np.empty((len(array1), len(array1[0]) + len(array2[0])))
  newarray[:,::2]  = array1
  newarray[:,1::2] = array2
  return newarray


def match_halos(matches, arrays):
  halos = matches.copy()
  for i, array in enumerate(arrays):
    if array != None:
      match_id_col = i % 2
      halos = sort_stack(halos, array, match_id_col)

  # interweave columns so that matching 2lpt/za columns are adjacent
  tmp_halos = halos
  halos = np.empty((len(tmp_halos), len(tmp_halos[0])))
  halos[:,:len(matches[0])] = matches
  startcol = len(matches[0])
  for i in range(0, len(arrays), 2):
    colrange1 = len(arrays[i][0])
    colrange2 = len(arrays[i+1][0])
    endcol = startcol + colrange1 + colrange2

    cols1 = tmp_halos[:,startcol:startcol+colrange1]
    cols2 = tmp_halos[:,startcol+colrange1:startcol+colrange1+colrange2]

    halos[:,startcol:endcol] = interweave_np_2d(cols1, cols2)
    startcol = endcol
  return halos


def sort_stack(halos, array, match_id_col):
  # add empty columns to halos to later fill with halo data
  rows = len(halos)
  origcols = len(halos[0])
  newcols = len(array[0])
  empty = np.empty((rows, newcols))
  empty[:] = np.nan
  halos = np.column_stack((halos, empty))

  # remove halos from array with no matches
  match_id = halos[:, match_id_col]
  array_id = array[:, id_col]
  array_mask = np.in1d(array_id, match_id)
  array = array[array_mask]

  # create mask so we only add lines for halos in array
  array_id = array[:, id_col]
  halo_mask = np.in1d(match_id, array_id)
  masked_halos = halos[halo_mask]

  # create masks to sort by halo id
  match_id_sort_mask = np.argsort(masked_halos[:, match_id_col])
  sorted_masked_halos = masked_halos[match_id_sort_mask]

  # sort array by halo id and copy to empty columns of view of halos
  array_id_sort_mask = np.argsort(array[:,id_col])
  sorted_masked_halos[:, origcols:] = array[array_id_sort_mask]

  # 'unmask' - put data back in original halos
  masked_halos[match_id_sort_mask] = sorted_masked_halos
  halos[halo_mask] = masked_halos

  return halos


def filter_halos(halos):
  #todo
  return halos


def write_results(output_file, header, halos):
  format = get_format(halos[0])
  with open(output_file, 'w') as fd:
    for line in header:
      fd.write(line + '\n')
    np.savetxt(fd, halos, fmt=format)


def get_format(line):
  format = ['%d' if col in int_cols else '%1.14g' for col in range(len(line))]
  format = ' '.join(format)
  return format


help_string = '''
Available options are:
    -h, --help
    -v, --verbose
    -o <outfile>, --outfile <outfile>
    -m <matchlist>, --match <matchlist>
    -d <densityprofile_file>, --density <densityprofile_file>
    -p <parents_file>, --parents <parents_file>
    -a <ascii_files>, --ascii <ascii_files> - must be last)
'''
shortopts = "hvo:m:d:p:a"
longopts  = ["help", "verbose", "outfile=", "matchfile=", "density=", "parents=", "ascii"]

lt_cols = []
lt_vals = []

gt_cols = []
gt_vals = []

eq_cols = []
eq_vals = []

ne_cols = []
ne_vals = []

#int_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8]
int_cols = []

match_id2_col             =  1
match_npart2_col          =  2
match_id1_col             =  4
match_npart1_col          =  5
match_ncommon_col         =  6
match_hnum1_col           =  3
match_hnum2_col           =  0

id_col                    =  0      #  col of each input file
sort_col                  = 47      #  col of final table - use None to turn off sorting
parents_col               = -1

filter_bad_matches        =  True
filter_duplicate_matches  =  False
reorder_match_columns     =  True
filter_halo_properties    =  False
min_npart                 =  20     #  use 0 or None to use all size halos
minperc_ncommon           =  0.05   #  a fraction, use 0 or None to use any match percent


if __name__ == '__main__':
  main()
