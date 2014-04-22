#!python_alias
usage = """ written to mitgate the disk footprint of catalog files by re-writing them with only the top X couplings retained. As long as X is a large number (say, X > 10,000) then we will never have to worry about running off the end of the list. this is NOT intended to be used with "minima" lists"""

import os
import sys
import glob

import mode_selection
from mode_selection import coupling_list
from ggg import ggg_coupling_list

####################################################################################################
#
#
#  main function
#
#
####################################################################################################
def clean(old_ctg_filename, new_ctg_filename, max_num_pairs, alpha, c, wo, k_hat, coupling_type="general", verbose=False):
  """
  reads in old_ctg_filename and writes a (shorter) list into new_ctg_filename
  """
  if verbose: print "reading in %d pairs from : %s" % (max_num_pairs, old_ctg_filename)

  # load using mode_selection.coupling_list()
  if coupling_type == "general":
    this_list = coupling_list().load_unsorted_mode_lists([old_ctg_filename], num_pairs=opts.max_num_pairs, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False)

  # load using ggg.ggg_coupling_list()
  elif coupling_type == "ggg":
    if "heuristic" in old_ctg_filename:
      metric = "heuristic"
    elif "Ethr" in old_ctg_filename:
      metric = "Ethr"
    else:
      metric = raw_input("could not determine \"metric\" from filename : %s\n\tplease enter metric [heuristic/Ethr]" % old_ctg_filename)

    this_list = ggg_coupling_list(alpha, c, wo, k_hat).load_unsorted_mode_lists(metric, [old_ctg_filename], num_pairs=max_num_pairs, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False)

  # write
  if verbose: print "\twriting %d pairs to : %s" % (max_num_pairs, new_ctg_filename)
  this_list.write_mode_list(new_ctg_filename)

  return True

####################################################################################################
#
#
#                                  MAIN
#
#
####################################################################################################
if __name__ == "__main__":

  from optparse import OptionParser

  parser = OptionParser(usage=usage)

  parser.add_option("-v", "--verbose", action="store_true")

  parser.add_option("", "--max-num-pairs", default=10000, type='int', help='the maxmimum number of couplings to be kept in any single *ctg file')

  parser.add_option("", "--catalog-dir", default=False, type="string", help="the name of a directory in which *ctg files are stored. If provided, the script will clean up all *ctg files in 'catalog_dir'")
  parser.add_option("", "--catalog-cachefilename", default=False, type="string", help="the name of a cachefile containing the names of *ctg files to be cleaned")

  parser.add_option("", "--output-dir", default=False, type="string", help="the output directory into which new *ctg files will be written. If NOT provided, new *ctg files will be written over the existing files.")

  parser.add_option("", "--coupling-type", default="general", type="string", help="the type of couplings stored in these *ctg files, used to determine how to load catalogs. If NOT provided, script attempts to read them using standard mode_selection.coupling_list() objects. eg: \"ggg\"")

  parser.add_option("", "--alpha", default=False, type="float")
  parser.add_option("", "--c", default=False, type="float")
  parser.add_option("", "--wo", default=False, type="float")
  parser.add_option("", "--k-hat", default=False, type="float")

  opts, args = parser.parse_args()

  ### ensure we want to over-write existing files!
  if not opts.output_dir:
    ans = raw_input("WARNING! --output-dir NOT supplied and new *ctg files will replace existing files. Are you sure you want to continue? [yes/NO] ")
    while ans not in ["yes", "NO"]:
      ans = raw_input("please enter either \"yes\" or \"NO\" ")

    if ans == "NO":
      sys.exit()
  elif not os.path.exists(opts.output_dir):
    os.makedirs(opts.output_dir)

  ### define the type of coupling list we wish to use
  if opts.coupling_type not in ["general", "ggg"]:
    raise ValueError("--coupling-type %s NOT understood" % opts.coupling_type)

  if opts.coupling_type == "ggg":
    if not opts.alpha:
      opts.alpha = float(raw_input("alpha = "))
    if not opts.c:
      opts.c = float(raw_input("c = "))
    if not opts.wo:
      opts.wo = float(raw_input("wo = "))
    if not opts.k_hat:
      opts.k_hat = float(raw_input("k_hat = "))

  ##################################################################################################
  #
  #
  #                           collect all old catalog filenames
  #
  #
  ##################################################################################################

  old_ctg_filenames = []

  if opts.catalog_dir:
    if opts.verbose: print "globbing all *ctg filenames in : "+opts.catalog_dir
    old_ctg_filenames += sorted(glob.glob(opts.catalog_dir+"/*ctg"))

  if opts.catalog_cachefilename:
    if opts.verbose: print "reading in *ctg filenames from : "+opts.catalog_cachefilename
    f = open(opts.catalog_cachefilename, "r")
    for line in f:
      old_ctg_filenames.append( line.strip() )
    f.close()

  ##################################################################################################
  #
  #
  #                            iterate and write new catalog files
  #
  #
  ##################################################################################################

  if opts.verbose: 
    counter = 1
    N = len(old_ctg_filenames)

  for old_ctg_filename in old_ctg_filenames:
    # define new filename
    if opts.output_dir:
      new_ctg_filename = opts.output_dir + "/" + old_ctg_filename.split("/")[-1]
    else:
      new_ctg_filename = old_ctg_filename

    if opts.verbose: print "file %d / %d" % (counter, N)

    # delegate actual cleaning
    clean(old_ctg_filename, new_ctg_filename, opts.max_num_pairs, opts.alpha, opts.c, opts.wo, opts.k_hat, coupling_type=opts.coupling_type, verbose=opts.verbose)

    if opts.verbose: counter += 1



