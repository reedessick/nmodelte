#!python_alias
usage="""written to convert between variables and take ffts. "c" is for "conver" """

import nmode_utils as nm_u

from optparse import *

####################################################################################################
#
#
#                        Parse input options
#
#
####################################################################################################

parser = OptionParser(usage=usage)

parser.add_option("-i", "--input-filename", default=False, type="string", help="the input filename.")
parser.add_option("", "--min", default=False, type="float", help="minimum value for either time or frequency")
parser.add_option("", "--max", default=False, type="float", help="maximum value for either time or frequency")
parser.add_option("", "--downsample", default=False, type="float", help="only read in one out of this many points from --input-filename")

parser.add_option("", "--convert-outfilename", default=False, type="string", help="the output filename for conversion between x <--> q")
parser.add_option("", "--fft-outfilename", default=False, type="string", help="the output filename for forward fft")
parser.add_option("", "--amp-fft-outfilename", default=False, type="string", help="the output filename for forward fft of amplitude |x|,|q|")

parser.add_option("-c", "--current", default=False, type="string", help="the current variable type")
parser.add_option("-t", "--target", default=False, type="string", help="the target variable type")
parser.add_option("-C", "--convert", default=False, action="store_true", help="convert current -> target in the time domain")
parser.add_option("-l", "--logfilename", default=False, type="string", help="log file to use with nmode_utils.convert()")

parser.add_option("-f", "--fft", default=False, action="store_true", help="perform a forward fft: a --> \\tilde{a}")

parser.add_option("", "--amp-fft", default=False, action="store_true", help="perform a forward fft: |a| --> \\tilde{|a|}")

parser.add_option("-v", "--verbose", default=False, action="store_true")

opts, args = parser.parse_args()

if not opts.input_filename:
  opts.input_filename = raw_input("input_filename = ")

if opts.convert and (not opts.current):
  opts.current = raw_input("current = ")

if opts.convert and (not opts.target):
  opts.target = raw_input("target = ")

if opts.convert and (not opts.logfilename):
  opts.logfilename = raw_input("logfilename = ")

####################################################################################################
#
#
#                       load data
#
#
####################################################################################################
### load data
if opts.verbose: print "loading data from %s" % opts.input_filename
t_P, q, N_m = nm_u.load_out(opts.input_filename, tmin=opts.min, tmax=opts.max, downsample=opts.downsample)

if opts.convert:
  if opts.verbose: print "loading network from %s" % opts.logfilename
  system = nm_u.load_log(opts.logfilename)
  network = system.network

####################################################################################################
#
#
#              convert current -> target
#
#
####################################################################################################
if opts.convert:
  if opts.verbose: print "converting variable types in time domain : %s -> %s"%(opts.current, opts.target)
  converted_t_P, converted_q, converted_current = nm_u.convert(t_P, q, opts.current, opts.target, system, system.Porb)

  if opts.convert_outfilename:
    if opts.verbose: print "writing converted data to %s" % opts.convert_outfilename
    outfile = open(opts.convert_outfilename, "w")
    for ind in range(len(converted_t_P)):
      print >> outfile, nm_u.report_func(converted_t_P[ind], [ l for Q in converted_q for l in Q[ind] ], P=1., phase=False)
    outfile.close()

  else:
    for ind in range(len(t_P)):
      print nm_u.report_func(t_P, [ l for Q in q for l in Q[ind] ], P=1., phase=False)

####################################################################################################
#
#
#            perform forward FFT
#
#
####################################################################################################

if opts.fft:
  if opts.verbose: print "computing forward FFT"
  freq, fq, N_m = nm_u.nmode_fft(t_P, q, verbose=opts.verbose)

  if opts.fft_outfilename:
    if opts.verbose: print "writing FFT data to %s" % opts.fft_outfilename
    outfile = open(opts.fft_outfilename, "w")
    for ind in range(len(freq)):
      print >>outfile, nm_u.report_func(freq[ind], [l for Q in fq for l in Q[ind]], P=1.)
    outfile.close()
  else:
    for ind in range(len(freq)):
      print nm_u.report_func(freq[ind], [l for Q in fq for l in Q[ind]], P=1.)
    
####################################################################################################
#
#
#       perform forward FFT of |a| (amplitudue abs value)
#
#
####################################################################################################

if opts.amp_fft:
  if opts.verbose: print "computing forward FFT of abs(a)"
  freq, fA, N_m = nm_u.nmode_fft(t_P, [ [ [nm_u.amp(l1,l2), 0] for l1, l2 in Q] for Q in q], verbose=opts.verbose)

  if opts.amp_fft_outfilename:
    if opts.verbose: print "writing abs(a) FFT data to %s" % opts.amp_fft_outfilename 
    outfile = open(opts.amp_fft_outfilename, "w")
    for ind in range(len(freq)):
      print >>outfile, nm_u.report_func(freq[ind], [l for Q in fA for l in Q[ind]], P=1.)
    outfile.close()
  else:
    for ind in range(len(freq)):
      print nm_u.report_func(freq[ind], [l for Q in fA for l in Q[ind]], P=1.)




