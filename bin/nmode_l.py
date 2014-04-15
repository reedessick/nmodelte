#!/usr/local/bin/python

usage="""
python nmode_l.py [--options] outfilename1 outfilename2

written to compare two sets of integration output. This module computes residuals in the time-domain and manipulates them.

It checks for basic network consistency (finds which modeNos correspond to each other between the networks, consistency of time-sampling, etc).

Assumes both outfiles represent the same integration variable (x vs q).

If less than 2 logfilenames are specified through the options, assumes networks are identical (can only check number of modes).

If plots are requested, it generates them in the specified output-dir.

If numerical values are requested, it prints them to the terminal ( > into a file if you want to store them)
"""

import sys
import nmode_plotting as nmp
import nmode_utils as nmu
import liapunov

from optparse import *

##################################################
parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")
parser.add_option("-t", "--tag", default=False, type="string", help="a tag that is inserted into figure names")

### options for loadin data

parser.add_option("", "--logfilename1", default=False, type="string", help="logfilename corresponding to the first outfilename supplied")
parser.add_option("", "--logfilename2", default=False, type="string", help="logfilename corresponding to the second outfilename supplied")

parser.add_option("", "--downsample1", default=False, type="int")
parser.add_option("", "--downsample2", default=False, type="int")

parser.add_option("", "--tmin", default=False, type="float")
parser.add_option("", "--tmax", default=False, type="float")

parser.add_option("", "--time-offset1", default=0.0, type="float", help="a time offset applied to the integration data from the first outfilename supplied. t1 -> t1 + time_offset1*Porb1")
parser.add_option("", "--time-offset2", default=0.0, type="float", help="a time offset applied to the integration data from the second outfilename supplied. t2 -> t2 + time_offset2*Porb2")

### options for computing statistics

parser.add_option("", "--liapunov", default=False, action="store_true", help="estimate the liapunov exponent from the observed residuals")
parser.add_option("", "--sum-sqr-err", default=False, action="store_true", help="compute the sum-square-errors from the residuals. This is normalized by the number of samples")
parser.add_option("", "--chi2", default=False, action="store_true", help="compute a chi2 goodness of fit test, normalized by the number of samples")
parser.add_option("", "--inner-product", default=False, action="store_true", help="compute the inner-product normalized by the amplitude of each time series.")

### options for plotting
parser.add_option("", "--res", default=False, action="store_true", help="plot residuals")
#parser.add_option("", "--amp-res", default=False, action="store_true", help="plot amplitudes and residuals")
parser.add_option("", "--cross-correlation", default=False, action="store_true", help="compute the cross correlation between the time series, plotted as both time and frequency series")

parser.add_option("-g", "--grid", default=False, action="store_true")
parser.add_option("", "--no-legend", default=False, action="store_true")

parser.add_option("", "--output-dir", default="./", type="string")

opts, args = parser.parse_args()

if len(args) != 2:
  sys.exit("must supply exactly 2 outfilenames as arguments")  

### flags for control sequences later on
stats = (opts.liapunov or opts.sum_sqr_err or opts.chi2 or opts.inner_product or opts.cross_correlation)
plots = (opts.res or opts.amp_res)

##################################################
#
#   load data 
#
##################################################

filename1, filename2 = args

if opts.verbose: print "loading data from", filename1
t_P1, q1, N_m1 = nmu.load_out(filename1, tmin=opts.tmin, tmax=opts.tmax, downsample=opts.downsample1)

if opts.logfilename1:
  if opts.verbose: print "loading system information from", opts.logfilename1
  system1 = nmu.load_log(opts.logfilename1)


if opts.verbose: print "loading data from ", filename2
t_P2, q2, N_m2 = nmu.load_out(filename2, tmin=opts.tmin, tmax=opts.tmax, downsample=opts.downsample2)

if opts.logfilename2:
  if opts.verbose: print "loading system information from", opts.logfilename2
  system2 = nmu.load_log(opts.logfilename2)

##################################################
#
#  check for network consistency/compatibility
#
##################################################

### check modes included in both networks ==> common_modes
if opts.logfilename1 and opts.logfilename2:

  common_nlms = nmu.common_nlms([m.get_nlms() for m in system1.network.modes], [m.get_nlms() for m in system2.network.modes])
  if len(common_nlms) == 0:
    sys.exit("no common modes between these networks")

  t1 = (np.array(t_P1) + opts.time_offset1)*system1.Porb
  t2 = (np.array(t_P2) + opts.time_offset2)*system2.Porb

else:

  if N_m1 != N_m2:
    sys.exit("networks are not the same size (network1: %d modes, network2: %d modes) and too few logfilenames were supplied to compute conversion between networks!" % (N_m1, N_m2))
  elif opts.logfilename1:
    common_nlms = [ (m.get_nlms(), i, i) for i, m in enumerate(system1.network.modes) ] # assumes both networks are identical
    t1 = (np.array(t_P1) + opts.time_offset1)*system1.Porb
    t2 = (np.array(t_P2) + opts.time_offset2)*system1.Porb
  elif opts.logfilename2:
    common_nlms = [ (m.get_nlms(), i, i) for i, m in enumerate(system2.network.modes) ] # assumes both networks are identical
    t1 = (np.array(t_P1) + opts.time_offset1)*system2.Porb
    t2 = (np.array(t_P2) + opts.time_offset2)*system2.Porb
  else:
    common_nlms = [ (False, i, i) for i in range(N_m1) ] # assumes both networks are identical
    t1 = np.array(t_P1) + opts.time_offset1
    t2 = np.array(t_P2) + opts.time_offset2
  
### check for overlapping time ranges
if max(t1[0], t2[0]) > min(t1[-1], t2[-1]):
  sys.exit("no time range overlap between data\n\tt1 : [%f, %f]\n\tt2 : [%f, %f]" % (t1[0], t1[-1], t2[0], t2[-1]) )

##################################################
#
#  compute residuals between integration data
#
##################################################

if opts.verbose: print "computing residuals"
rt, r, rq1, rq2 = liapunov.compute_residuals(t1, t2, q1, q2, common_nlms) # ain't delegation a bitch

##################################################
#
#  generate summary information about residuals
#
##################################################

if stats:

  if opts.verbose: print "generating statistics"

  print "filename1 = ", filename1
  print "Dt_1/Porb1 = %f" % opts.time_offset1

  print "filename1 = ", filename2
  print "Dt_2/Porb2 = %f" % opts.time_offset2

  print "\n%f < t < %f" % rt[0], rt[-1]
  print "dt = %f" % rt[1]-rt[0]
  print "num samples = %d" % len(rt)

##############################################################################################
#
# compute and report individual statistics
#
##############################################################################################

  print "Network Statistics:"

  if opts.liapunov: # compute the numerical value of the liapunov exponent: essentially fit residuals to an exponential
    (_, Ty_fit, Tchi2_red), (_, y_fit, chi2_red) = liapunov.nmode_liapunov_exp(rt, r)
    #     network totals             for each mode
    print "\ty_fit = %fe%d" % nmu.float_to_scientific(Ty_fit)
    print "\tchi2_red = %f" % Tchi2_red

  if opts.sum_sqr_err: # compute the sum-square-errors statistic  
    tot_SSE, SSE = liapunov.sum_sqr_err(r)
    # these are normalized by the number of sample points
    print "\tsum{r**2}/len(r) = %fe%d" % nmu.float_to_scientific(tot_SSE)

  if opts.chi2: # compute the chi2 goodness of fit statistic
    Tchi2s, chi2s = liapunov.nmode_chi2(rq1, rq2)
    # returns the reduced chi2 goodness of fit test between the two time series
    print "\tchi2{rq1,rq2}/len(rq1) = %fe%d" % nmu.float_to_scientific(Tchi2s)

  if opts.inner_product: # compute the inner product of the data, normalized by the inner products of the separate time streams
    Tinner_prod, inner_prods = liapunov.nmode_inner_product(rq1, rq2)
    #  standard definition of inner product in time domain
    print "\t(rq1|rq2) = %fe%d" % nmu.float_to_scientific(Tinner_prod)

  # for individual modes
  for ind, (nlms, modeNo1, modeNo2) in enumerate(common_nlms):
    print "modeNo1 = %d\nmodeNo2 = %d" % (modeNo1, modeNo2)
    if nlm:
      print "\tn, l, m, s = %d, %d, %d, %d" % nlms

    if opts.liapunov: # compute the numerical value of the liapunov exponent: essentially fit residuals to an exponential
      print "\t(y_fit)_i = %fe%d" % nmu.float_to_scientific(y_fit[ind])
      print "\t(chi2_red)_i = %f" % chi2_red[ind]

    if opts.sum_sqr_err: # compute the sum-square-errors statistic  
      print "\tsum{(r_i)**2}/len(r_i) = %fe%d" % nmu.float_to_scientific(SSE[ind])

    if opts.chi2: # compute the chi2 goodness of fit statistic
      print "\tchi2{rq1_i, rq2_i}/len(rq1_i) = %fe%d" % nmu.float_to_scientific(chi2s[ind])

    if opts.inner_product: # compute the inner product of the data, normalized by the inner products of the separate time streams
      print "\t(rq1_i|rq2_i) = %fe%d" % nmu.float_to_scientific(inner_prods[ind])

##################################################
#
#  plotting
#
##################################################

if plots:

  if opts.tag:
    opts.tag += "_"

  rt_P1 = [t/system1.Porb for t in rt]
  n_l_m = [(n,l,m) for (n,l,m,s), i, j in common_nlms]

  if opts.verbose: print "generating plots"

  ###
  if opts.res: # plot the residuals
    if opts.verbose: print "\tres"
    fig, ax = nmp.amp_plot(rt_P1, r, n_l_m=n_l_m)
    ax.set_xlabel(r"$t/P_{\mathrm{orb}}^1$")
    ax.set_ylabel(r"$|r_i|$")

    ax.grid(opts.grid)
   
    if not opts.no_legend:
      ax.legend()

    figname = "%s/%sres.png" % (opts.output_dir, opts.tag)
    if opts.verbose: print "saving %s" % figname
    fig.savefig(figname)
    nmp.plt.close(fig)

  ###
  if opts.cross_correlation: # compute the cross-correlation between the data. This should be like function of the relative -time between sample points
    if opts.verbose: print "\tcross_correlation"
    tcc, fcc = liapunov.nmode_cross_correlation(rq1, rq2)

    rfreq1 = nmu.fft_freq(rt_P1)

    # time domain
    fig, ax = nmp.amp_phs_plot(rt_P1, tcc, n_l_m=n_l_m)
    ax.set_xlabel(r"$t/P_{\mathrm{orb}}^1$")
    ax.set_ylabel(r"$cross-correlation$")

    ax.grid(opts.grid)
    if opts.no_legend:
      ax.legend()

    figname = "%s/%scross_correlation.time.png" % (opts.output_dir, opts.tag)
    if opts.verbose: print "saving %s" % figname
    fig.savefig(figname)
    nmp.plt.close(fig)

    # freq domain
    fig, ax = nmp.amp_phs_plot(rfreq1, fcc, n_l_m=n_l_m)
    ax.set_xlabel(r"$f/f_{\mathrm{orb}}^1$")
    ax.set_ylabel(r"$cross-correlation$")

    ax.grid(opts.grid)
    if opts.no_legend:
      ax.legend()

    figname = "%s/%scross_correlation.freq.png" % (opts.output_dir, opts.tag)
    if opts.verbose: print "saving %s" % figname
    fig.savefig(figname)
    nmp.plt.close(fig)



  ###
#  if opts.amp_res: # plot the amplitudes and residuals
#    print "WRITE amp_res"
#    # in lib/liapunov.py , as a specialty projection
