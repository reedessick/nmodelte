#!/usr/local/bin/python

usage="""an executable to load, interpret and plot integration data. "p" is for "plotting" """

import numpy as np

import nmode_plotting as nm_p
import nmode_utils as nm_u
import nmode_state as nm_s
from nmode_plotting import plt

from optparse import *

####################################################################################################
#
#
#                        Parse input options
#
#
####################################################################################################

parser=OptionParser(usage=usage)

parser.add_option("-v", "--verbose", dest="verbose", action="store_true", help="print updates on the progress of the script")
parser.add_option("", "--eig-verbose", action="store_true", help="print updates when computing eigenvalues. Also, prints eigenvectors to the terminal")
parser.add_option("-t", "--tag", dest="tag", default="", type="string", help="a tag that is appended to the end of figures' filenames.")
parser.add_option("-d", "--downsample", default=False, type="int", help="plot only 1 point out of this many")

parser.add_option("-l", "--logfilename", default=False, type="string")

# time domain options
parser.add_option("-F", "--filename", dest="filename", default=False, type="string", help="the name of the file you want to plot")

parser.add_option("", "--amp-time", action="store_true")
parser.add_option("", "--phs-time", action="store_true")
parser.add_option("", "--amp-phs-time", action="store_true", help="generate time-domain data: amplitude and phase.")
parser.add_option("", "--real-time", action="store_true")
parser.add_option("", "--imag-time", action="store_true")
parser.add_option("", "--real-imag-time", action="store_true", help="generate time-domain data: real and imaginary parts.")
parser.add_option("", "--multi-gen-amp-time", action="store_true")
parser.add_option("", "--Hns-time", action="store_true")

parser.add_option("", "--lin-amp-time", action="store_true", help="sets all amplitude scales to \"linear\" in time-domain plots. Also controls Hns-time's yscale")

### phase portraits in the time domain
parser.add_option("", "--amp-portrait", action="store_true", help="generate a phase portrait of the amplitude of each mode: |x|=|q|")
parser.add_option("", "--multi-gen-amp-portrait", action="store_true", help="generate a separate amp_portrait for each generation")

parser.add_option("", "--real-imag-portrait", action="store_true", help="generate a phase portrait of the Real part vs the Imaginary part of each mode")
parser.add_option("", "--multi-gen-real-imag-portrait", action="store_true", help="generate a separate real_imag_portrait for each generation")

### eigenvalue maps
parser.add_option("", "--max-eig-Lamp-phs", default=False, action="store_true", help="plot the maximum eigenvalue of the system over time with linear amplitdue mode shapes")
parser.add_option("", "--max-eig-amp-phs", default=False, action="store_true", help="plot the maximum eigenvalue of the system over time with logarithmic amplitude mode shapes")
parser.add_option("", "--max-eig-real-imag", default=False, action="store_true", help="plot the maximum eigenvalue of the system over time with linear real/imag mode shapes")
parser.add_option("", "--eigvals", default=False, action="store_true", help="plot trajectories of eigenvalues of system in complex plane over time")
parser.add_option("", "--eig-minr", default=False, type="float")
parser.add_option("", "--eig-maxr", default=False, type="float")
parser.add_option("", "--eig-mini", default=False, type="float")
parser.add_option("", "--eig-maxi", default=False, type="float")

### growth rates
parser.add_option("", "--growth-rates", default=False, action="store_true", help="compute logarithmic derivatives of integration variable and plot real,imag parts vs. time")

parser.add_option("", "--log-time", action="store_true", help="plot time-domain data with logarithmic time axis.")
parser.add_option("", "--lin-time", action="store_true", help="generate time-domain figures with linear time axis.")

parser.add_option("", "--time-ymin", default=False, type="float", help="the lower limit on the y-axis in time domain figures.")
parser.add_option("", "--time-ymax", default=False, type="float", help="the upper limit on the y-axis in time domain figures.")
parser.add_option("", "--time-xmin", default=False, type="float", help="the lower limit on the x-axis in time domain figures.")
parser.add_option("", "--time-xmax", default=False, type="float", help="the upper limit on the x-axis in time domain figures.")

parser.add_option("", "--tcurrent", dest="tc", default=False, type="string", help="the current time variable type: x or q")

# frequency domain options
parser.add_option("", "--freqfilename", default=False, type="string", help="write frequency data into this file.")

parser.add_option("", "--amp-freq", action="store_true")
parser.add_option("", "--phs-freq", action="store_true")
parser.add_option("", "--amp-phs-freq", action="store_true", help="generate frequency-domain data: amplitude and phase")
parser.add_option("", "--real-freq", action="store_true")
parser.add_option("", "--imag-freq", action="store_true")
parser.add_option("", "--real-imag-freq", action="store_true", help="generate frequency-domain data: real and imaginary parts.")
parser.add_option("", "--multi-gen-amp-freq", action="store_true")

parser.add_option("", "--lin-amp-freq", action="store_true", help="sets all amplitude scales to \"linear\" in time-domain plots.")

parser.add_option("", "--log-freq", action="store_true", help="plot frequency-domain data with logarthmic frequency axis.")
parser.add_option("", "--lin-freq", action="store_true", help="generate frequency-domain figures with linear frequency axis.")

parser.add_option("", "--freq-ymin", default=False, type="float", help="the lower limit on the y-axis in freq domain figures.")
parser.add_option("", "--freq-ymax", default=False, type="float", help="the upper limit on the y-axis in freq domain figures.")
parser.add_option("", "--freq-xmin", default=False, type="float", help="the lower limit on the x-axis in freq domain figures.")
parser.add_option("", "--freq-xmax", default=False, type="float", help="the upper limit on the x-axis in freq domain figures.")

parser.add_option("", "--fcurrent", dest="fc", default=False, type="string", help="the current freq variable type: x or q . You can also specify A if using the FFT of the amplitude's abs value.")

# legend options
parser.add_option("", "--no-legend", action="store_true", default=False, help="do not plot a legend on any figure")
parser.add_option("", "--legend-col", default=1, type="int", help="the number of columns for your legend")
parser.add_option("", "--legend-font-size", default = 12, type="int", help="the font size for text in the legend.")
parser.add_option("", "--legend-loc", default="best", type="string", help="location for legend.")

# plot only a few mode numbers
parser.add_option("", "--mode-nums", default=False, type="string", help="a string containing all the mode numbers you want to plot.")

# put grids on all plots
parser.add_option("-g", "--grid", default=False, action="store_true", help="put grid-lines on all plots.")

opts, args = parser.parse_args()

##################################################

### make time domain plots?
if (opts.amp_time or opts.phs_time or opts.amp_phs_time or opts.real_time or opts.imag_time or opts.real_imag_time or opts.multi_gen_amp_time or opts.Hns_time or opts.growth_rates):
  time_domain = True
else:
  time_domain = False

### make freq domain plots?
if (opts.amp_freq or opts.phs_freq or opts.amp_phs_freq or opts.real_freq or opts.imag_freq or opts.real_imag_freq or opts.multi_gen_amp_freq):
  freq_domain = True
else:
  freq_domain = False

### make multi-generation plots?
if (opts.multi_gen_amp_time or opts.multi_gen_amp_freq or opts.multi_gen_amp_portrait or opts.multi_gen_real_imag_portrait):
  multi_gen = True
else:
  multi_gen = False

### make phase portraits?
if (opts.amp_portrait or opts.real_imag_portrait or opts.multi_gen_amp_portrait or opts.multi_gen_real_imag_portrait):
  portrait = True
else:
  portrait = False

### compute Hns?
if opts.Hns_time: 
  compute_Hns = True
else:
  compute_Hns = False

### make eigenvalue plots?
if (opts.max_eig_amp_phs or opts.max_eig_Lamp_phs or opts.max_eig_real_imag or opts.eigvals):
  max_eig = True
else:
  max_eig = False

### other stuff
if opts.tag != "":
  opts.tag = "." + opts.tag

if opts.mode_nums:
  opts.mode_nums = [int(l) for l in opts.mode_nums.split()]

### ensure we have all the required data
if (multi_gen or portrait or compute_Hns or max_eig) and (not opts.logfilename):
  opts.logfilename = raw_input("logfilename = ")

if time_domain and (not (opts.lin_time or opts.log_time) ):
  opts.lin_time = (raw_input("plot lin_time (y/n)?") == "y")
  opts.log_time = (raw_input("plot log_time (y/n)?") == "y")
  if not (opts.lin_time or opts.log_time):
    import sys
    sys.exit("nothing to do")

if freq_domain and (not (opts.lin_freq or opts.log_freq) ):
  opts.lin_freq = (raw_input("plot lin_freq (y/n)?") == "y")
  opts.log_freq = (raw_input("plot log_freq (y/n)?") == "y")
  if not (opts.lin_freq or opts.log_freq):
    import sys
    sys.exit("nothing to do")

if (portrait or compute_Hns or max_eig) and (not time_domain): 
  time_domain = True # set this to be true after you check for time axis scaling

## variable types
if time_domain and (not opts.tc):
  opts.tc = raw_input("current time variable type = ")

if freq_domain and (not opts.fc):
  opts.fc = raw_input("current freq variable type = ")

####################################################################################################
#
#
#                         params and network specifics
#
#
####################################################################################################
n_l_m = False
if opts.logfilename:
  if opts.verbose: print "loading system from "+opts.logfilename
  system = nm_u.load_log(opts.logfilename)
  n_l_m = system.network.nlm

  if multi_gen:
    if opts.verbose: print "determining generational structure"
    gens, _ = system.network.gens()
#    gi = system.network.find_G0()
#    gens = []
#    N_m = len(n_l_m)
#    Nm = len(gi)
#    gens.append(gi)
#    while N_m > Nm: # continue until all modes are included
#      gi, _ = system.network.find_Gip1(gi)
#      Nm += len(gi)
#      gens.append(gi)

####################################################################################################
#
#
#                            time domain
#
#
####################################################################################################

if time_domain:

  if opts.verbose: print "building time-domain figures"

  if not opts.filename:
    opts.filename = raw_input("filename = ")
  if opts.verbose: print "loading "+opts.filename

  t_P, q, N_m = nm_u.load_out(opts.filename, tmin=opts.time_xmin, tmax=opts.time_xmax, downsample=opts.downsample)

  if compute_Hns:
    if opts.tc == "q":
      Hns = nm_s.Hns_q(q, system.network, Eo=1.)
    elif opts.tc=="x":
      Hns = nm_s.Hns_x(t_P, q, system, Eo=1.)

  ################################################
  ### Hns_time
  if opts.Hns_time:
    if opts.verbose: print "\tHns_time"
    fig = plt.figure()
    ax = plt.subplot(1,1,1)

    ax.plot(t_P, Hns, label=r"$H_{\mathrm{NS}}$")

    ax.set_ylabel(r"$H_{\mathrm{NS}}$")    
    ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')
    if not opts.lin_amp_time:
      ax.set_yscale('log')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_ymin:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], opts.time_ymin, _ax[3]])
    if opts.time_ymax:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], _ax[2], opts.time_ymax])
    if opts.time_xmin:
      _ax = ax.axis()
      ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.time_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".LHns.lin_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".Hns.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".LHns.log_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".Hns.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.time_xmin:
        _ax = ax.axis()
        ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.time_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ### amp_time
  if opts.amp_time:
    if opts.verbose: print "\tamp_time"
    fig, ax = nm_p.amp_plot(t_P, q, n_l_m=n_l_m, mode_nums=opts.mode_nums)
    
    ax.set_ylabel(r'$|'+opts.tc+'_i|$')
    ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')
    if not opts.lin_amp_time:
      ax.set_yscale('log')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_ymin:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], opts.time_ymin, _ax[3]])
    if opts.time_ymax:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], _ax[2], opts.time_ymax])
    if opts.time_xmin:
      _ax = ax.axis()
      ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.time_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".Lamp.lin_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".amp.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".Lamp.log_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".amp.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.time_xmin:
        _ax = ax.axis()
        ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.time_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### phs_time
  if opts.phs_time:
    if opts.verbose: print "\tphs_time"
    fig, ax = nm_p.phs_plot(t_P, q, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax.set_ylabel(r'$\mathrm{arg}['+opts.tc+'_i]/2\pi$')
    ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_xmin:
      _ax = ax.axis()
      ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.time_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_time:
      figname = opts.filename + ".phs.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      figname = opts.filename + ".phs.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.time_xmin:
        _ax = ax.axis()
        ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.time_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)
    
  ################################################
  ### amp_phs_time
  if opts.amp_phs_time:
    if opts.verbose: print "\tamp_phs_time"
    fig, ax1, ax2 = nm_p.amp_phs_plot(t_P, q, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r'$|'+opts.tc+'_i|$')
    ax2.set_xlabel(r'$t/P_{\mathrm{orb}}$')
    ax2.set_ylabel(r'$\mathrm{arg}['+opts.tc+'_i]/2\pi$')

    if not opts.lin_amp_time:
      ax1.set_yscale('log')

    if not opts.no_legend:
      ax1.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_ymin:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], opts.time_ymin, ax[3]])
    if opts.time_ymax:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], ax[2], opts.time_ymax])
    if opts.time_xmin:
      ax = ax1.axis()
      ax1.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
    if opts.time_xmax:
      ax = ax1.axis()
      ax1.axis([ax[0], opts.time_xmax, ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([ax[0], opts.time_xmax, ax[2], ax[3]])

    # grid lines
    ax1.grid(opts.grid)
    ax2.grid(opts.grid)

    if opts.lin_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".Lamp_phs.lin_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".amp_phs.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".Lamp_phs.log_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".amp_phs.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax1.set_xscale('log')
      ax2.set_xscale('log')
      if opts.time_xmin:
        ax = ax1.axis()
        ax1.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
      if opts.time_xmax:
        ax = ax1.axis()
        ax1.axis([ax[0], opts.time_xmax, ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([ax[0], opts.time_xmax, ax[2], ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### real_time
  if opts.real_time:
    if opts.verbose: print "\treal_time"
    fig, ax = nm_p.real_plot(t_P, q, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax.set_ylabel(r'$\mathbb{R}\{'+opts.tc+'_i\}$')
    ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_ymin:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], opts.time_ymin, _ax[3]])
    if opts.time_ymax:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], _ax[2], opts.time_ymax])
    if opts.time_xmin:
      _ax = ax.axis()
      ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.time_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_time:
      figname = opts.filename + ".real.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      figname = opts.filename + ".real.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.time_xmin:
        _ax = ax.axis()
        ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.time_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)


  ################################################
  ### imag_time
  if opts.imag_time:
    if opts.verbose: print "\timag_time"
    fig, ax = nm_p.imag_plot(t_P, q, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax.set_ylabel(r'$\mathbb{I}\{'+opts.tc+'_i\}$')
    ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_ymin:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], opts.time_ymin, _ax[3]])
    if opts.time_ymax:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], _ax[2], opts.time_ymax])
    if opts.time_xmin:
      _ax = ax.axis()
      ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.time_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_time:
      figname = opts.filename + ".imag.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      figname = opts.filename + ".imag.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.time_xmin:
        _ax = ax.axis()
        ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.time_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### real_imag_time
  if opts.real_imag_time:
    if opts.verbose: print "\treal_imag_time"
    fig, ax1, ax2 = nm_p.real_imag_plot(t_P, q, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r'$\mathbb{R}\{'+opts.tc+'_i\}$')
    ax2.set_xlabel(r'$t/P_{\mathrm{orb}}$')
    ax2.set_ylabel(r'$\mathbb{I}\{'+opts.tc+'_i\}$')

    if not opts.no_legend:
      ax1.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_ymin:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], opts.time_ymin, ax[3]])
      ax = ax2.axis()
      ax2.axis([ax[0], ax[1], opts.time_ymin, ax[3]])
    if opts.time_ymax:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], ax[2], opts.time_ymax])
      ax = ax2.axis()
      ax2.axis([ax[0], ax[1], ax[2], opts.time_ymax])
    if opts.time_xmin:
      ax = ax1.axis()
      ax1.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
    if opts.time_xmax:
      ax = ax1.axis()
      ax1.axis([ax[0], opts.time_xmax, ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([ax[0], opts.time_xmax, ax[2], ax[3]])

    # grid lines
    ax1.grid(opts.grid)
    ax2.grid(opts.grid)

    if opts.lin_time:
      figname = opts.filename + ".real_imag.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      figname = opts.filename + ".real_imag.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax1.set_xscale('log')
      ax2.set_xscale('log')
      if opts.time_xmin:
        ax = ax1.axis()
        ax1.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
      if opts.time_xmax:
        ax = ax1.axis()
        ax1.axis([ax[0], opts.time_xmax, ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([ax[0], opts.time_xmax, ax[2], ax[3]])

      plt.savefig(figname)

    plt.close(fig)

##################################################
#
#      multi-gen plot time-domain plots
#
##################################################
  if opts.multi_gen_amp_time:
    if opts.verbose: print "\tmulti_gen_amp_time"
    fig, axs = nm_p.multi_gen_amp_plot(t_P, q, gens, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    for ax in axs:
      ax.set_ylabel(r'$|'+opts.tc+'_i|$')
      if not opts.lin_amp_time:
        ax.set_yscale('log')

      if not opts.no_legend:
        ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

      if opts.time_ymin:
        _ax = ax.axis()
        ax.axis([_ax[0], _ax[1], opts.time_ymin, _ax[3]])
      if opts.time_ymax:
        _ax = ax.axis()
        ax.axis([_ax[0], _ax[1], _ax[2], opts.time_ymax])
      if opts.time_xmin:
        _ax = ax.axis()
        ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.time_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

      # grid lines
      ax.grid(opts.grid)

    ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')

    if opts.lin_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".multi-gen_Lamp.lin_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".multi-gen_amp.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      if opts.lin_amp_time:
        figname = opts.filename + ".multi-gen_Lamp.log_time" + opts.tag + ".png"
      else:
        figname = opts.filename + ".multi-gen_amp.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      for ax in axs:
        ax.set_xscale('log')
        if opts.time_xmin:
          _ax = ax.axis()
          ax.axis([opts.time_xmin, _ax[1], _ax[2], _ax[3]])
        if opts.time_xmax:
          _ax = ax.axis()
          ax.axis([_ax[0], opts.time_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

##################################################
#
#      phase portraits in the time-domain
#
##################################################
  
  ################################################
  ### amplitude phase portrait
  if opts.amp_portrait:
    if opts.verbose: print "\tamp_portrait"
    fig, ax = nm_p.amp_phase_portrait([t*system.Porb for t in t_P], q, system, opts.tc, n_l_m=n_l_m, mode_nums=opts.mode_nums, verbose=opts.verbose)

    ax.set_xlabel(r'$|'+opts.tc+'_i|$')
    ax.set_ylabel(r'$\partial_t |'+opts.tc+'_i| P_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    ax.grid(opts.grid)

    figname = opts.filename+".amp_portrait"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### mulit-gen amplitude phase portrait
  if opts.multi_gen_amp_portrait:
    if opts.verbose: print "\tmulti_gen_amp_portrait"
    fig, axs = nm_p.multi_gen_amp_phase_portrait([t*system.Porb for t in t_P], q, system, opts.tc, gens, n_l_m=n_l_m, mode_nums=opts.mode_nums, verbose=opts.verbose)

    for ax in axs:
      ax.set_ylabel(r'$\partial_t |'+opts.tc+'_i| P_{\mathrm{orb}}$')

      if not opts.no_legend:
        ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

      ax.grid(opts.grid)

    ax.set_xlabel(r'$|'+opts.tc+'_i|$')

    figname = opts.filename+".multi-gen_amp_portrait"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### real-imag phase portraits
  if opts.real_imag_portrait:
    if opts.verbose: print "\treal_imag_portrait"
    fig, ax = nm_p.real_imag_phase_portrait(q, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax.set_xlabel(r'$\mathbb{R}\{'+opts.tc+'_i\}$')
    ax.set_ylabel(r'$\mathbb{I}\{'+opts.tc+'_i\}$')
    
    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    ax.grid(opts.grid)

    figname = opts.filename+".real_imag_portrait"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### multi-gen real-imag phase portraits
  if opts.multi_gen_real_imag_portrait:
    if opts.verbose: print "\tmulti_gen_real_imag_portrait"
    fig, axs = nm_p.multi_gen_real_imag_phase_portrait(q, gens, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    for ax in axs:
      ax.set_ylabel(r'$\mathbb{I}\{'+opts.tc+'_i\}$')

      if not opts.no_legend:
        ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

      # grid lines
      ax.grid(opts.grid)

    ax.set_xlabel(r'$\mathbb{R}\{'+opts.tc+'\}$')

    figname = opts.filename+".multi-gen_real_imag_portrait"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    plt.savefig(figname)

    plt.close(fig)

##################################################
#
#                eigevalue maps
#
##################################################
### growth rates (logarithmic derivatives)
  if opts.growth_rates:
    if opts.verbose: print "growth_rates"
    fig, ax1, ax2 = nm_p.real_imag_plot(t_P, nm_p.growth_rates([t*system.Porb for t in t_P], q, system, opts.tc, verbose=opts.verbose), n_l_m=n_l_m, mode_nums=opts.mode_nums)
  
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r'$\mathbb{R}\{\partial_{t} '+opts.tc+'_i / '+opts.tc+'_i \} \cdot P_{\mathrm{orb}}$')
    ax2.set_xlabel(r'$t/P_{\mathrm{orb}}$')
    ax2.set_ylabel(r'$\mathbb{I}\{\partial_{t} '+opts.tc+'_i / '+opts.tc+'_i \} / \Omega_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax1.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.time_ymin:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], opts.time_ymin, ax[3]])
      ax = ax2.axis()
      ax2.axis([ax[0], ax[1], opts.time_ymin, ax[3]])
    if opts.time_ymax:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], ax[2], opts.time_ymax])
      ax = ax2.axis()
      ax2.axis([ax[0], ax[1], ax[2], opts.time_ymax])
    if opts.time_xmin:
      ax = ax1.axis()
      ax1.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
    if opts.time_xmax:
      ax = ax1.axis()
      ax1.axis([ax[0], opts.time_xmax, ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([ax[0], opts.time_xmax, ax[2], ax[3]])

    # grid lines
    ax1.grid(opts.grid)
    ax2.grid(opts.grid)

    if opts.lin_time:
      figname = opts.filename + ".growth_rates.lin_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_time:
      figname = opts.filename + ".growth_rates.log_time" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax1.set_xscale('log')
      ax2.set_xscale('log')
      if opts.time_xmin:
        ax = ax1.axis()
        ax1.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([opts.time_xmin, ax[1], ax[2], ax[3]])
      if opts.time_xmax:
        ax = ax1.axis()
        ax1.axis([ax[0], opts.time_xmax, ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([ax[0], opts.time_xmax, ax[2], ax[3]])

      plt.savefig(figname)

    plt.close(fig)
 
### linearization and diagonalizatioin
  if max_eig:

    if not (opts.eigvals):
      if opts.verbose: print "computing eigenvalues, eigenvectors"
      if opts.tc == "x":
        eigvals, eigvecs = nm_p.nmode_eig_x(t_P, q, system, max_eig=True, verbose=opts.eig_verbose)
      elif opts.tc == "q":
        eigvals, eigvecs = nm_p.nmode_eig_q(q, system, max_eig=True, verbose=opts.eig_verbose)

#      if opts.eig_verbose: print eigvals, eigvecs

    elif (opts.max_eig_amp_phs or opts.max_eig_real_imag):
      if opts.verbose: print "computing eigenvalues, eigenvectors"
      if opts.tc == "x":
        eigvals, eigvecs = nm_p.nmode_eig_x(t_P, q, system, max_eig=False, verbose=opts.eig_verbose)
      elif opts.tc == "q":
        eigvals, eigvecs = nm_p.nmode_eig_q(q, system, max_eig=False, verbose=opts.eig_verbose)

#      if opts.eig_verbose: print eigvals, eigvecs

    else: # only eigvals needed
      if opts.verbose: print "computing eigenvalues"
      if opts.tc == "x":
        eigvals = nm_p.nmode_eigval_x(t_P, q, system, max_eig=False, verbose=opts.eig_verbose)
      elif opts.tc == "q":
        eigvals = nm_p.nmode_eigval_q(q, system, max_eig=False, verbose=opts.eig_verbose)

#      if opts.eig_verbose: print eigvals

    #######################
    if opts.eigvals:
      if opts.verbose: print "\teigvals"

      fig, ax1, ax2 = nm_p.eigval_plot(t_P, eigvals, system)

      ax1.grid(opts.grid)
      ax2.grid(opts.grid)

      if opts.eig_minr:
        ax1.set_ylim(ymin=opts.eig_minr)
      if opts.eig_maxr:
        ax1.set_ylim(ymax=opts.eig_maxr)
      if opts.eig_mini:
        ax2.set_ylim(ymin=opts.eig_mini)
      if opts.eig_maxi:
        ax2.set_ylim(ymax=opts.eig_maxi)

      if opts.time_xmin:
        ax1.set_xlim(xmin=opts.time_xmin)
        ax2.set_xlim(xmin=opts.time_xmin)
      if opts.time_xmax:
        ax1.set_xlim(xmax=opts.time_xmax)
        ax2.set_xlim(xmax=opts.time_xmax)

      figname = opts.filename+".eigvals"+opts.tag+".png"
      if opts.verbose: print "saving "+figname
      nm_p.plt.savefig(figname)
      nm_p.plt.close(fig)

      if opts.max_eig_amp_phs or opts.max_eig_real_imag: # set up eigvals, eigvecs for the rest of the methods
        eigenvalues = []
        eigenvectors = []
        for ind, eva in enumerate(eigvals):
          ev = [(np.real(e), ev_ind) for ev_ind, e in enumerate(eva)]
          ev.sort(key=lambda l: l[0], reverse=True)
          max_ev_ind = ev[0][1]
          eigenvalues.append(eva[max_ev_ind])
          eigenvectors.append(eigvecs[ind][:,max_ev_ind])
        eigvals = eigenvalues
        eigvecs = eigenvectors
 
    #######################
    if opts.max_eig_Lamp_phs or opts.max_eig_amp_phs:

      vecs_amp = np.zeros((len(t_P),N_m))
      vecs_phs = np.zeros((len(t_P),N_m))

      for ind in range(len(t_P)):
        v = eigvecs[ind]
        for modeNo in range(N_m):
          r, i = v[2*modeNo:2*modeNo+2]
          vecs_amp[ind][modeNo] = nm_u.amp(r,i)
          vecs_phs[ind][modeNo] = nm_u.phs(r,i)/(2*np.pi)

      if opts.max_eig_Lamp_phs:
        if opts.verbose: print "\tmax_eig_Lamp_phs"

        fig, [ax1,ax1p,ax2,ax3], [cb2,cb3] = nm_p.eig_plot(t_P, eigvals, vecs_amp, vecs_phs, system, cmap=nm_p.plt.cm.YlOrRd)

        if opts.eig_minr:
          ax1.set_ylim(ymin=opts.eig_minr)
        if opts.eig_maxr:
          ax1.set_ylim(ymax=opts.eig_maxr)
        if opts.eig_mini:
          ax1p.set_ylim(ymin=opts.eig_mini)
        if opts.eig_maxi:
          ax1p.set_ylim(ymax=opts.eig_maxi)

        if opts.time_xmin:
          ax1.set_xlim(xmin=opts.time_xmin)
        if opts.time_xmax:
          ax1.set_xlim(xmax=opts.time_xmax)

        cb2.set_label(r'$|'+opts.tc+'_i|$')
        cb3.set_label(r'$\mathrm{arg}['+opts.tc+'_i]/2\pi$')
        ax1.grid(opts.grid)

        figname = opts.filename+".max_eig_Lamp_phs"+opts.tag+".png"
        if opts.verbose: print "saving "+figname
        nm_p.plt.savefig(figname)
        nm_p.plt.close(fig)

      if opts.max_eig_amp_phs:
        if opts.verbose: print "\tmax_eig_amp_phs"

        # convert vecs_amp --> log(vecs_amp) for scaling purposes

        fig, [ax1,ax1p,ax2,ax3], [cb2,cb3] = nm_p.eig_plot(t_P, eigvals, np.log10(vecs_amp), vecs_phs, system, cmap=nm_p.plt.cm.YlOrRd)

        if opts.eig_minr:
          ax1.set_ylim(ymin=opts.eig_minr)
        if opts.eig_maxr:
          ax1.set_ylim(ymax=opts.eig_maxr)
        if opts.eig_mini:
          ax1p.set_ylim(ymin=opts.eig_mini)
        if opts.eig_maxi:
          ax1p.set_ylim(ymax=opts.eig_maxi)

        if opts.time_xmin:
          ax1.set_xlim(xmin=opts.time_xmin)
        if opts.time_xmax:
          ax1.set_xlim(xmax=opts.time_xmax)

        cb2.set_label(r'$\log_{10} |'+opts.tc+'_i|$')
        cb3.set_label(r'$\mathrm{arg}['+opts.tc+'_i]/2\pi$')
        ax1.grid(opts.grid)

        figname = opts.filename+".max_eig_amp_phs"+opts.tag+".png"
        if opts.verbose: print "saving "+figname
        nm_p.plt.savefig(figname)
        nm_p.plt.close(fig)

    #######################
    if opts.max_eig_real_imag:
      if opts.verbose: print "\tmax_eig_real_imag"

      vecs_real = np.zeros((len(t_P),N_m))
      vecs_imag = np.zeros((len(t_P),N_m))

      for ind in range(len(t_P)):
        v = eigvecs[ind]
        for modeNo in range(N_m):
          vecs_real[ind][modeNo] = v[2*modeNo]
          vecs_imag[ind][modeNo] = v[2*modeNo+1]

      fig, [ax1,ax1p,ax2,ax3], [cb2,cb3] = nm_p.eig_plot(t_P, eigvals, vecs_real, vecs_imag, system, cmap=nm_p.plt.cm.coolwarm)

      if opts.eig_minr:
        ax1.set_ylim(ymin=opts.eig_minr)
      if opts.eig_maxr:
        ax1.set_ylim(ymax=opts.eig_maxr)
      if opts.eig_mini:
        ax1p.set_ylim(ymin=opts.eig_mini)
      if opts.eig_maxi:
        ax1p.set_ylim(ymax=opts.eig_maxi)

      if opts.time_xmin:
        ax1.set_xlim(xmin=opts.time_xmin)
      if opts.time_xmax:
        ax1.set_xlim(xmax=opts.time_xmax)

      cb2.set_label(r'$\mathbb{R}\{'+opts.tc+'_i\}$')
      cb3.set_label(r'$\mathbb{I}\{'+opts.tc+'_i\}$')
      ax1.grid(opts.grid)

      figname = opts.filename+".max_eig_real_imag"+opts.tag+".png"
      if opts.verbose: print "saving "+figname
      nm_p.plt.savefig(figname)
      nm_p.plt.close(fig)


####################################################################################################
#
#
#                            frequency domain
#
#
####################################################################################################

if freq_domain:

  if opts.verbose: print "building freq-domain figures"

  if not opts.freqfilename:
    opts.freqfilename = raw_input("freqfilename = ")
  if opts.verbose: print "loading "+opts.freqfilename

  freq, fq, N_m = nm_u.load_out(opts.freqfilename, tmin=opts.freq_xmin, tmax=opts.freq_xmax, downsample=opts.downsample)

  ################################################
  ### amp_freq
  if opts.amp_freq:
    if opts.verbose: print "\tamp_freq"
    fig, ax = nm_p.amp_plot(freq, fq, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax.set_ylabel(r'$|\widetilde{'+opts.fc+'}_i|$')
    ax.set_xlabel(r'$f/f_{\mathrm{orb}}$')

    if not opts.lin_amp_freq:
      ax.set_yscale('log')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.freq_ymin:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], opts.freq_ymin, _ax[3]])
    if opts.freq_ymax:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], _ax[2], opts.freq_ymax])
    if opts.freq_xmin:
      _ax = ax.axis()
      ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.freq_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_freq:
      if opts.lin_amp_freq:
        figname = opts.freqfilename + ".Lamp.lin_freq" + opts.tag + ".png"
      else:
        figname = opts.freqfilename + ".amp.lin_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_freq:
      if opts.lin_amp_freq:
        figname = opts.freqfilename + ".Lamp.log_freq" + opts.tag + ".png"
      else:
        figname = opts.freqfilename + ".amp.log_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.freq_xmin:
        _ax = ax.axis()
        ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.freq_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### phs_freq  
  if opts.phs_freq:
    if opts.verbose: print "\tphs_freq"
    fig, ax = nm_p.phs_plot(freq, fq, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax.set_ylabel(r'$\mathrm{arg}[\widetilde{'+opts.fc+'}_i]/2\pi$')
    ax.set_xlabel(r'$f/f_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.freq_xmin:
      _ax = ax.axis()
      ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.freq_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_freq:
      figname = opts.freqfilename + ".phs.lin_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_freq:
      figname = opts.freqfilename + ".phs.log_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.freq_xmin:
        _ax = ax.axis()
        ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.freq_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### amp_phs_freq
  if opts.amp_phs_freq:
    if opts.verbose: print "\tamp_phs_freq"
    fig, ax1, ax2 = nm_p.amp_phs_plot(freq, fq, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r'$|\widetilde{'+opts.fc+'}_i|$')
    ax2.set_xlabel(r'$f/f_{\mathrm{orb}}$')
    ax2.set_ylabel(r'$\mathrm{arg}[\widetilde{'+opts.fc+'}_i]/2\pi$')

    if not opts.lin_amp_freq:
      ax1.set_yscale('log')

    if not opts.no_legend:
      ax1.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.freq_ymin:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], opts.freq_ymin, ax[3]])
    if opts.freq_ymax:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], ax[2], opts.freq_ymax])
    if opts.freq_xmin:
      ax = ax1.axis()
      ax1.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
    if opts.freq_xmax:
      ax = ax1.axis()
      ax1.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])

    # grid lines
    ax1.grid(opts.grid)
    ax2.grid(opts.grid)

    if opts.lin_freq:
      if opts.lin_amp_freq:
        figname = opts.freqfilename + ".Lamp_phs.lin_freq" + opts.tag + ".png"
      else:
        figname = opts.freqfilename + ".amp_phs.lin_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_freq:
      if opts.lin_amp_freq:
        figname = opts.freqfilename + ".Lamp_phs.log_freq" + opts.tag + ".png"
      else:
        figname = opts.freqfilename + ".amp_phs.log_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax1.set_xscale('log')
      ax2.set_xscale('log')
      if opts.freq_xmin:
        ax = ax1.axis()
        ax1.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
      if opts.freq_xmax:
        ax = ax1.axis()
        ax1.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### real_freq
  if opts.real_freq:
    if opts.verbose: print "\treal_freq"
    fig, ax = nm_p.real_plot(freq, fq, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax1.set_ylabel(r'$\mathbb{R}\{\widetilde{'+opts.fc+'}_i\}$')
    ax.set_xlabel(r'$f/f_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.freq_ymin:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], opts.freq_ymin, _ax[3]])
    if opts.freq_ymax:
      _ax = ax.axis()
      ax.axis([_ax[0], _ax[1], _ax[2], opts.freq_ymax])
    if opts.freq_xmin:
      _ax = ax.axis()
      ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.freq_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_freq:
      figname = opts.freqfilename + ".real.lin_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_freq:
      figname = opts.freqfilename + ".real.log_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.freq_xmin:
        _ax = ax.axis()
        ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.freq_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### imag_freq  
  if opts.imag_freq:
    if opts.verbose: print "\timag_freq"
    fig, ax = nm_p.imag_plot(freq, fq, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    ax.set_ylabel(r'$\mathbb{I}\{\widetilde{'+opts.fc+'}_i\}$')
    ax.set_xlabel(r'$f/f_{\mathrm{orb}}$')

    if not opts.no_legend:
      ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.freq_xmin:
      _ax = ax.axis()
      ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
    if opts.freq_xmax:
      _ax = ax.axis()
      ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

    # grid lines
    ax.grid(opts.grid)

    if opts.lin_freq:
      figname = opts.freqfilename + ".imag.lin_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_freq:
      figname = opts.freqfilename + ".imag.log_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax.set_xscale('log')
      if opts.freq_xmin:
        _ax = ax.axis()
        ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.freq_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

  ################################################
  ### real_imag_freq
  if opts.real_imag_freq:
    if opts.verbose: print "\treal_imag_freq"
    fig, ax1, ax2 = nm_p.real_imag_plot(freq, fq, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylabel(r'$\mathbb{R}\{\widetilde{'+opts.fc+'}_i\}$')
    ax2.set_xlabel(r'$f/f_{\mathrm{orb}}$')
    ax2.set_ylabel(r'$\mathbb{I}\{\widetilde{'+opts.fc+'}_i\}$')

    if not opts.no_legend:
      ax1.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

    if opts.freq_ymin:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], opts.freq_ymin, ax[3]])
    if opts.freq_ymax:
      ax = ax1.axis()
      ax1.axis([ax[0], ax[1], ax[2], opts.freq_ymax])
    if opts.freq_xmin:
      ax = ax1.axis()
      ax1.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
    if opts.freq_xmax:
      ax = ax1.axis()
      ax1.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])
      ax = ax2.axis()
      ax2.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])

    # grid lines
    ax1.grid(opts.grid)
    ax2.grid(opts.grid)

    if opts.lin_freq:
      figname = opts.freqfilename + ".real_imag.lin_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_freq:
      figname = opts.freqfilename + ".real_imag.log_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      ax1.set_xscale('log')
      ax2.set_xscale('log')
      if opts.freq_xmin:
        ax = ax1.axis()
        ax1.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([opts.freq_xmin, ax[1], ax[2], ax[3]])
      if opts.freq_xmax:
        ax = ax1.axis()
        ax1.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])
        ax = ax2.axis()
        ax2.axis([ax[0], opts.freq_xmax, ax[2], ax[3]])

      plt.savefig(figname)

    plt.close(fig)

##################################################
#
#      multi-gen plot freq-domain plots
#
##################################################
  if opts.multi_gen_amp_freq:
    if opts.verbose: print "\tmulti_gen_amp_freq"
    fig, axs = nm_p.multi_gen_amp_plot(freq, fq, gens, n_l_m=n_l_m, mode_nums=opts.mode_nums)

    for ax in axs:
      ax.set_ylabel(r'|$\widetilde{'+opts.fc+'}_i|$')
      if not opts.lin_amp_freq:
        ax.set_yscale('log')

      if not opts.no_legend:
        ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

      if opts.freq_ymin:
        _ax = ax.axis()
        ax.axis([_ax[0], _ax[1], opts.freq_ymin, _ax[3]])
      if opts.freq_ymax:
        _ax = ax.axis()
        ax.axis([_ax[0], _ax[1], _ax[2], opts.freq_ymax])
      if opts.freq_xmin:
        _ax = ax.axis()
        ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
      if opts.freq_xmax:
        _ax = ax.axis()
        ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

      # grid lines
      ax.grid(opts.grid)

    ax.set_xlabel(r'$f/f_{\mathrm{orb}}$')

    if opts.lin_freq:
      if opts.lin_amp_freq:
        figname = opts.freqfilename + ".multi-gen_Lamp.lin_freq" + opts.tag + ".png"
      else:
        figname = opts.freqfilename + ".multi-gen_amp.lin_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      plt.savefig(figname)

    if opts.log_freq:
      if opts.lin_amp_freq:
        figname = opts.freqfilename + ".multi-gen_Lamp.log_freq" + opts.tag + ".png"
      else:
        figname = opts.freqfilename + ".multi-gen_amp.log_freq" + opts.tag + ".png"
      if opts.verbose: print "saving "+figname
      for ax in axs:
        ax.set_xscale('log')
        if opts.freq_xmin:
          _ax = ax.axis()
          ax.axis([opts.freq_xmin, _ax[1], _ax[2], _ax[3]])
        if opts.freq_xmax:
          _ax = ax.axis()
          ax.axis([_ax[0], opts.freq_xmax, _ax[2], _ax[3]])

      plt.savefig(figname)

    plt.close(fig)

