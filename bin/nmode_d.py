#!python_alias
usage="""
diagnostic.py [--options]

a module written to help interpret network structure and properities
"""

thoughts="""
X  make scatter plots of mode parameters

   make coupling/network visualization via an adjacency bundled tree structure

   compute probability distributions for different dynamical quantities: E_i, Hns, etc.

   correlation matricies between different modes?

   cross correlations between modes -> look for non-random behavior

X   we can generate scatter plots for 3mode detunings, 3mode Ethrs, heuristics for each coupling
      ==> these can be histogramed and the natural break down is according to generation (associate all couplings in which the parent mode is of a generation with that generation)
"""

import sys, os

import numpy as np

import nmode_utils as nmu
import nmode_state as nms
import nmode_diagnostic as nmd
#from nmode_plotting import E_distrib, H_distrib, multi_gen_E_distrib
  
from optparse import *

####################################################################################################

parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-F", "--outfilename", default=False, type="string")
parser.add_option("", "--tmin", default=False, type="float")
parser.add_option("", "--tmax", default=False, type="float")
parser.add_option("", "--downsample", default=False, type="int")
parser.add_option("", "--tcurrent", default=False, type="string")

parser.add_option("-l", "--logfilename", default=False, type="string")

### stacked histogram options
parser.add_option("", "--stacked-hist", default=False, type="string", help="generate a stacked histogram. input should be a list of parameters over which we make histograms. eg: \"w l\"")
parser.add_option("", "--multi-gen-stacked-hist", default=False, action="store_true", help="generates a stacked histogram with curves stacked by generation. MUST BE USED WITH --stacked-hist")
parser.add_option("", "--log-stacked-hist", default=False, action="store_true", help="set stacked histogram y scale to logarithmic. MUST BE USED WITH --stacked-hist")
parser.add_option("", "--stacked-hist-bin-width", default=False, type="string", help="bin widths for stacked histogram. input a string with the same number of entries as stacked hist")

### conc diagrams
parser.add_option("", "--Ei-conc", default=False, action="store_true", help="generate a \"conc\" diagram that maps the fraction of network energy (all energy associated with individual modes) as a function of the number of modes")
parser.add_option("", "--disp-conc", default=False, action="store_true", help="generate a \"conc\" diagram that maps the fraction of network dissipation as a function of the number of modes")

#parser.add_option("", "--multi-gen-Ei-conc", default=False, action="store_true", help="generate a \"conc\" diagram split by generation that maps the fraction of network energy (all energy associated with individual modes) as a function of the number of modes")
#parser.add_option("", "--multi-gen-disp-conc", default=False, action="store_true", help="generate a \"conc\" diagram split by generation that maps the fraction of network dissipation as a function of the number of modes")

parser.add_option("", "--annotate-conc", default=1, type="int", help="the number of leading modes to annotate on conc diagrams.")
parser.add_option("", "--annotate-fontsize", default=8, type="float")
parser.add_option("", "--conc-ymin", default=False, type="float", help="put a lower bound on conentration diagrams")
parser.add_option("", "--Cconc-ymin", default=False, type="float", help="put a lowe bound on cumulative concentration diagrams")

parser.add_option("", "--conc-fit", default=False, action="store_true", help="fit a broken power law to concentration diagrams")
parser.add_option("", "--conc-fit-start", default=0, type="int")
parser.add_option("", "--conc-fit-stop", default=np.infty, type="int") 
parser.add_option("", "--conc-fit-iters", default=100, type="int")
parser.add_option("", "--conc-fit-rtol", default=1e-8, type="float")

parser.add_option("", "--conc-verbose", default=False, action="store_true")

### coupling diagram options
parser.add_option("", "--coupling-diagram", default=False, type="string", help="the type of coupling diagram to generate. Supply at least one of \"w\", \"nlm\", \"nl-placement\", \"nl-sibling\", \"nl-shared_parent\", \"nl-shared_child\", \"nl-triples\"")
parser.add_option("", "--coupling-diagram-coloration", default="", type="string", help="coloration schema for coupling diagrams that support it (nl)")
parser.add_option("", "--coupling-diagram-colormap", default="copper", type="string", help="colormap for mode_colors in coupling diagrams that support it")
parser.add_option("", "--coupling-diagram-genNos", default="", type="string", help="a string of the genNos to include in \"nl-*\" coupling diagrams")
parser.add_option("", "--coupling-diagram-modeNo", default=-1, type=int, help="a mode number that will be highlighted in the coupling diagram.")

parser.add_option("", "--coupling-hist", default=False, action="store_true")
parser.add_option("", "--coupling-hist-num-bins", default=50, type="int")
parser.add_option("", "--coupling-hist-log", default=False, action="store_true")
parser.add_option("", "--coupling-hist-genNos", default=False, type="string")

### scatter plots of triples parameters
parser.add_option("", "--scat-detuning", default=False, action="store_true", help="make a scatter plot of energy in a coupling vs the three-mode detuning")
parser.add_option("", "--scat-heuristic", default=False, action="store_true", help="make a scatter plot of energy in a coupling vs the ranking heuristic")
parser.add_option("", "--scat-Ethr", default=False, action="store_true", help="make a scatter plot of energy in a coupling vs the exact 3mode threshold energy")

### distributions of mode energy, dissipation, etc
parser.add_option("", "--Ei-distrib", action="store_true")
parser.add_option("", "--Hns-distrib", action="store_true")
parser.add_option("", "--multi-gen-Ei-distrib", action="store_true")

parser.add_option("", "--num-bins", default=50, type="int", help="the number of bins used for *-distrib plots")
parser.add_option("", "--log-bins", default=False, action="store_true", help="logarithmically space the bins for *-distrib plots")
parser.add_option("", "--minE", default=1e-40, type="float", help="the minimum energy used for *-distrib plots")
parser.add_option("", "--maxE", default=1e-10, type="float", help="the maximum energy used for *-distrib plots")

parser.add_option("", "--no-legend", action="store_true", default=False, help="do not plot a legend on any figure")
parser.add_option("", "--legend-col", default=1, type="int", help="the number of columns for your legend")
parser.add_option("", "--legend-font-size", default = 12, type="int", help="the font size for text in the legend.")
parser.add_option("", "--legend-loc", default="best", type="string", help="location for legend.")

parser.add_option("-g", "--grid", default=False, action="store_true")
parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

###

if opts.tag != "":
  opts.tag = "." + opts.tag

if opts.coupling_hist_genNos:
  opts.coupling_hist_genNos = [int(l) for l in opts.coupling_hist_genNos.split()]

opts.coupling_diagram_coloration = opts.coupling_diagram_coloration.split()


## summary logicals
scat3 = (opts.scat_detuning or opts.scat_heuristic or opts.scat_Ethr)
distrib = (opts.Ei_distrib or opts.Hns_distrib or opts.multi_gen_Ei_distrib)
multi_gen = (opts.multi_gen_Ei_distrib or opts.multi_gen_stacked_hist or opts.coupling_hist) #or opts.multi_gen_Ei_conc or opts.multi_gen_disp_conc)
conc = (opts.Ei_conc or opts.disp_conc) #or opts.multi_gen_Ei_conc or opts.multi_gen_disp_conc)

if (opts.stacked_hist or scat3 or distrib or conc or (opts.coupling_diagram and opts.coupling_diagram_coloration)) and (not opts.outfilename):
  opts.outfilename = raw_input("outfilename = ")

if (opts.outfilename) and (not opts.tcurrent):
  opts.tcurrent = raw_input("current time variable type = ")

if (opts.stacked_hist or opts.coupling_diagram or scat3 or distrib or multi_gen or conc) and (not opts.logfilename):
  opts.logfilename = raw_input("logfilename = ")

if opts.stacked_hist:
  opts.stacked_hist = opts.stacked_hist.split()

  if opts.stacked_hist_bin_width:
    stacked_hist_bin_width = [float(l) for l in opts.stacked_hist_bin_width.split()]
    opts.stacked_hist = [(opts.stacked_hist[l], stacked_hist_bin_width[l]) for l in range(len(opts.stacked_hist))]
  else:
    opts.stacked_hist = [(l, False) for l in opts.stacked_hist]

### i'm too lazy to fix this right now
if opts.multi_gen_stacked_hist:
  print "\nWARNING: multi_gen variability has NOT been implemented for multi_gen_stacked_hist. INSTEAD, variability measure are for the entire network rather than individual generations\n"

####################################################################################################
#
#
#       load data
#
#
####################################################################################################

if opts.outfilename:
  if opts.verbose: print "loading integration data from %s" % opts.outfilename
  t_P, q, N_m = nmu.load_out(opts.outfilename, tmin=opts.tmin, tmax=opts.tmax, downsample=opts.downsample)

if opts.logfilename:
  if opts.verbose: print "loading system parameters from %s" % opts.logfilename
  system = nmu.load_log(opts.logfilename)
  network = system.network

  if conc:
    opts.conc_fit_stop = min(opts.conc_fit_stop, len(network))

if multi_gen:
  if opts.verbose: print "determining generational structure"
  gens, _ = system.network.gens()

####################################################################################################
#
#
#           stacked histograms
#
#
####################################################################################################
if opts.stacked_hist:
  if opts.verbose: print "\tstacked_hist"

  ### prepare data
  E = nms.compute_E(q, Eo=1.)
  mE = [nms.sample_mean(e) for e in E]
#  sE = [nms.sample_var(e, xo=mE[ind])**0.5 for ind, e in enumerate(E)]

  mdE = [2*network.wyU[ind][1]*e for ind, e in enumerate(mE)]
#  sdE = [2*network.wyU[ind][1]*e for ind, e in enumerate(sE)]

  data = [(1, len(network.K[ind]), mE[ind], mdE[ind]) for ind in range(N_m)] # number of modes, number of connections, mean Energy, mean dissipation

  ### make and save each plot
  for xvar, bin_width in opts.stacked_hist:

    if xvar == "w":
      if opts.verbose: print "\t\tw"
      if not bin_width:
        bin_width = 0.05
      bins = np.arange(-2.5-bin_width/2, 2.5+bin_width/2+bin_width, bin_width)

      if opts.multi_gen_stacked_hist:
        fig, axs = nmd.generational_stacked_histogram(network, bins, [w/system.Oorb for w,y,U in network.wyU], data, log=opts.log_stacked_hist)
      else:
        fig, axs = nmd.stacked_histogram(bins, [w/system.Oorb for w,y,U in network.wyU], data, log=opts.log_stacked_hist)
      axs[-1].set_xlabel(r"$\omega/\Omega_{\mathrm{orb}}$")

      # variability measures
      ax2_ymin = axs[2].get_ylim()[0]
      ax3_ymin = axs[3].get_ylim()[0]

      binned_modeNos = nmd.bin_by_w(bins*system.Oorb, network)
      for binNo, bin in enumerate(binned_modeNos):
        me = sum([mE[i] for i in bin])
        se = nms.sample_var([ sum([E[i][ind] for i in bin]) for ind in range(len(E[0])) ], xo=me)**0.5

        mde = sum([mdE[i] for i in bin])
        sde = nms.sample_var([ sum([2*network.wyU[i][1]*E[i][ind] for i in bin]) for ind in range(len(E[0]))], xo=mde)**0.5

        if opts.log_stacked_hist:
          axs[2].fill_between(bins[binNo:binNo+2], (me+se)*np.ones((2,)), max((me-se), ax2_ymin)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)
          axs[3].fill_between(bins[binNo:binNo+2], (mde+sde)*np.ones((2,)), max((mde-sde), ax3_ymin)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)
        else:
          axs[2].fill_between(bins[binNo:binNo+2], (me+se)*np.ones((2,)), (me-se)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)
          axs[3].fill_between(bins[binNo:binNo+2], (mde+sde)*np.ones((2,)), (mde-sde)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)

      axs[2].set_ylim(ymin=ax2_ymin)
      axs[3].set_ylim(ymin=ax3_ymin)

      for ax in axs:
        ax.set_xlim(xmin=-2.5, xmax=2.5)
        ax.grid(opts.grid, which="both")

    if xvar == "l":
      if opts.verbose: print "\t\tl"
      if not bin_width:
        bin_width = 1.0
      ls = [l for n,l,m in network.nlm]
      bins = np.arange(0.5, max(ls)+0.5+bin_width, bin_width)
      if opts.multi_gen_stacked_hist:
        fig, axs = nmd.generational_stacked_histogram(network, bins, ls, data, log=opts.log_stacked_hist)
      else:
        fig, axs = nmd.stacked_histogram(bins, ls, data, log=opts.log_stacked_hist)
      axs[-1].set_xlabel(r"l")

      # variability measures
      ax2_ymin = axs[2].get_ylim()[0]
      ax3_ymin = axs[3].get_ylim()[0]

      binned_modeNos = nmd.bin_by_l(bins*system.Oorb, network)
      for binNo, bin in enumerate(binned_modeNos):
        me = sum([mE[i] for i in bin])
        se = nms.sample_var([ sum([E[i][ind] for i in bin]) for ind in range(len(E[0])) ], xo=me)**0.5

        mde = sum([mdE[i] for i in bin])
        sde = nms.sample_var([ sum([2*network.wyU[i][1]*E[i][ind] for i in bin]) for ind in range(len(E[0]))], xo=mde)**0.5

        if opts.log_stacked_hist:
          axs[2].fill_between(bins[binNo:binNo+2], (me+se)*np.ones((2,)), max((me-se), ax2_ymin)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)
          axs[3].fill_between(bins[binNo:binNo+2], (mde+sde)*np.ones((2,)), max((mde-sde), ax3_ymin)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)
        else:
          axs[2].fill_between(bins[binNo:binNo+2], (me+se)*np.ones((2,)), (me-se)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)
          axs[3].fill_between(bins[binNo:binNo+2], (mde+sde)*np.ones((2,)), (mde-sde)*np.ones((2,)), facecolor='b', edgecolor="none", alpha=0.2)

      axs[2].set_ylim(ymin=ax2_ymin)
      axs[3].set_ylim(ymin=ax3_ymin)

      for ax in axs:
        ax.set_xlim(xmin=0, xmax=np.ceil(bins[-1]))
        ax.grid(opts.grid, which="both")

    axs[0].set_ylabel(r"No. modes")
    axs[1].set_ylabel(r"No. $\kappa$")
    axs[2].set_ylabel(r"$A_i^2$")
    axs[3].set_ylabel(r"$2 \gamma_i A_i^2$")

    ### save
    insert = "stacked_hist."
    if opts.log_stacked_hist: insert = "log_"+insert
    if opts.multi_gen_stacked_hist: insert = "multi_gen_"+insert
    figname = opts.outfilename + "." + insert + xvar + opts.tag + ".png"
    if opts.verbose: print "saving "+figname
    fig.savefig(figname)
    nmd.plt.close(fig)

####################################################################################################
#
#
#                             network structure maps
#
#
####################################################################################################
if opts.coupling_hist:
  if opts.verbose: print "coupling_hist"

  fig, ax = nmd.coupling_hist(network, gens=gens, num_bins=opts.coupling_hist_num_bins, log=opts.coupling_hist_log, genNos=opts.coupling_hist_genNos)

  ax.grid(opts.grid)

  if not opts.no_legend:
    ax.legend(loc=opts.legend_loc, ncol=opts.legend_col, prop={'size':opts.legend_font_size})

  figname = opts.logfilename+".coupling_hist"+opts.tag+".png"
  if opts.verbose: print "saving "+figname
  fig.savefig(figname)
  nmd.plt.close(fig)

##################################################
if opts.coupling_diagram:

  for diagram_type in opts.coupling_diagram.split():
    if diagram_type == "w":
      if opts.verbose: print "\tbuilding coupling_diagram \"w\""
      fig, ax, gens, coups, edges = nmd.coupling_tree_w(system, verbose=opts.verbose) # delegate

      if opts.coupling_diagram_modeNo >= 0:
        if opts.verbose: print "\thighlighting modeNo: %d" % opts.coupling_diagram_modeNo
        nmd.couling_tree_highlight(opts.coupling_diagram_modeNo, gens, coups, edges)
 
      ax.grid(opts.grid)

      figname = opts.logfilename+".coupling_diagram-"+diagram_type+opts.tag+".png"
      if opts.verbose: print "saving "+figname
      fig.savefig(figname)
      nmd.plt.close(fig)

    elif diagram_type == "nlm":
      if opts.verbose: print "\tbuilding coupling_diagram \"nlm\""
      fig, ax, gens, coups, edges = nmd.coupling_tree_nlm(system, verbose=opts.verbose) # delegate

      if opts.coupling_diagram_modeNo >= 0:
        if opts.verbose: print "\thighlighting modeNo: %d" % opts.coupling_diagram_modeNo
        nmd.couling_tree_highlight(opts.coupling_diagram_modeNo, gens, coups, edges)
 
      ax.grid(opts.grid)

      figname = opts.logfilename+".coupling_diagram-"+diagram_type+opts.tag+".png"
      if opts.verbose: print "saving "+figname
      fig.savefig(figname)
      nmd.plt.close(fig)

    elif "nl-" in diagram_type: ### a class of diagrams
      if opts.verbose: print "\t building coupling_diagaram \"%s\"" % diagram_type
      if opts.coupling_diagram_coloration:

        ### define colormap
        colormap = nmd.plt.get_cmap(opts.coupling_diagram_colormap)

        ### compute colors for each mode based on coloration scheme
        for coloration in opts.coupling_diagram_coloration:
          if coloration == "A":
            if opts.verbose: print "\t\t coloration: A"
            mA = np.array([np.mean(_) for _ in nms.compute_A(q, Eo=1.0)])
            mode_colors = [colormap(_) for _ in mA/max(mA)]
            mode_order = mA.argsort() # smallest values plotted first

          elif coloration == "E":
            if opts.verbose: print "\t\t coloration: E"
            mE = np.array([np.mean(_) for _ in nms.compute_E(q, Eo=1.0)])
            mode_colors = [colormap(_) for _ in mE/max(mE)]
            mode_order = mE.argsort()

          elif coloration == "disp":
            if opts.verbose: print "\t\t coloration: disp"
            myE = -np.array([np.mean(_) for _ in nms.viscous_disp(q, network, Eo=1.0)[-1]])
            mode_colors = [colormap(_) for _ in myE/max(myE)]
            mode_order = myE.argsort()

          else:
            if opts.verbose: print "\t\tcoloration \"%s\" not understood. skipping..." % coloration
            continue

          fig = nmd.plt.figure()
          ax = fig.add_axes([0.1, 0.1, 0.750, 0.8])
          cbax = fig.add_axes([0.875, 0.1, 0.025, 0.8])

          fig, ax = nmd.coupling_tree_nl(system, tree_type=diagram_type.split("-")[-1], genNos=[int(l) for l in opts.coupling_diagram_genNos.split()], verbose=opts.verbose, mode_colors=mode_colors, mode_order=mode_order, fig_ax=(fig,ax))
          ax.grid(opts.grid)

          ### add color bar!
#          colorbar = fig.colorbar(cax=cbax, orientation="vertical")
          colorbar = nmd.matplotlib.colorbar.ColorbarBase(cbax, cmap=colormap, orientation='vertical')
          if coloration == "A":
            colorbar.ax.set_ylabel(r"$A_i/\mathrm{max}\{A_i\}$")
          elif coloration == "E":
            colorbar.ax.set_ylabel(r"$A_i^2/\mathrm{max}\{A_i^2\}$")
          elif coloration == "disp":
            colorbar.ax.set_ylabel(r"$\gamma_i A_i^2/\mathrm{max}\{\gamma_i A_i^2\}$")

          ax.grid(opts.grid)

          figname = opts.logfilename+".coupling_diagram-"+diagram_type+"-"+coloration+opts.tag+".png"
          if opts.verbose: print "saving "+figname
          fig.savefig(figname)
          nmd.plt.close(fig)

      else:
        fig, ax = nmd.coupling_tree_nl(system, tree_type=diagram_type.split("-")[-1], genNos=[int(l) for l in opts.coupling_diagram_genNos.split()], verbose=opts.verbose)
        ax.grid(opts.grid)
        figname = opts.logfilename+".coupling_diagram-"+diagram_type+opts.tag+".png"
        if opts.verbose: print "saving "+figname
        fig.savefig(figname)
        nmd.plt.close(fig)

    else:
      if opts.verbose: print "\tcoupling_diagram \"%s\" not understood. skipping..." % diagram_type

####################################################################################################
#
#
#                             energies associated with each mode
#
#
####################################################################################################

### mode energy prob distributions
if opts.Ei_distrib:
  if opts.verbose: print "\tEi_distrib"
  fig, ax = nmd.E_distrib(q, minE=opts.minE, maxE=opts.maxE, num_bins=opts.num_bins, log=False, log_bins=opts.log_bins, n_l_m=network.nlm, mode_nums=False)

  ax.set_xlabel(r"$A_i^2$")
  ax.set_ylabel(r"$p(A_i^2)$")

  ax.grid(opts.grid) # grid lines

  figname = opts.outfilename+".Ei_distrib"+opts.tag+".png"
  if opts.verbose: print "saving "+figname
  nmd.plt.savefig(figname)
  nmd.plt.close(fig)


### mulit-gen energy prob distribution
if opts.multi_gen_Ei_distrib:
  if opts.verbose: print "\tmulti_gen_Ei_distrib"
  fig, axs = nmd.multi_gen_E_distrib(q, gens, minE=opts.minE, maxE=opts.maxE, num_bins=opts.num_bins, log=False, log_bins=opts.log_bins, n_l_m=network.nlm, mode_nums=False)
  
  for ax in axs:
    ax.set_ylabel(r'$p(A_i^2)$')
    ax.grid(opts.grid) # grid lines

  ax.set_xlabel(r'$A_i^2$')

  figname = opts.outfilename+".multi-gen_Ei_distrib"+opts.tag+".png"
  if opts.verbose: print "saving "+figname
  nmd.plt.savefig(figname)
  nmd.plt.close(fig)

####################################################################################################
#
#
#                             energies associated with each coupling
#
#
####################################################################################################
"""
probability distributions (sample over time). 

        WILL THIS BE USEFUL????
"""

### scatter plots
if scat3:
  if opts.verbose: print "computing Hns_coup"
  if opts.tcurrent == "q":
    Hns_coup = nms.Hns_coup_q(q, network, coup=False, Eo=1.0)
  elif opts.tcurrent == "x":
    Hns_coup = nms.Hns_coup_x(t_P, q, system, coup=False, Eo=1.0)

  ### detuning
  if opts.scat_detuning:
    if opts.verbose: print "\tscat_detuning"
    fig, ax, ax_Hns, ax_det = nmd.Hns_coup_detuning(Hns_coup, system)

    ax.grid(opts.grid)
    ax_Hns.grid(opts.grid)
    ax_det.grid(opts.grid)

    figname = opts.outfilename+".scat_detuning"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    fig.savefig(figname)
    nmd.plt.close(fig)

  ### heuristic
  if opts.scat_heuristic:
    if opts.verbose: print "\tscat_heuristic"
    fig, ax, ax_Hns, ax_heuristic = nmd.Hns_coup_heuristic(Hns_coup, system)

    ax.grid(opts.grid)
    ax_Hns.grid(opts.grid)
    ax_heuristic.grid(opts.grid)

    figname = opts.outfilename+".scat_heuristic"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    fig.savefig(figname)
    nmd.plt.close(fig)

  ### Ethr
  if opts.scat_Ethr:
    if opts.verbose: print "\tscat_Ethr"
    fig, ax, ax_Hns, ax_Ethr = nmd.Hns_coup_Ethr(Hns_coup, system)

    ax.grid(opts.grid)
    ax_Hns.grid(opts.grid)
    ax_Ethr.grid(opts.grid)

    figname = opts.outfilename+".scat_Ethr"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    fig.savefig(figname)
    nmd.plt.close(fig)

####################################################################################################
#
#
#                             energies associated with entire network
#
#
####################################################################################################


if opts.Hns_distrib:

  if opts.verbose: print "\tHns_distrib"
  if opts.tcurrent == "q":
    fig, ax = nmd.H_distrib(nms.Hns_q(q, network, Eo=1.) , minE=opts.minE, maxE=opts.maxE, num_bins=opts.num_bins, log=False, log_bins=opts.log_bins)
  elif opts.tcurrent=="x":
    fig, ax = nmd.H_distrib(nms.Hns_x(t_P, q, system, Eo=1.) , minE=opts.minE, maxE=opts.maxE, num_bins=opts.num_bins, log=False, log_bins=opts.log_bins)


  ax.set_xlabel(r"$H_{*}$")
  ax.set_ylabel(r"$p(H_{*})$")

  ax.grid(opts.grid)

  figname = opts.outfilename+".Hns_distrib"+opts.tag+".png"
  if opts.verbose: print "saving "+figname
  nmd.plt.savefig(figname)
  nmd.plt.close(fig)


####################################################################################################
#
#
#                                      conc diagrams
#
#
####################################################################################################
if opts.Ei_conc:# or opts.multi_gen_Ei_conc:

#  E = nms.compute_E(q, Eo=1.)
#  mE = [nms.sample_mean(x) for x in E]
#  sE = [nms.sample_var(x, xo=mx)**0.5 for x, mx in zip(E, mE)]
  mE = [nms.sample_mean(x) for x in nms.compute_E(q,Eo=1.)]

  fitparams=False
  if opts.conc_fit:
    if opts.verbose: print "computing fit for Ei_conc"
#    data = zip(mE, sE)
    data = zip(mE, [1.0 for l in mE])
    data.sort(key=lambda l:l[0], reverse=True)

    fitparams, covar, red_chi2 =  nms.broken_PowLaw_fitter(range(opts.conc_fit_start+1, opts.conc_fit_stop+1), np.array([l[0] for l in data[opts.conc_fit_start:opts.conc_fit_stop]])/sum(mE), [l[1] for l in data[opts.conc_fit_start:opts.conc_fit_stop]], max_iters=opts.conc_fit_iters, rtol=opts.conc_fit_rtol, verbose=opts.conc_verbose) # set errors on mE estimates to 1.0 so they don't affect fit
    fitparams = list(fitparams)

  ### Ei_conc
  if opts.Ei_conc:
    if opts.verbose: print "\tEi_conc"
    ##############################################
    ### cumulative
    fig, [ax_ul, ax_ll, ax_lr] = nmd.Cconc(mE, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, fitparams=fitparams, fitrange=[opts.conc_fit_start, opts.conc_fit_stop])

    ax_ll.set_xlabel("No. of modes")
    ax_lr.set_xlabel("No. of modes")
    ax_ul.set_ylabel(r"fraction of $\sum_i A_i^2$")
    ax_ll.set_ylabel(r"fraction of $\sum_i A_i^2$")

    if opts.Cconc_ymin:
      ax_ul.set_ylim(ymin=opts.Cconc_ymin)
      ax_ll.set_ylim(ymin=opts.Cconc_ymin)
      ax_lr.set_ylim(ymin=opts.Cconc_ymin)

    ax_ul.grid(opts.grid)
    ax_ll.grid(opts.grid)
    ax_lr.grid(opts.grid)

    if opts.conc_fit:
      annotation = r"""$ A_n^2 = C \cdot n^{-\alpha} \cdot \left(\frac{1+ e^{n/\beta}}{2}\right)^{-\gamma}$

$C = %.6f \pm %.6f$
$\alpha  = %.6f \pm %.6f$
$\beta = %.6f \pm %.6f$
$\gamma = %.6f \pm %.6f$

$\chi^2_{\mathrm{reduced}} = %.6f$""" % (fitparams[0], covar[0]**0.5, fitparams[1], covar[1]**0.5, fitparams[2], covar[2]**0.5, fitparams[3], covar[3]**0.5, red_chi2) 
      nmd.plt.figtext(0.90, 0.75, annotation , ha='right', va='center', color='r')

    figname = opts.outfilename+".Ei_Cconc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)
    nmd.plt.close(fig)

    ##############################################
    ### differential
    fig, [ax_ul, ax_ll, ax_lr] = nmd.conc(mE, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, fitparams=fitparams, fitrange=[opts.conc_fit_start, opts.conc_fit_stop])

    ax_ll.set_xlabel("No. of modes")
    ax_lr.set_xlabel("No. of modes")
    ax_ll.set_ylabel(r"fraction of $\sum_i A_i^2$")
    ax_ul.set_ylabel(r"fraction of $\sum_i A_i^2$")

    if opts.conc_ymin:
      ax_ul.set_ylim(ymin=opts.conc_ymin)
      ax_ll.set_ylim(ymin=opts.conc_ymin)
      ax_lr.set_ylim(ymin=opts.conc_ymin)

    ax_ul.grid(opts.grid)
    ax_ll.grid(opts.grid)
    ax_lr.grid(opts.grid)

    if opts.conc_fit:
      annotation = r"""$ A_n^2 = C \cdot n^{-\alpha} \cdot \left(\frac{1+ e^{n/\beta}}{2}\right)^{-\gamma}$

$C = %.6f \pm %.6f$
$\alpha  = %.6f \pm %.6f$
$\beta = %.6f \pm %.6f$
$\gamma = %.6f \pm %.6f$

$\chi^2_{\mathrm{reduced}} = %.6f$""" % (fitparams[0], covar[0]**0.5, fitparams[1], covar[1]**0.5, fitparams[2], covar[2]**0.5, fitparams[3], covar[3]**0.5, red_chi2)
      nmd.plt.figtext(0.90, 0.75, annotation , ha='right', va='center', color='r')

    figname = opts.outfilename+".Ei_conc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)
    nmd.plt.close(fig)

  '''
  ### multi_gen_Ei_conc
  if opts.multi_gen_Ei_conc:
    if opts.verbose: print "\tmulti_gen_Ei_conc"
    ### cumulative
    fig, axs = nmd.multi_gen_Cconc(mE, gens, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, log=opts.log_conc)

    for genNo, ax in enumerate(axs):
      ax.set_ylabel('fraction of\n'+r'$\sum_{i\in'+str(genNo)+'} A_i^2$')
      ax.grid(opts.grid) # grid lines

      if opts.conc_ymin:
        ax.set_ylim(ymin=opts.conc_ymin)

    ax.set_xlabel(r'fraction of modes')

    if opts.log_conc:
      figname = opts.outfilename+".multi-gen_Ei_LCconc"+opts.tag+".png"
    else:
      figname = opts.outfilename+".multi-gen_Ei_Cconc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)

    nmd.plt.close(fig)

    ### differential
    fig, axs = nmd.multi_gen_conc(mE, gens, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, log=opts.log_conc)

    for genNo, ax in enumerate(axs):
      ax.set_ylabel('fraction of\n'+r'$\sum_{i\in'+str(genNo)+'} A_i^2$')
      ax.grid(opts.grid) # grid lines

      if opts.conc_ymin:
        ax.set_ylim(ymin=opts.conc_ymin)

    ax.set_xlabel(r'fraction of modes')

    if opts.log_conc:
      figname = opts.outfilename+".multi-gen_Ei_Lconc"+opts.tag+".png"
    else:
      figname = opts.outfilename+".multi-gen_Ei_conc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)

    nmd.plt.close(fig)
  '''

##################################################
if opts.disp_conc: #or opts.multi_gen_disp_conc:

#  yE = -1.0*np.array(nms.viscous_disp(q, network, Eo=1.)[1]) # make these positive...
#  myE = [nms.sample_mean(x) for x in yE]
#  syE = [nms.sample_var(x, xo=mx)**0.5 for x, mx in zip(yE, myE)]
  myE = [-1.0*nms.sample_mean(x) for x in nms.viscous_disp(q, network, Eo=1.)[1]]

  fitparams=False
  if opts.conc_fit:
    if opts.verbose: print "computing fit for disp_conc"
    data = zip(myE, [1.0 for l in myE]) # set errors on mE estimates to 1.0 so they don't affect fit
    data.sort(key=lambda l:l[0], reverse=True)

    fitparams, covar, red_chi2 =  nms.broken_PowLaw_fitter(range(opts.conc_fit_start+1, opts.conc_fit_stop+1), np.array([l[0] for l in data[opts.conc_fit_start:opts.conc_fit_stop]])/sum(myE), [l[1] for l in data[opts.conc_fit_start:opts.conc_fit_stop]], max_iters=opts.conc_fit_iters, rtol=opts.conc_fit_rtol, verbose=opts.conc_verbose)
    fitparams = list(fitparams)

  ### disp_conc
  if opts.disp_conc:
    if opts.verbose: print "\tdisp_conc"
    ##############################################
    ### cumulative
    fig, [ax_ul, ax_ll, ax_lr] = nmd.Cconc(myE, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, fitparams=fitparams, fitrange=[opts.conc_fit_start, opts.conc_fit_stop])

    ax_ll.set_xlabel("No. of modes")
    ax_lr.set_xlabel("No. of modes")
    ax_ul.set_ylabel(r"fraction of $2 \sum_i \gamma_i A_i^2$")
    ax_ll.set_ylabel(r"fraction of $2 \sum_i \gamma_i A_i^2$")

    if opts.Cconc_ymin:
      ax_ul.set_ylim(ymin=opts.Cconc_ymin)
      ax_ll.set_ylim(ymin=opts.Cconc_ymin)
      ax_lr.set_ylim(ymin=opts.Cconc_ymin)

    ax_ul.grid(opts.grid)
    ax_ll.grid(opts.grid)
    ax_lr.grid(opts.grid)

    if opts.conc_fit:
      annotation = r"""$2\gamma_{n} A_{n}^2 = C \cdot n^{-\alpha} \cdot \left(\frac{1+ e^{n/\beta}}{2}\right)^{-\gamma}$

$C = %.6f \pm %.6f$
$\alpha = %.6f \pm %.6f$ 
$\beta = %.6f \pm %.6f$
$\gamma = %.6f \pm %.6f$

$\chi^2_{\mathrm{reduced}} = %.6f$""" % (fitparams[0], covar[0]**0.5, fitparams[1], covar[1]**0.5, fitparams[2], covar[2]**0.5, fitparams[3], covar[3]**0.5, red_chi2)
      nmd.plt.figtext(0.90, 0.75, annotation, ha='right', va='center', color='r')

    figname = opts.outfilename+".disp_Cconc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)
    nmd.plt.close(fig)

    ##############################################
    ### differential
    fig, [ax_ul, ax_ll, ax_lr] = nmd.conc(myE, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, fitparams=fitparams, fitrange=[opts.conc_fit_start, opts.conc_fit_stop])

    ax_ll.set_xlabel("No. of modes")
    ax_lr.set_xlabel("No. of modes")
    ax_ll.set_ylabel(r"fraction of $2 \sum_i \gamma_i A_i^2$")
    ax_ul.set_ylabel(r"fraction of $2 \sum_i \gamma_i A_i^2$")

    if opts.conc_ymin:
      ax_ul.set_ylim(ymin=opts.conc_ymin)
      ax_ll.set_ylim(ymin=opts.conc_ymin)
      ax_lr.set_ylim(ymin=opts.conc_ymin)

    ax_ul.grid(opts.grid)
    ax_ll.grid(opts.grid)
    ax_lr.grid(opts.grid)

    if opts.conc_fit:
      annotation = r"""$2\gamma_{n} A_{n}^2 = C \cdot n^{-\alpha} \cdot \left(\frac{1+ e^{n/\beta}}{2}\right)^{-\gamma}$

$C = %.6f \pm %.6f$
$\alpha = %.6f \pm %.6f$ 
$\beta = %.6f \pm %.6f$
$\gamma = %.6f \pm %.6f$

$\chi^2_{\mathrm{reduced}} = %.6f$""" % (fitparams[0], covar[0]**0.5, fitparams[1], covar[1]**0.5, fitparams[2], covar[2]**0.5, fitparams[3], covar[3]**0.5, red_chi2)
      nmd.plt.figtext(0.90, 0.75, annotation , ha='right', va='center', color='r')

    figname = opts.outfilename+".disp_conc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)
    nmd.plt.close(fig)

  '''
  ### multi_gen_disp_conc
  if opts.multi_gen_disp_conc:
    if opts.verbose: print "\tmulti_gen_disp_conc"
    ### cumulative
    fig, axs = nmd.multi_gen_Cconc(myE, gens, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, log=opts.log_conc)      

    for genNo, ax in enumerate(axs):
      ax.set_ylabel('fraction of\n'+r'$2\sum_{i\in'+str(genNo)+'} \gamma_i A_i^2$')
      ax.grid(opts.grid) # grid lines

      if opts.conc_ymin:
        ax.set_ylim(ymin=opts.conc_ymin)

    ax.set_xlabel(r'fraction of modes')

    if opts.log_conc:
      figname = opts.outfilename+".multi-gen_disp_LCconc"+opts.tag+".png"
    else:
      figname = opts.outfilename+".multi-gen_disp_Cconc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)
    nmd.plt.close(fig)

    ### differential
    fig, axs = nmd.multi_gen_conc(myE, gens, network, annotate=opts.annotate_conc, annotate_fontsize=opts.annotate_fontsize, log=opts.log_conc)

    for genNo, ax in enumerate(axs):
      ax.set_ylabel('fraction of\n'+r'$2\sum_{i\in'+str(genNo)+'} \gamma_i A_i^2$')
      ax.grid(opts.grid) # grid lines

      if opts.conc_ymin:
        ax.set_ylim(ymin=opts.conc_ymin)

    ax.set_xlabel(r'fraction of modes')

    if opts.log_conc:
      figname = opts.outfilename+".multi-gen_disp_Lconc"+opts.tag+".png"
    else:
      figname = opts.outfilename+".multi-gen_disp_conc"+opts.tag+".png"
    if opts.verbose: print "saving "+figname
    nmd.plt.savefig(figname)
    nmd.plt.close(fig)
  '''
