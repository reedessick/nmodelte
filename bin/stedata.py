#!/usr/local/bin/python

usage = "written to generate a large table with all pertinant ste data in it"

from optparse import OptionParser
import nmode_utils as nmu
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

####################################################################################################

parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--cachefilename", default=False, type="string", help="the filename of a cachefile containing paths to all *.ste.pickle files of interest")
parser.add_option("-o", "--output-dir", default="./", type="string", help="the directory into which the resulting table will be written. DEFAULT=\"./\"")

parser.add_option("", "--D1-plots", default=False, action="store_true", help="attempts to build linear plots with annotations for all available statistics")
parser.add_option("", "--D2-plots", default=False, action="store_true", help="attempts to build plots of statistics vs. Nmodes, Ngens, Ntriples")

parser.add_option("", "--unit-system", default="CGS", type="string", help="system of units to use. \"SI\" or \"CGS\" (default)")

parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if opts.tag != "":
  opts.tag = "."+opts.tag

nmu.set_units(system=opts.unit_system)
units_time = nmu.units["time"]
units_energy = nmu.units["energy"]

####################################################################################################
#
# find and read in data
#
####################################################################################################
### find list of logs
if opts.cachefilename:
  if opts.verbose: print "reading in filenames from "+opts.cachefilename
  cachefile = open(opts.cachefilename, "r")
  filenames = sorted( [l.strip() for l in cachefile.readlines() ])
  cachefile.close()
else:
  import os
  import glob
  if opts.verbose: print "reading in all *ste.pickle filenames from "+os.getcwd()
  filenames = sorted( glob.glob("./*.ste.pickle") )

dat = []

### read data from logs
for filename in filenames:
  if opts.verbose: print "reading system from "+filename

  try:
    sdata, _ = nmu.load_ste(filename)
    dat.append( (filename, sdata) )
  except:
    print "\tERROR when reading "+filename
    continue

####################################################################################################
#
#  plot data
#
####################################################################################################
known_y  = [ 
("Hns" ,        r"$\left<H_{\ast}\right>$" ,                                                                                       "mean{Hns/|Eorb|}" ,                "stdv{Hns/|Eorb|}" ), 
("E" ,          r"$\left<\sum A_i^2\right> /E_{\mathrm{orb}}$" ,                                                                   "mean{sum{E}/|Eorb|}" ,             "stdv{sum{E}/|Eorb|}" ),
("yE" ,         r"$\left<2\sum \gamma_i A_i^2 \right> \cdot P_{\mathrm{orb}}/E_{\mathrm{orb}} $" ,                                 "mean{|sum{Edot}|*(Porb/|Eorb|)}",  "stdv{|sum{Edot}|*(Porb/|Eorb|)}" ),
("1-gini_E" ,   r"$1-\mathrm{GINI}\ \left<A_i^2\right>$" ,                                                                         "Gini_index{mean{sum{E}}}"     ,    False ), 
("1-gini_yE" ,  r"$1-\mathrm{GINI}\ \left<2\gamma_i A_i^2\right>$" ,                                                               "Gini_index{mean{sum{Edot}}}"  ,    False ),
("dEdt",        r"$\left< \mathrm{d}E/\mathrm{d}t \right> \cdot P_{\mathrm{orb}} / E_{\mathrm{orb}}$",                             "mean{sum{dE/dt*Porb}/|Eorb|}",     "stdv{sum{dE/dt*Porb}/|Eorb|}"),
#("1-gini_dEdt", r"$1-\mathrm{GINI}\ \left<\mathrm{d}E_i/\mathrm{d}t\right>$",                            "Gini_index{mean{sum{dE/dt}}}",     False),
("HintHns",     r"$\left<H_{int}+H_{\ast}\right>$",                                                                                "mean{Hint+Hns}/|Eorb|",            "stdv{Hint+Hns}/|Eorb|"),
("dHintHns",    r"$\left< \mathrm{d}\left(H_{int}+H_{\ast}\right)/\mathrm{d}t \right> \cdot P_{\mathrm{orb}} / E_{\mathrm{orb}}$", "mean{d(Hint+Hns)/dt}*Porb/|Eorb|", "stdv{d(Hint+Hns)/dt}*Porb/|Eorb|")
]
N_y = len(known_y)

miny = 1e-20

known_x  = [
("Nmodes", "$N_{\mathrm{modes}}$"),
("Ngens", "$N_{\mathrm{gens}}$"),
("Ntriples", "$N_{\mathrm{triples}}$")
]
N_x = len(known_x)

##################################################
if opts.D1_plots:
  if opts.verbose: print "1D_plots"
  axes_loc = [0.1,0.15,0.8,0.7]
  colors = ['b','r','g','c','m','y','k']
  
  for ystr, ylabel, ykey, yvkey in known_y:
    if opts.verbose: print "\t",ystr

    fig = plt.figure()
    ax = fig.add_axes(axes_loc)
    NgensD = {}
    for stefilename, sdata in dat:
      try:
        y = sdata["stats"][ykey]
      except KeyError:
        print "could not find data for %s in %s... skipping" % (ykey, stefilename)
        continue
      Ngens = sdata["system"]["Ngens"]

      if "Gini_index" in ykey:
        y = 1-np.array(y)

      if not NgensD.has_key(Ngens):
        ax.plot([y,y], [0,1], color=colors[Ngens], alpha=0.75, label="$N_{\mathrm{gens}}=%d$"%Ngens)
        NgensD[Ngens] = 1
      else:
        ax.plot([y,y], [0,1], color=colors[Ngens], alpha=0.75)

      if yvkey:
        yv = sdata["stats"][yvkey]
        ax.fill_between([max(y-yv, miny),y+yv], [0.25,0.25], [0.75,0.75], facecolor=colors[Ngens], edgecolor="none", alpha=0.05)

    ax.set_xlabel(ylabel)
    ax.set_xscale('log')
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.set_yticks([])
    ax.set_ylim(ymin=-0.1, ymax=1.1)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower center", ncol=len(NgensD.keys()), frameon=False)

    plt.setp(fig, figwidth=10, figheight=5)

    figname = "%s/%s%s.png" % (opts.output_dir, ystr, opts.tag)
    if opts.verbose: print "\t\tsaving ", figname
    fig.savefig(figname)
    plt.close(fig)

##################################################
if opts.D2_plots:
  if opts.verbose: print "2D_plots"
  axes_loc = [0.15,0.15,0.8,0.7]
  colors = ['b','r','g','c','m','y','k']

  for ystr, ylabel, ykey, yvkey in known_y:
    if opts.verbose: print "\t",ystr
    for xkey, xlabel in known_x:
      if opts.verbose: print "\t\t",xkey      

      fig = plt.figure()
      ax = fig.add_axes(axes_loc)
      NgensD = {}
      for stefilename, sdata in dat:
        try:
          y = sdata["stats"][ykey]
        except KeyError:
          print "could not find data for %s in %s... skipping" % (ykey, stefilename)
          continue

        if "Gini_index" in ykey:
          y = 1-np.array(y)

        Ngens = sdata["system"]["Ngens"]
        x = sdata["system"][xkey]

        if not NgensD.has_key(Ngens):
          ax.plot(x, y, color=colors[Ngens], markeredgecolor=colors[Ngens], markerfacecolor="none", marker="o", alpha=0.75, label="$N_{\mathrm{gens}}=%d$"%Ngens)
          NgensD[Ngens] = 1
        else:
          ax.plot(x, y, color=colors[Ngens], markeredgecolor=colors[Ngens], markerfacecolor="none", marker="o", alpha=0.75)

        if yvkey:
          yv = sdata["stats"][yvkey]
          ax.plot([x,x],[max(y-yv,miny),y+yv], color=colors[Ngens], alpha=0.25)

      ax.set_xlabel(xlabel)
      if xkey != "Ngens":
        ax.set_xscale('log')
      else: 
        maxNgens = max(NgensD.keys())
        ax.set_xticks(range(1,maxNgens+1))
        ax.set_xlim(xmin=0.75, xmax=maxNgens+0.25)
      ax.set_ylabel(ylabel)
      ax.set_yscale('log')
      ax.grid(True)
      ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower center", numpoints=1, frameon=False, ncol=max(NgensD.keys()))

      plt.setp(fig, figwidth=10, figheight=8)

      figname = "%s/%s-%s%s.png" % (opts.output_dir, xkey, ystr, opts.tag)
      if opts.verbose: print "\t\t\tsaving ", figname
      fig.savefig(figname)
      plt.close(fig)

####################################################################################################
#
#  write table >> file
#
####################################################################################################
### document strings
ldoc_preamble=r"""\documentclass[10pt]{article}
\usepackage{fullpage}
\usepackage{multirow}
\usepackage[landscape, top=0.1cm, bottom=0.1cm, left=0.1cm, right=0.1cm]{geometry}
\usepackage{floatrow}
\DeclareFloatFont{fontsize}{\small}
\floatsetup[table]{font=fontsize}

\begin{document}
"""
pdoc_preamble=r"""\documentclass[10pt]{article}
\usepackage{fullpage}
\usepackage{multirow}
\usepackage[top=0.1cm, bottom=0.1cm, left=0.1cm, right=0.1cm]{geometry}
\usepackage{floatrow}
\DeclareFloatFont{fontsize}{\small}
\floatsetup[table]{font=fontsize}

\begin{document}
"""

doc_suffix=r"""
\end{document}
"""

### table strings
table_preamble=r"""
\begin{table}
\caption{%s}
\begin{center}
\begin{tabular}{*{%d}{c}}
"""

table_suffix = r"""
\end{tabular}
\end{center}
\end{table}
"""

### column headings
stats_headings = r"""
\cline{1-%d}
\multirow{3}{*}{data set} & $M_{\mathrm{prim}}/M_{\odot}$ & $P_{\mathrm{orb}}$ [%s] & No. modes   & $t/P_{\mathrm{orb}}[0]$         &  $\left<\sum_i A_i^2\right>$ [%s]          & $\left<2\sum_i \gamma_i A_i^2\right>$ [%s/%s]   & \multirow{2}{*}{$\left<\sum_i \mathrm{d}E_i/\mathrm{d}t\right>$ [%s/%s]}   & \multirow{2}{*}{$\left<H_\ast\right>$ [%s]} & \multirow{2}{*}{$\left<H_{int}+H_\ast\right>$ [%s]} & \multirow{2}{*}{$\left<\mathrm{d}(H_{int}+H_\ast)/\mathrm{d}t\right>$ [%s/%s]} \\
\multirow{3}{*}{}         & $R_{\mathrm{prim}}/R_{\odot}$ & $\epsilon$              & No. gens    & $\mathrm{d}t/P_{\mathrm{orb}}$  &  $\sigma_{\sum_i A_i^2}$ [%s]              & $\sigma_{2\sum_i \gamma_i A_i^2}$ [%s/%s]       & \multirow{2}{*}{}                                                          & \multirow{2}{*}{}                           & \multirow{2}{*}{}                                   & \multirow{2}{*}{} \\ 
\multirow{3}{*}{}         & $M_{\mathrm{comp}}/M_{J}$     & $E_{\mathrm{orb}}$ [%s] & No. triples & $t/P_{\mathrm{orb}}[-1]$        &  gini\{ $\left<\sum_i A_i^2\right>$ \}     & gini\{ $\left<2\sum_i \gamma_i A_i^2\right>$ \} & $\sigma_{\sum_i \mathrm{d}E_i/\mathrm{d}t}$ [%s/%s]                        & $\sigma_{H_\ast}$ [%s]                      & $\sigma_{H_{int}+H_\ast}$                            & $\sigma_{\mathrm{d}(H_{int}+H_\ast)/\mathrm{d}t}$ [%s/%s] \\
\cline{1-%d}
"""
#stats_headings = r"""
#\cline{1-%d}
#\multirow{3}{*}{data set} & $M_{\mathrm{prim}}/M_{\odot}$ & $P_{\mathrm{orb}}$ [%s] & No. modes   & $t/P_{\mathrm{orb}}[0]$         &  $\left<\sum_i A_i^2\right>$ [%s]          & $\left<2\sum_i \gamma_i A_i^2\right>$ [%s/%s]   & $\left<\sum_i \mathrm{d}E_i/\mathrm{d}t\right>$ [%s/%s]   & \multirow{2}{*}{$\left<H_\ast\right>$ [%s]} & \multirow{2}{*}{$\left<H_{int}+H_\ast\right>$ [%s]} & \multirow{2}{*}{$\left<\mathrm{d}(H_{int}+H_\ast)/\mathrm{d}t\right>$ [%s/%s]} \\
#\multirow{3}{*}{}         & $R_{\mathrm{prim}}/R_{\odot}$ & $\epsilon$              & No. gens    & $\mathrm{d}t/P_{\mathrm{orb}}$  &  $\sigma_{\sum_i A_i^2}$ [%s]              & $\sigma_{2\sum_i \gamma_i A_i^2}$ [%s/%s]       & $\sigma_{\sum_i \mathrm{d}E_i/\mathrm{d}t}$ [%s/%s]       & \multirow{2}{*}{}                           & \multirow{2}{*}{}                                   & \multirow{2}{*}{} \\ 
#\multirow{3}{*}{}         & $M_{\mathrm{comp}}/M_{J}$     & $E_{\mathrm{orb}}$ [%s] & No. triples & $t/P_{\mathrm{orb}}[-1]$        &  gini\{ $\left<\sum_i A_i^2\right>$ \}     & gini\{ $\left<2\sum_i \gamma_i A_i^2\right>$ \} & gini\{ $\left<\sum_i \mathrm{d}E_i/\mathrm{d}t\right>$ \} & $\sigma_{H_\ast}$ [%s]                      & $\sigma_{H_{int}+H_\ast}$                            & $\sigma_{\mathrm{d}(H_{int}+H_\ast)/\mathrm{d}t}$ [%s/%s] \\
#\cline{1-%d}
#"""
no_stats_cols = 11

lookup_headings = r"""
\cline{1-%d}
\multirow{3}{*}{data set} & logfilename \\
\multirow{3}{*}{}         & outfilename \\
\multirow{3}{*}{}         & stefilename \\
\cline{1-%d}
"""
no_lookup_cols = 2


sum_tablename = opts.output_dir+"/ste-sum%s.tex" % opts.tag
if opts.verbose: print "building summary tables"

### set up strings which will be written to file
if opts.tag != "":
  stats_table = r"""%s
%s""" % (table_preamble % ("Statistics: "+opts.tag[1:].replace("_","\_"), no_stats_cols), stats_headings % (no_stats_cols, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, no_stats_cols) )
else:
  stats_table = r"""%s
%s""" % (table_preamble % ("Statistics", no_stats_cols), stats_headings % (no_stats_cols, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, no_stats_cols) )

stats_row_No = 0

if opts.tag != "":
  lookup_table = r"""%s
%s""" % (table_preamble % ("Filenames: "+opts.tag[1:].replace("_","\_"), no_lookup_cols), lookup_headings % (no_lookup_cols, no_lookup_cols))
else:
  lookup_table = r"""%s
%s""" % (table_preamble % ("Filenames", no_lookup_cols), lookup_headings % (no_lookup_cols, no_lookup_cols))

lookup_row_No = 0

### iterate over data and build strings
for data_set, (stefilename, sdata) in enumerate(dat):
  if opts.verbose: print "working on "+stefilename

  ### extract relevant data
  nmu.set_units(system=sdata["unit_system"])

  logfilename = sdata["system"]["logfilename"]
  Mprim = sdata["system"]["Mprim/Msun"]
  Rprim = sdata["system"]["Rprim/Rsun"]
  Mcomp = sdata["system"]["Mcomp/Mjup"]
  Porb  = nmu.convert_time(sdata["system"]["Porb"], nmu.units["time"], units_time)
  ecc   = sdata["system"]["ecc"]
  Eorb  = nmu.convert_energy(sdata["system"]["Eorb"], nmu.units["energy"], units_energy) ; absEorb = abs(Eorb) ; Eorb = "%.6fe%d" % nmu.float_to_scientific( Eorb )
  N_m   = sdata["system"]["Nmodes"]
  N_g   = sdata["system"]["Ngens"]
  N_t   = sdata["system"]["Ntriples"]

  if sdata.has_key("time_domain"):
    outfilename = sdata["time_domain"]["outfilename"]
    start = sdata["time_domain"]["t_P[0]"]
    stop  = sdata["time_domain"]["t_P[-1]"]
    step  = sdata["time_domain"]["dt_P"]
  else:
    outfilename = start = stop = step = "--" 

  if sdata["stats"].has_key("mean{Hns/|Eorb|}"): # stellar hamiltonian
    mHns = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["mean{Hns/|Eorb|}"]*absEorb )
    sHns = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["stdv{Hns/|Eorb|}"]*absEorb )
  else:
    mHns = sHns = "--"

  if sdata["stats"].has_key("mean{Hint+Hns}/|Eorb|"): # stellar+interaction hamiltonian
    mHintHns = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["mean{Hint+Hns}/|Eorb|"]*absEorb )
    sHintHns = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["stdv{Hint+Hns}/|Eorb|"]*absEorb )
  else:
    mHintHns = sHintHns = "--"

  if sdata["stats"].has_key("mean{d(Hint+Hns)/dt}*Porb/|Eorb|"): # time rate of change of stellar+interaction hamiltonian
    mddtHintHns = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["mean{d(Hint+Hns)/dt}*Porb/|Eorb|"]*absEorb/Porb )
    sddtHintHns = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["stdv{d(Hint+Hns)/dt}*Porb/|Eorb|"]*absEorb/Porb )
  else:
    mddtHintHns = sddtHintHns = "--"

  if sdata["stats"].has_key("mean{|sum{Edot}|*(Porb/|Eorb|)}"): # viscous dissipation
    myE = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["mean{|sum{Edot}|*(Porb/|Eorb|)}"]*(-absEorb/Porb) )
    syE = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["stdv{|sum{Edot}|*(Porb/|Eorb|)}"]*(-absEorb/Porb) )
    gyE = "%.6f" % sdata["stats"]["Gini_index{mean{sum{Edot}}}"]
  else:
    myE = syE = gyE = "--"

  if sdata["stats"].has_key("mean{sum{E}/|Eorb|}"): # mode amplitudes squared
    mE = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["mean{sum{E}/|Eorb|}"]*absEorb )
    sE = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["stdv{sum{E}/|Eorb|}"]*absEorb )
    gE = "%.6f" % sdata["stats"]["Gini_index{mean{sum{E}}}"]
  else:
    mE = sE = gE = "--"

  if sdata["stats"].has_key("mean{sum{dE/dt*Porb}/|Eorb|}"): # time rate of change of mode energies
    mdE = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["mean{sum{dE/dt*Porb}/|Eorb|}"]*absEorb/Porb )
    sdE = "$%.6f\cdot10^{%d}$" % nmu.float_to_scientific( sdata["stats"]["stdv{sum{dE/dt*Porb}/|Eorb|}"]*absEorb/Porb )
    gdE = "%.6f" % sdata["stats"]["Gini_index{mean{sum{dE/dt}}}"]
  else:
    mdE = sdE = gdE = "--"

  ### build lookup table
  lookup_table += r"""
\multirow{3}{*}{%d} & %s \\
\multirow{3}{*}{}   & %s \\
\multirow{3}{*}{}   & %s \\""" % (data_set, logfilename.replace("_","\_"), outfilename.replace("_","\_"), stefilename.replace("_","\_"))
  lookup_row_No += 1

  lookup_table += r"""
\cline{1-%d}""" % no_lookup_cols

  ### build stats table
  stats_table += r"""
\multirow{3}{*}{%d} & %.3f & %.1f & %d & %.1f & %s & %s & %s & \multirow{2}{*}{ %s } & \multirow{2}{*}{ %s } & \multirow{2}{*}{ %s } \\
\multirow{3}{*}{}   & %.3f & %.3f & %d & %.1f & %s & %s & %s & \multirow{2}{*}{}     & \multirow{2}{*}{ }    & \multirow{2}{*}{ }    \\
\multirow{3}{*}{}   & %.3f & %s   & %d & %.1f & %s & %s & %s & %s                    & %s                    & %s                    \\""" % (data_set, Mprim, Porb, N_m, start, mE, myE, mdE, mHns, mHintHns, mddtHintHns, Rprim, ecc, N_g, step, sE, syE, Mcomp, Eorb, N_t, stop, gE, gyE, sdE, sHns, sHintHns, sddtHintHns)
#  stats_table += r"""
#\multirow{3}{*}{%d} & %.3f & %.1f & %d & %.1f & %s & %s & %s & \multirow{2}{*}{ %s } & \multirow{2}{*}{ %s } & \multirow{2}{*}{ %s } \\
#\multirow{3}{*}{}   & %.3f & %.3f & %d & %.1f & %s & %s & %s & \multirow{2}{*}{}     & \multirow{2}{*}{ }    & \multirow{2}{*}{ }    \\
#\multirow{3}{*}{}   & %.3f & %s   & %d & %.1f & %s & %s & %s & %s                    & %s                    & %s                    \\""" % (data_set, Mprim, Porb, N_m, start, mE, myE, mdE, mHns, mHintHns, mddtHintHns, Rprim, ecc, N_g, step, sE, syE, sdE, Mcomp, Eorb, N_t, stop, gE, gyE, gdE, sHns, sHintHns, sddtHintHns)
  stats_row_No += 1

  stats_table +=  r"""
\cline{1-%d}""" % no_stats_cols

  if (stats_row_No+1)%14 == 0:
    stats_table += table_suffix
    if opts.tag != "":
      stats_table += r"""%s
%s""" % (table_preamble % ("Statistics: "+opts.tag[1:].replace("_","\_"), no_stats_cols), stats_headings % (no_stats_cols, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, no_stats_cols) )
    else:
      stats_table += r"""%s
%s""" % (table_preamble % ("Statistics", no_stats_cols), stats_headings % (no_stats_cols, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, units_energy, units_energy, units_time, units_energy, units_time, units_energy, units_energy, units_energy, units_time, no_stats_cols) )

  if (lookup_row_No+1)%14 == 0:
    lookup_table += table_suffix
    if opts.tag != "":
      lookup_table += r"""%s
%s""" % (table_preamble % ("Filenames: "+opts.tag[1:].replace("_","\_"), no_lookup_cols), lookup_headings % (no_lookup_cols, no_lookup_cols))
    else:
      lookup_table += r"""%s
%s""" % (table_preamble % ("Filenames", no_lookup_cols), lookup_headings % (no_lookup_cols, no_lookup_cols))


### finish tables
lookup_table += table_suffix
stats_table += table_suffix

### write document
if opts.verbose: print "writing summary tables to:\n\t%s" % (sum_tablename)

sum_tablefile = open(sum_tablename, "w")
print >> sum_tablefile, ldoc_preamble
print >> sum_tablefile, stats_table
print >> sum_tablefile, lookup_table
print >> sum_tablefile, doc_suffix
sum_tablefile.close()


