#!python_alias
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
parser.add_option("", "--D2-bokeh", default=False, action="store_true")

parser.add_option("", "--special-plots", default=False, action="store_true", help="generate plots that do not neatly fall in the D1, D2 catagories")

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

      if NgensD.keys() != []: ### only produce the plot if it isn't stupid to do so

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

##################################################
if opts.D2_bokeh:
  colors = ['blue','red','green','cyan','magenta','yellow','black']

  import bokeh
  from bokeh import plotting as bplt
  from bokeh.objects import HoverTool
  from collections import OrderedDict

  TOOLS="pan,box_zoom,reset,hover"
  plot_width=600
  plot_height=300 

  if opts.verbose: print "2D_bokeh"

  ### set up html page
  filename = opts.cachefilename + "%s.html"%opts.tag
  bplt.output_file(filename)

  ### loop over all pairs
  for ystr, ylabel, ykey, yvkey in known_y:
    if opts.verbose: print "\t", ystr
    for xkey, xlabel in known_x:
      if opts.verbose: print "\t\t", xkey

      ### loop and pull out data
      all_x = []
      all_y = []
      all_color = []

      all_yv = []

      all_stefilename = []
      all_Nmodes = []
      all_Ntriples = []
      all_Ngens = []

      for stefilename, sdata in dat:
        try:
          y = sdata["stats"][ykey]
        except KeyError:
          print "could not find data for %s in %s... skipping" % (ykey, stefilename)
          continue

        all_stefilename.append( stefilename )

        if "Gini_index" in ykey:
          y = 1-np.array(y)

        Ngens = sdata["system"]["Ngens"]
        Nmodes = sdata["system"]["Nmodes"]
        Ntriples = sdata["system"]["Ntriples"]

        x = sdata["system"][xkey]

        all_x.append( x )
        all_y.append( y )
        all_color.append( colors[Ngens] )
        all_Ngens.append( Ngens )
        all_Nmodes.append( Nmodes )
        all_Ntriples.append( Ntriples )

        if yvkey:
          all_yv.append( sdata["stats"][yvkey] )

      ### check for data
      if not len(all_y):
        continue

      fig = bplt.figure()
      bplt.hold()

      ### build data source for hover tool
      source = bplt.ColumnDataSource(data=dict(x=all_x, y=all_y, stefilename=all_stefilename, Nmodes=all_Nmodes, Ngens=all_Ngens, Ntriples=all_Ntriples))

      ### plot circle glyphs
      bplt.circle(all_x, all_y, source=source, tools=TOOLS, fill_color=None, fill_alpha=0.6, line_color=all_color, Title="%s vs %s"%(xlabel,  ylabel), plot_width=plot_width, plot_height=plot_height)
#      bplt.circle(all_x, all_y, radius=radii, source=source, tools=TOOLS, fill_color=None, fill_alpha=0.6, line_color=all_color, Title="%s vs %s"%(xlabel,  ylabel))

      ### annotate circle glyphs
#      text(x, y, text=inds, alpha=0.5, text_font_size="5pt", text_baseline="middle", text_align="center", angle=0)

      ### find hover tool, and tell it what to look for
      hover = [t for t in bplt.curplot().tools if isinstance(t, HoverTool)][0]
      hover.tooltips = OrderedDict([ 
                                    #("index", "$index"), 
                                    ("(x,y)", "($x, $y)"),
                                    ("stefilename", "@stefilename"),
                                    ("Nmodes","@Nmodes"),
                                    ("Ntriples","@Ntriples"),
                                    ("Ngens","@Ngens"),
                                   ])

      ### label axes
      bplt.xaxis().axis_label=xlabel
      bplt.yaxis().axis_label=ylabel
      

#  bplt.show() 

#=================================================
if opts.special_plots:

  axes_loc = [0.15,0.15,0.7,0.8]
  cb_loc = [0.9, 0.15, 0.05, 0.8]

  ### Ng1 vs Ng2
  if opts.verbose: print "Ng1 vs Ng2"

  Ng1 = []
  Ng2 = []
  Ntriples = []
  mEdot = []
  vEdot = []

  ### pull out data
  for stefilename, sdata in dat:
    ### pull out damping estimates
    try:
      mEdot.append( sdata["stats"]["mean{|sum{Edot}|*Porb/|Eorb|}"] )
      vEdot.append( sdata["stats"]["stdv{|sum{Edot}|*Porb/|Eorb|}"] )
    except KeyError:
      print "could not find data for mean{|sum{Edot}|*Porb/|Eorb|} in %s... skipping" % (stefilename)
      continue

    ### pull out number of modes in each generation
    Ngens = sdata["system"]["Ngens"]
    Ngi = sdata["system"]["Ngi"]
    if Ngens < 2:
      Ng1.append( 0 )
      Ng2.append( 0 )
    elif Ngens < 3:
      Ng1.append( Ngi[1] )
      Ng2.append( 0 )
    else:
      Ng1.append( Ngi[1] )
      Ng2.append( Ngi[2] )

    ### pull out Ntriples
    Ntriples.append( sdata["system"]["Ntriples"] )

  ### convert mEdot into sizes
  cmap_norm = matplotlib.colors.Normalize(vmin=min(mEdot), vmax=max(mEdot))
  min_size = 1
  max_size = 10
  sizes = [ min_size + (max_size-min_size)*cmap_norm( Edot ) for Edot in mEdot ]

  ### convert Ntriples into colors
  Ntmax = max(Ntriples)
  Ntmin = min(Ntriples)
  cmap_norm = matplotlib.colors.Normalize(vmin=min(Ntriples), vmax=max(Ntriples))
  cmap = plt.get_cmap("jet")
  colors = [ cmap( cmap_norm(Nt) ) for Nt in Ntriples]

  ### generate plot
  fig = plt.figure()
  ax = fig.add_axes(axes_loc)
  ax_cb = fig.add_axes(cb_loc)

  ax.plot(Ng1, Ng2, marker="o", markeredgecolor=colors, markerfacecolor=colors, markersize=sizes)
  cb = matplotlib.colorbar.ColorbarBase(ax_cb, cmap=cmap, norm=cmap_norm, orientation='vertical')

  ax.set_xlabel("$N_{G1}$")
  ax.set_ylabel("$N_{G2}$")
  cb.set_label(r"$\frac{\sum 2\gamma_i A_i^2}{|E_\mathrm{orb}|/P_\mathrm{orb}}$")
 
  ax.grid(True)
 
  figname = "%s/Ng1-Ng2%s.png" % (opts.output_dir, opts.tag)
  if opts.verbose: print "\t\t\tsaving ", figname
  fig.savefig(figname)
  plt.close(fig)

  
