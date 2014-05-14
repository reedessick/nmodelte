usage="""a module for network diagnostics and visualizations"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from collections import defaultdict

import nmode_utils as nmu
import nmode_state as nms
import mode_selection as ms

####################################################################################################
#
#
#           binning/sorting functions
#
#
####################################################################################################
def bin_by_x(bins, x):
  """ 
  helper function for other binning functions. requires x to have the form: [(modeNo_i, x_i), ...]
  also requires x to be sorted
  """
  Nm = len(x)

  ind = 0
  while (x[ind][1] < bins[0]) and (ind < Nm):
    ind += 1

  binned = []
  for end in bins[1:]:
    this_bin = []
    while (ind < Nm) and (x[ind][1] < end):
      this_bin.append(x[ind][0])
      ind += 1
    binned.append(this_bin)

  return binned

##################################################
#def bin_by_w(bins, network):
#  """
#  places modes into bins using their natural frequencies.
#  bins is a list of bin edges
#  returns a list of mode numbers corresponding to each bin
#  """
#  ws = [(modeNo, mode.w) for modeNo, mode in enumerate(network.modes)]
#  ws.sort(key=lambda l: l[1])
# 
#  return bin_by_x(bins, ws)
#
##################################################
#def bin_by_l(bins, network):
#  """
#  places modes into bins using their angular quantum number
#  bins is a list of bin edges
#  returns a list of mode numbers corresponding to each bin
#  """
#  ls = [(modeNo, mode.l) for modeNo, mode in enumerate(network.modes)]
#  ls.sort(key=lambda l: l[1])
#
#  return bin_by_x(bins, ls)
#
####################################################################################################
#
#
#  stacked histograms
#
#
####################################################################################################
def stacked_histogram(bins, xdata, data, log=False):
  """
  generates a stacked histogram of data
  data must have the following form: [(datum, datum, datum)), ...]
  where 
    len(data) = N_m
    xdata is the data for the x axis
    len(data[i]) = number of levels in stacked histogram + 1
  """

  depth = len(data[0]) # number of stacked plots

  fig = plt.figure()
  axs = []

  ax_width = 0.8
  ax_height = 0.8/depth
  buffer = 0.01

  # iterate and build plots
  for d in range(depth):
    ax = fig.add_axes([0.15, 0.95-(d+1)*ax_height, ax_width, ax_height-buffer])
    n, _, _ = ax.hist(xdata, bins, weights=[l[d] for l in data], histtype="step", log=log)

    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_ylim(ymax=1.1*max(n)) ### auto-scaling appears to mess up occassionally
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs

##################################################
def generational_stacked_histogram(network, bins, xdata, data, log=False):
  """
  generates as stacked histogram of data, separating data by generation
  data must have the following form: [(datum, datum, datum)), ...]
  where 
    len(data) = N_m
    xdata is the data for the x axis
    len(data[i]) = number of levels in stacked histogram + 1

  IMPORTANTLY: assumes the order of modeNos in network is the same as that in xdata, data

  """

  ### divide network by generation.
  N_m = len(network)

  gens, coups = network.gens()

  ### build figure
  depth = len(data[0]) # number of stacked plots

  fig = plt.figure()
  axs = []

  ax_width = 0.8
  ax_height = 0.8/depth
  buffer = 0.01

  # iterate and build plots
  xdata_by_gen = []
  for gen in gens:
    this_gen = []
    for modeNo in gen:
      x = xdata[modeNo]
      if x < bins[0]: # in case binning is poor...
        x = bins[0]
      elif x > bins[-1]: # in case binning is poor...
        x = bins[-1]
      this_gen.append( x )
    xdata_by_gen.append( this_gen )
#  xdata_by_gen = [ [xdata[modeNo] for modeNo in gen] for gen in gens ]

  for d in range(depth):
    ax = fig.add_axes([0.15, 0.95-(d+1)*ax_height, ax_width, ax_height-buffer])
    ax.hist(xdata_by_gen, bins, weights=[[data[modeNo][d] for modeNo in gen] for gen in gens], histtype="step", log=log)

    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs


####################################################################################################
#
#
#  coupling diagram
#
#
####################################################################################################
def coupling_hist(network, gens=False, num_bins=50, log=False, genNos=False):
  """
  histogram of the number of couplings belonging to any given mode, broken down by generation
  """
  if not gens:
    gens, _ = network.gens()
  if not genNos:
    genNos = range(len(gens))
  n_k = []
  labels = []
  for genNo in genNos:
    n_k.append( [len(network.K[modeNo]) for modeNo in gens[genNo]] )
    labels.append( "gen %d" % genNo )

  fig = plt.figure()
  ax = plt.subplot(1,1,1)
  ax.hist(n_k, num_bins, histtype="step", stacked=True, log=log, label=labels)

  if log:
    ax.set_ylim(ymin=0.90*ax.get_ylim()[0])

  ax.set_xlabel("No. couplings")
  ax.set_ylabel("count")

  return fig, ax

##################################################
def coupling_tree_nl(system, tree_type="placement", genNos=[], verbose=False, mode_colors="b", mode_order=None, conn_colors="r", fig_ax=None, mode_nums=[]):
  """
  coupling diagram in the n-l plane via delegation through __coupling_tree_nl.
    tree_type determines the types of connections drawn between modes
      currently supports : placement, siblings, shared_parent, shared_child, triples
    genNos determines which generations are included (empyt list means we plot everything)
  """

  if isinstance(mode_colors, str):
    mode_colors = [mode_colors for _ in system.network.modes]

  ## break down into generations
  if verbose: print "\t\tbreaking network into generations"
  network = system.network
  gens, coups = network.gens()
  N_g = len(gens)
  if genNos == []: # empty list means we plot all generations
    genNos = range(N_g)

  modeNos = []
  conns = []
  ### determine the type of plot
  if verbose: print "\t\tdetermining which modes/conns to plot based on tree_type=%s" % tree_type

  if tree_type == "placement": ## only the location of selected modes
    for genNo in genNos:
      modeNos += gens[genNo]

  elif tree_type == "triples":
    old_modeNos = set()
    old_conns = set()
    for genNo in genNos:
      if genNo < N_g-1:
        for o,i,j,_ in coups[genNo]: # genNo modes are parents --> add everything they touch
          old_modeNos.add(o)
          old_modeNos.add(i)
          old_modeNos.add(j)
          conns += [(o,i), (o,j), (i,j)]
      if (genNo>0) and (genNo-1 not in genNos):
        for o,i,j,_ in coups[genNo-1]: # genNo modes are children --> add everything they touch
          old_modeNos.add(o)
          old_modeNos.add(i)
          old_modeNos.add(j)
          old_conns.add((o,i))
          old_conns.add((o,j))
          old_conns.add((i,j))
    modeNo_map = {}
    for new_modeNo, old_modeNo in enumerate(old_modeNos): # build modeNo_map and fill out modes
      modeNo_map[old_modeNo] = new_modeNo
      modeNos.append( old_modeNo )
    for old_modeNo1, old_modeNo2 in old_conns: # convert old_conns with new_modeNo
      conns.append( (modeNo_map[old_modeNo1], modeNo_map[old_modeNo2]) )

  elif tree_type == "siblings": # modes in the same generation directly connected to eachother.
    modeNo_map = {}
    new_modeNo = 0
    for genNo in genNos:
      for modeNo in gens[genNo]:
        modeNos.append( modeNo )
        modeNo_map[modeNo] = new_modeNo
        new_modeNo += 1
      if genNo > 0: # gen0 has no siblings by this definition
        for _,i,j,_ in coups[genNo-1]:
          conns.append( (modeNo_map[i], modeNo_map[j]) )

  elif tree_type == "shared_parent": # modes in the same generation that share a parent
    modeNo_map = {}
    new_modeNo = 0
    for genNo in genNos:
      for modeNo in gens[genNo]:
        modeNos.append( modeNo )
        modeNo_map[modeNo] = new_modeNo
        new_modeNo += 1
      if genNo > 0: # gen0 has no parent
        parent_map = defaultdict( set )
        for o,i,j,_ in coups[genNo-1]: # the couplings that have this genNo as children
          parent_map[o].add(i) # add children to a set defined by the parent
          parent_map[o].add(j)
        for children in parent_map.values(): # map the lists of children which share parents into connections
          children = list(children)
          N_children = len(children)
          for ind1 in range(N_children):
            new_modeNo1 = modeNo_map[children[ind1]]
            for ind2 in range(ind1+1, N_children): # don't connect modes to themselves because that's stupid in this situation
              conns.append( (new_modeNo1, modeNo_map[children[ind2]]) )

  elif tree_type == "shared_child": # modes in the same generation that share a child
    modeNo_map = {}
    new_modeNo = 0
    for genNo in genNos:
      for modeNo in gens[genNo]:
        modeNos.append( modeNo )
        modeNo_map[modeNo] = new_modeNo
        new_modeNo += 1
      if genNo < N_g-1: # gen(N_g-1) has no children
        child_map = defaultdict( set )
        for o,i,j,_ in coups[genNo]: # these couplings are when genNo is the parent
          child_map[i].add(o) # add the parents to sets labeled by the daughters
          child_map[j].add(o)
        for parents in child_map.values(): # map the list of parents which share children into connections
          parents = list(parents)
          N_parents = len(parents)
          for ind1 in range(N_parents):
            new_modeNo1 = modeNo_map[parents[ind1]]
            for ind2 in range(ind1+1, N_parents): # don't connect modes to themselves because that's stupid in this situation
              conns.append( (new_modeNo1, modeNo_map[parents[ind2]]) )

  else:
    raise ValueError("unkown tree_type=%s in nmode_diagnostic.coupling_tree_nl" % tree_type)

  if mode_order != None: # plot modes in a certain order
    modeNos = [modeNo for modeNo in mode_order if (modeNo in modeNos)] # rearange modeNos by mode_order
    modeNo_map = dict( (modeNo, ind) for ind, modeNo in enumerate(modeNos) )
    conns = [(modeNo_map[modeNo1], modeNo_map[modeNo2]) for modeNo1, modeNo2 in conns]

  if mode_nums: ### down-select which modes are plotted
    modeNos = [modeNo for modeNo in mode_nums if (modeNo in modeNos)]
    modeNo_map = dict( (modeNo, ind) for ind, modeNo in enumerate(modeNos) )
    conns = [(modeNo_map[modeNo1], modeNo_map[modeNo2]) for modeNo1, modeNo2 in conns if (modeNo_map.has_key(modeNo1) and modeNo_map.has_key(modeNo2))]

  return __coupling_tree_nl([system.network.modes[modeNo] for modeNo in modeNos], mode_colors=[mode_colors[modeNo] for modeNo in modeNos], conns=conns, conn_colors=conn_colors, verbose=verbose, fig_ax=fig_ax)
#  return __coupling_tree_nl(modes, mode_colors=__mode_colors, conns=conns, conn_colors=conn_colors, verbose=verbose, fig_ax=fig_ax)

##################################################
def __coupling_tree_nl(modes, mode_colors="b", conns=[], conn_colors="r", verbose=False, fig_ax=None, mode_alpha=0.75, edge_alpha=0.50):
  """ 
  a plotting function for coupling_tree_nl. It plots the modes on the n-l plane (with offsets from l described by m) and draws connections between them described in conns)
  Because we can come up with many different lines to connect the modes, this method should be useful through delegation 

    modes = [mode1, mode2, mode3,...] where each modeN is an instance of networks.mode
    conns = [(modeNo1, modeNo2), ...] where modeNoN refers to the order in which modes live in "modes"

    fig_ax = (fig, ax) will allow you to plot on an existing axis (belonging to fig)

  warning! this plot does not distinguish between frequency signs (+/-), so there may be some confusion. This should be minimized by choosing only one frequency sign for the network's matriarch
  """
  if isinstance(mode_colors, str):
    mode_colors = [mode_colors for _ in modes]
  if isinstance(conn_colors, str):
    conn_colors = [conn_colors for _ in conns]

  if not fig_ax:
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
  else:
    fig, ax = fig_ax

  ### compute positions of each mode in modes
  if verbose: print "\t\tcomputing mode locations"
  x = [] # l+(m/2*l)
  y = [] # n
  for mode in modes:
    x.append( mode.l + 0.25*mode.m/mode.l )
    y.append( mode.n )

  ### draw lines between points
  if verbose: print "\t\tdrawing edges"
  for (modeNo1, modeNo2), color in zip(conns, conn_colors):
    _x, _y = nl_edge(x[modeNo1], y[modeNo1], x[modeNo2], y[modeNo2])
#    print _x, _y
    ax.plot(_x, _y, color=color, alpha=edge_alpha)

  ### draw modes
  if verbose: print "\t\tdrawing modes"
  for (_x, _y, color) in zip(x,y,mode_colors):
    ax.plot(_x, _y, marker="o", markersize=2, markeredgecolor=color, markerfacecolor=color, linestyle="none", alpha=mode_alpha)

  ax.set_xlabel(r"$l+m/4l$")
  ax.set_ylabel(r"$n$")

  if x: ### we're actually plotting something
    ax.set_xlim(xmin=min(x)-0.25, xmax=max(x)+0.25)
    ax.set_ylim(ymin=min(y)*0.90, ymax=max(y)*1.05)
  else:
    ax.set_xlim(xmin=-0.25, xmax=+0.25)
    ax.set_ylim(ymin=0, ymax=1)

  return fig, ax

#########################
def nl_edge(x1,y1, x2,y2, num_samples=101):
  """ returns a set of points for the edge between x1,y1 and x2,y2 for coupling_tree_nl """
  if (x1 == x2) and (y1 == y2): # trivial...
    return [x1], [y1] # only return something for syntactic reasons

  elif y1 == y2: # we use an arc segment when y1==y2
    return __nl_arc(x1,x2, y1, num_samples, 0.1)

  elif x1 == x2: # we use an arc segment when x1==x2
    y, x = __nl_arc(y1,y2, x1, num_samples, 0.1)
    return x, y

  else: ## simple 3rd order spline with zero derivatives in the x-direction 
    return __nl_spline(x1,x2, y1,y2, num_samples)

### 
def __nl_spline(x1,x2,y1,y2,num_samples):
  dx = np.linspace(0,1,num_samples)
  x = x1 + (x2-x1)*dx
  return x, y1 + 3*(y2-y1)*dx**2 - 2*(y2-y1)*dx**3

###
#def __nl_arc(x1,x2,y1, num_samples, f):
#    Dx = x2-x1
#    R = Dx*(f**2+0.25)/(2*f)
#    x = np.linspace(x1,x2,num_samples)
#    return x, y1+(R**2 - (x-(x1+x2)/2)**2)**0.5 - (R**2 - Dx**2/4.0)**0.5

def __nl_arc(x1,x2, y1, num_samples, f):
    Sx = (x1+x2)/2.0
    Dx = (x2-x1)/2.0
    r = Dx/Sx
    R = r*(f**2+0.25)/f
    x = np.linspace(x1,x2,num_samples)
    dy = (R**2 - (x/Sx-1)**2)**0.5 - (R**2 - r**2)**0.5
    return x, y1*(1+dy)


##################################################
def coupling_tree_w(system, verbose=False):
  """
  generate a diagram representing the network's structure

  places modes on a line by generation, and then connects modes with lines (representing couplings)
  """
  network = system.network
  Oorb = system.Oorb

  if verbose: print "\tsplitting network into generations"
  # divide network by generation.
  N_m = len(network)

  gens, coups = network.gens()
  maxNm = max([len(gen) for gen in gens]) # maximum number of modes in any one generation
  
  N_g = len(gens) # number of generations

  if verbose: print "\tdefining edges"
  ### define sets of connections with line weights (number of connections they represent)
  edges = []
  for coup in coups:
    pc = {}
    cc = {}
    for p,c1,c2,k in coup:
      # parent child1 relationship
      if pc.has_key((p,c1)):
        pc[(p,c1)] += 1
      else:
        pc[(p,c1)] = 1
      # parent child2 relationship
      if pc.has_key((p,c2)):
        pc[(p,c2)] += 1
      else:
        pc[(p,c2)] = 1
      # child1 child2 relationship
      if cc.has_key((c1,c2)):
        cc[(c1,c2)] += 1
      else:
        cc[(c1,c2)] = 1
    edges.append( [pc, cc] )       

  if verbose: print "\tdefining bounds for each generation"
  ### define boundaries for each generation: w -> position
  bounds = []
  flip = False # controls the orientation of the axes
  for gen in gens:
    mi = np.infty
    ma = -np.infty
    for modeNo in gen:
      w = network.modes[modeNo].w
      if mi > w:
        mi = w
      if ma < w:
        ma = w

    s = (ma+mi)/2.
    d = (ma-mi)/2.
    if d > 0:
      f = 1.10 # 5% buffer on ranges
      if flip:
        bounds.append( (s+f*d, s-f*d) )
      else:
        bounds.append( (s-f*d, s+f*d) )
    else:
      if flip:
        bounds.append( (1.5*s, 0.5*s) )
      else:
        bounds.append( (0.5*s, 1.5*s) )
    flip = (not flip) 

  if verbose: print "\tbuilding figure"
  ### build figure
  fig = plt.figure()
  ax_width = 0.95
  ax_height = 0.95
  xmin, xmax = 0, 1
  ymin, ymax = -0.2, N_g-0.9

  figheight=min(max([N_g*3, maxNm/3., 20]), 100)
  figwidth =min(max(maxNm*0.1, 40), 300)

  ax = fig.add_axes([(1-ax_width)/2., (1-ax_height)*3/5., ax_width, ax_height])
  ax.axison = False

  if verbose: print "\t\taxes"
  ### draw and annotate axes
  n_ticks = max(5, maxNm/80)
  axis_color='k'
  axis_alpha=0.75
  major_tick_width=0.04
  label_fontsize = min(4*N_g, 15)
  label_ha = 'center'
  label_va = 'top'
  alabel_va = 'bottom'

  for genNo in range(N_g):
    ax_x = N_g-genNo-1
#    ax.plot([ymin, ymax], ax_x*np.ones((2,)), color=axis_color, alpha=axis_alpha) # draw axis
    ### annotate !
#    tick_xmin = ax_x - major_tick_width/2 # xlimits for ticks
#    tick_xmax = ax_x + major_tick_width/2
    label_x = ax_x - major_tick_width
    alabel_x = ax_x + major_tick_width
    alabel_y = 0.01

    amin, amax = bounds[genNo] # compute placement of ticks
    s = (amin+amax)/2. # center
    d = abs(amin-amax)/2. # range/2
    dd = 2*d/n_ticks # tick spacing
    pow = int(np.floor(np.log10(dd/Oorb)))
    dd = 10**( pow )*Oorb
    if d/dd > n_ticks:
      dd *= 2
   
    nice_s = round(s/Oorb, 1)*Oorb
    t = nice_s
    while t < s+d:
      y = interp(t, amin, amax)
#      ax.plot([y,y], [tick_xmin, tick_xmax], color=axis_color, alpha=axis_alpha) # plot tick
      ax.text(y, label_x, repr(round(t/Oorb,-pow)), ha=label_ha, va=label_va, fontsize=label_fontsize) # tick value
      t += dd
    t = nice_s - dd
    while t > s-d:
      y = interp(t, amin, amax)
#      ax.plot([y,y], [tick_xmin, tick_xmax], color=axis_color, alpha=axis_alpha)
      ax.text(y, label_x, repr(round(t/Oorb,-pow)), ha=label_ha, va=label_va, fontsize=label_fontsize)
      t -= dd

    ax.text(alabel_y, alabel_x, r"$\omega/\Omega_{\mathrm{orb}}$", ha=label_ha, va=alabel_va, fontsize=label_fontsize) # total axis label

  if verbose: print "\t\tedges"
  ### draw edges
  linewidth_prefactor=0.2
  edge_color='b'
  edge_alpha=0.2
  for genNo, [pc, cc] in enumerate(edges): # genNo refers to parent generation
    p_x = N_g-genNo-1
    c_x = N_g-genNo-2
    pmin, pmax = bounds[genNo] # boundaries for axes
    cmin, cmax = bounds[genNo+1]

    for (p,c), N_k in pc.items(): # parent-child edges
      yp = interp(network.modes[p].w, pmin, pmax) # look up positions on axes
      yc = interp(network.modes[c].w, cmin, cmax)
      x, y = pc_edge(p_x, yp, c_x, yc) # generate edge (spline)
      pc[(p,c)] = [N_k, ax.plot(y,x,linewidth=linewidth_prefactor*N_k, color=edge_color, alpha=edge_alpha)] # plot and record handle for line

    for (c1,c2), N_k in cc.items(): # child-child edges
      yc1 = interp(network.modes[c1].w, cmin, cmax) # look up positions on axis
      yc2 = interp(network.modes[c2].w, cmin, cmax)
      x, y = cc_edge(c_x, yc1, yc2) # generate edge
      cc[(c1,c2)] = [N_k, ax.plot(y,x,linewidth=linewidth_prefactor*N_k, color=edge_color, alpha=edge_alpha)]
       
  ax.set_xlim(xmin=xmin, xmax=xmax)
  ax.set_ylim(ymin=ymin, ymax=ymax)

  plt.setp(fig, figwidth=figwidth, figheight=figheight)

  return fig, ax, gens, coups, edges

####################################################################################################
def coupling_tree_nlm(system, verbose=False):
  """
  generate a diagram representing the network's structure

  places modes on a line by generation, and then connects modes with lines (representing couplings)
  """
  network = system.network

  if verbose: print "\tsplitting network into generations"
  # divide network by generation.
  N_m = len(network)

  gens, coups = network.gens()
  maxNm = max([len(gen) for gen in gens]) # maximum number of modes in any one generation

  N_g = len(gens) # number of generations

  if verbose: print "\tdefining edges"
  ### define sets of connections with line weights (number of connections they represent)
  edges = []
  for coup in coups:
    pc = {}
    cc = {}
    for p,c1,c2,k in coup:
      # parent child1 relationship
      if pc.has_key((p,c1)):
        pc[(p,c1)] += 1
      else:
        pc[(p,c1)] = 1
      # parent child2 relationship
      if pc.has_key((p,c2)):
        pc[(p,c2)] += 1
      else:
        pc[(p,c2)] = 1
      # child1 child2 relationship
      if cc.has_key((c1,c2)):
        cc[(c1,c2)] += 1
      else:
        cc[(c1,c2)] = 1
    edges.append( [pc, cc] )

  if verbose: print "\tbuilding figure"
  ### build figure
  fig = plt.figure()
  ax_width = 0.95
  ax_height = 0.95
  xmin, xmax = 0, 1
  ymin, ymax = -0.2, N_g-0.9

  figheight=min(max([N_g*3, maxNm/3., 20]), 100)
  figwidth =min(max(maxNm*0.1, 40), 300)

  ax = fig.add_axes([(1-ax_width)/2., (1-ax_height)*3/5., ax_width, ax_height])
  ax.axison = False

  if verbose: print "\t\taxes"
  ### draw and annotate axes
  n_ticks = max(5, maxNm/80)
  axis_color='k'
  axis_alpha=0.75
  major_tick_width=0.04
  label_fontsize = min(4*N_g, 15)
  label_ha = 'center'
  label_va = 'top'
  alabel_va = 'bottom'

  ### draw axes, annotate points and build loc maps
  loc_maps = []
  for genNo in range(N_g):
    gen = gens[genNo]
    no_modes = len(gen)
    ax_x = N_g-genNo-1 # central location for annotations of this axis
    _data = []
    for modeNo in gen:
      n,l,m = network.nlm[modeNo] # pull data
      absW = abs(network.wyU[modeNo][0])
      _data.append( (modeNo, absW, l, m) ) # add to _data in order to build loc_map

    # build loc_map
    _data.sort(key=lambda l: l[0]) # sort by abs(w)
    _data.sort(key=lambda l: l[3]) # sort by m
    _data.sort(key=lambda l: l[2]) # sort by l
    loc_map = dict( [ (modeNo, 1.0*(order+1)/(no_modes+1)) for order, (modeNo,_,_,_) in enumerate(_data) ] )

    for modeNo in gen: # annotate using loc_map
      y = loc_map[modeNo]
      ax.text(y, ax_x, "%d\n%d\n%d" % network.nlm[modeNo], ha='center', va='center', fontsize=label_fontsize, alpha=axis_alpha, color=axis_color)

    loc_maps.append( dict( [ (modeNo, 1.0*(order+1)/(no_modes+1)) for order, (modeNo,_,_,_) in enumerate(_data) ] ) ) # add to list

  if verbose: print "\t\tedges"
  ### draw edges
  linewidth_prefactor=0.2
  edge_color='b'
  edge_alpha=0.5
  for genNo, [pc, cc] in enumerate(edges): # genNo refers to parent generation
    p_x = N_g-genNo-1
    c_x = N_g-genNo-2

    p_locmap = loc_maps[genNo] # maps for location of modes
    c_locmap = loc_maps[genNo+1]

    for (p,c), N_k in pc.items(): # parent-child edges
      yp = p_locmap[p]
      yc = c_locmap[c]

      x, y = pc_edge(p_x, yp, c_x, yc) # generate edge (spline)
      pc[(p,c)] = [N_k, ax.plot(y,x,linewidth=linewidth_prefactor*N_k, color=edge_color, alpha=edge_alpha)] # plot and record handle for line

    for (c1,c2), N_k in cc.items(): # child-child edges
      yc1 = c_locmap[c1]
      yc2 = c_locmap[c2]

      x, y = cc_edge(c_x, yc1, yc2) # generate edge
      cc[(c1,c2)] = [N_k, ax.plot(y,x, linewidth=linewidth_prefactor*N_k, color=edge_color, alpha=edge_alpha)]

  ax.set_xlim(xmin=xmin, xmax=xmax)
  ax.set_ylim(ymin=ymin, ymax=ymax)

  plt.setp(fig, figwidth=figwidth, figheight=figheight)

  return fig, ax, gens, coups, edges

##################################################
def interp(x, xmin, xmax):
  """ simple linear interpolation. gives value of y(x) assuming y(xmin)=0, y(xmax)=1 """
  if xmin == xmax:
    raise ValueError("xmin and xmax cannot be equal")
  return (x-xmin)/(xmax-xmin)

##################################################
def pc_edge(x1,y1,x2,y2, num_samples=101):
  """ returns a series of points (x,y) that define the spline interpolation between (x1,y1) <--> (x2,y2) such that dy/dx=0 at x1,x2 """
  dx = np.linspace(0,1,num_samples)
  x = x1 + (x2-x1)*dx
  return x, y1 + 3*(y2-y1)*dx**2 - 2*(y2-y1)*dx**3

##################################################
def cc_edge(x1,y1,y2, num_samples=101):
  """ returns a series of points (x,y) that define the edge connecting two children (a circle) """
  if y1 == y2:
    return cc_self_edge(x1,y1)

  if y1 > y2:
    return cc_edge(x1,y2,y1)

  s = (y1+y2)/2.
  d = abs(y1-y2)/2.

  dx = np.linspace(0,1,num_samples)

  top = (1-dx**2)**0.5
  bot = -(1-(dx-1)**2)**0.5

  return x1-d*np.concatenate([-dx, -1+dx]), s+d*np.concatenate([top,bot])

##################################################
def cc_self_edge(x1,y1, num_samples=101):
  """ allows for a special "self edge" to be defined for self-coupled daughters """
  return cc_edge(x1, y1, y1*1.0001, num_samples=num_samples)

##################################################
def couling_tree_highlight(modeNo, gens, coups, edges):
  """
  alters the figure in place (via line objects in edges) and highlights a specified mode
  """
  highlight_color = 'r'
  highlight_alpha = 1.0
  highlight_zorder = 10

  background_color = 'b'
  background_alpha = 0.1
  background_zorder = 1

  for genNo, gen in enumerate(gens): # select the correct generation
    if modeNo in gen:
      break

  ptriples = [] # select relevant triples
  for p,i,j,k in coups[genNo]: # the mode will ONLY be a parent in these couplings
    if p == modeNo:
      ptriples.append(sorted((i,j)))
  if genNo > 0: # if this mode has a parent, we find those couplings too
    ctriples = []
    for p,i,j,k in coups[genNo-1]:
      if (i == modeNo):
        ctriples.append((p,j))
      elif (j == modeNo):
        ctriples.append((p,i))

  # iterate through and change edge properties
  for edgeNo, (pc, cc) in enumerate(edges):

    if edgeNo == genNo: # mode will only be a parent here
      for (p,c), [N_k, [edge]] in pc.items(): # parent-child connection
        if p == modeNo: # takes care of both p-c branches
          edge.set_color(highlight_color)
          edge.set_alpha(highlight_alpha)
          edge.set_zorder(highlight_zorder)
        else:
          edge.set_color(background_color)
          edge.set_alpha(background_alpha)
          edge.set_zorder(background_zorder)

      for (c1,c2), [N_k, [edge]] in cc.items(): # corresponding child-child connections
        if (sorted((c1,c2)) in ptriples): # takes care of c-c braches corresponding to this parent mode
          edge.set_color(highlight_color)
          edge.set_alpha(highlight_alpha)
          edge.set_zorder(highlight_zorder)
        else:
          edge.set_color(background_color)
          edge.set_alpha(background_alpha)
          edge.set_zorder(background_zorder)

    elif edgeNo == genNo-1: # mode will be a child here
      for (p,c), [N_k, [edge]] in pc.items(): # parent-child connection
        if c == modeNo: # takes care of p-c when this mode is the child
          edge.set_color(highlight_color)
          edge.set_alpha(highlight_alpha)
          edge.set_zorder(highlight_zorder)
        elif (p,c) in ctriples: # takes care of p-c when this mode is the other child
          edge.set_color(highlight_color)
          edge.set_alpha(highlight_alpha)
          edge.set_zorder(highlight_zorder)
        else:
          edge.set_color(background_color)
          edge.set_alpha(background_alpha)
          edge.set_zorder(background_zorder)
      for (c1,c2), [N_k, [edge]] in cc.items():
        if (c1 == modeNo) or (c2 == modeNo): # takes care of c-c edges
          edge.set_color(highlight_color)
          edge.set_alpha(highlight_alpha)
          edge.set_zorder(highlight_zorder)
        else:
          edge.set_color(background_color)
          edge.set_alpha(background_alpha)
          edge.set_zorder(background_zorder)

    else: # this cannot touch the highlighted mode, so we set everything to background
      for (p,c), [N_k, [edge]] in pc.items():
        edge.set_color(background_color)
        edge.set_alpha(background_alpha)
        edge.set_zorder(background_zorder)
      for (c1,c2), [N_k, [edge]] in cc.items():
        edge.set_color(background_color)
        edge.set_alpha(background_alpha)
        edge.set_zorder(background_zorder)

  # no return statement because we alter things in place

####################################################################################################
#
#
#                          scatters of energies associated with each coupling
#
#
####################################################################################################
def __compute_bounds(x, s="none", f=1.05):
  """ used to compute the lower and upper bounds on axes because the automatic scaling appears to be poor"""
  len_x = len(x)
  if s=="none":
    s = np.zeros((len_x,))

  xmin = min([x[i] - s[i] for i in range(len_x)])
  xmax = max([x[i] + s[i] for i in range(len_x)])

  S = (xmax+xmin)/2.
  D = (xmax-xmin)/2.

  return S-f*D, S+f*D

#########################
def __compute_bins(x, n=50, log=False):
  minx, maxx = __compute_bounds(x)
  if log:
    minx = max(minx, 0)
    return np.logspace(max(np.log10(minx),-100), np.log10(maxx), n+1)
  else:
    return np.linspace(minx, maxx, n+1)

##################################################
def Hns_coup_detuning(Hns_coup, system):
  """ 
  generates a scatter plot of the energy assoicated with each coupling as a function of the total detuning
    \Delta_a + \Delta_b + \Delta_c = \omega_a + \omega_b + \omega_c

  requires Hns_coup to have the form returned by nmode_state.Hns_coup()
  """
  network = system.network

  ### generate detunings, mean and stdv of Hns_coup
  mH = []
  sH = []
  detunings = []
  for (i,j,k,kappa), h_ns_coup in Hns_coup.items():
    mh = nms.sample_mean(h_ns_coup)
    mH.append( mh )
    sH.append( nms.sample_var(h_ns_coup, xo=mh)**0.5 )
    detunings.append( (network.modes[i].w + network.modes[j].w + network.modes[k].w)/system.Oorb )

  mH = abs(np.array(mH))
  sH = abs(np.array(sH))
  detunings = np.array(detunings)

  fig = plt.figure()
  ax = fig.add_axes([0.15, 0.15, 0.6, 0.6])
  ax_Hns = fig.add_axes([0.775, 0.15, 0.175, 0.6])
  ax_detuning = fig.add_axes([0.15, 0.775, 0.6, 0.175])


  ax.semilogy(detunings, mH, marker="o", markersize=2, markerfacecolor="none", markeredgecolor="b", linestyle="none")
  for det, m, s in zip(detunings, mH, sH):
    ax.plot([det,det], [m-s, m+s], color='b', alpha=0.2)

  ax_Hns.hist(mH, __compute_bins(mH, log=True), histtype="step", orientation="horizontal", log=True)
  ax_Hns.set_yscale('log')
  ax_detuning.hist(detunings, __compute_bins(detunings), histtype="step", orientation="vertical", log=True)

  xmin, xmax = __compute_bounds(detunings)
  ax.set_xlim(xmin=xmin, xmax=xmax)
  ax_detuning.set_xlim(ax.get_xlim())

  ymin, ymax = __compute_bounds(mH, s=sH)
  ymin = max(ymin, 0)
  ax.set_ylim(ymin=ymin, ymax=ymax)
  ax_Hns.set_ylim(ax.get_ylim())

  ax.set_xlabel(r"$(\omega_a + \omega_b + \omega_c)/\Omega_{\mathrm{orb}}$")
  ax.set_ylabel(r"$|H_{abc}|$")

  ax_Hns.set_xlabel(r"count")
  plt.setp(ax_Hns.get_yticklabels(), visible=False)
  ax_detuning.set_ylabel(r"count")
  plt.setp(ax_detuning.get_xticklabels(), visible=False)

  return fig, ax, ax_Hns, ax_detuning

##################################################
def Hns_coup_heuristic(Hns_coup, system):
  """ 
  generates a scatter plot of the energy assoicated with each coupling as a function of the heuristic ranking scheme
    h = (yb + yc)**2 + (Db + Dc)**2 = (yb + yc)**2 + (Oa + wb + wc)**2  <== delegated to mode_selection.compute_heuristic()

  requires Hns_coup to have the form returned by nmode_state.Hns_coup()
  """
  network = system.network

  ### find three mode frequencies
  freqs = system.compute_3mode_freqs()
#  freqs.sort(key=lambda l: l[1])
#  freqs = [w for w, modeNo in freqs]

  ### generate detunings, mean and stdv of Hns_coup
  mH = []
  sH = []
  heuristics = []
  for (i,j,k,kappa), h_ns_coup in Hns_coup.items():
    mh = nms.sample_mean(h_ns_coup)
    mH.append( mh )
    sH.append( nms.sample_var(h_ns_coup, xo=mh)**0.5 )
    wi = network.modes[i].w
    wj = network.modes[j].w
    wk = network.modes[k].w
    if (abs(wi) > abs(wj)) and (abs(wi) > abs(wk)):
      heuristics.append( ms.compute_heuristic(freqs[i], wj, wk, network.modes[j].y, network.modes[k].y) )
    elif (abs(wj) > abs(wk)):
      heuristics.append( ms.compute_heuristic(freqs[j], wi, wk, network.modes[i].y, network.modes[k].y) )
    else:
      heuristics.append( ms.compute_heuristic(freqs[k], wi, wj, network.modes[i].y, network.modes[k].y) )
  
  mH = abs(np.array(mH))
 
  fig = plt.figure()
  ax = fig.add_axes([0.15, 0.15, 0.6, 0.6])
  ax_Hns = fig.add_axes([0.775, 0.15, 0.175, 0.6])
  ax_heuristic = fig.add_axes([0.15, 0.775, 0.6, 0.175])

  ax.semilogy(heuristics, mH, marker="o", markersize=2, markerfacecolor="none", markeredgecolor="b", linestyle="none")
  for heu, m, s in zip(heuristics, mH, sH):
    ax.plot([heu,heu], [m-s, m+s], color='b', alpha=0.2)

  ax_Hns.hist(mH, __compute_bins(mH, log=True), histtype="step", orientation="horizontal", log=True)
  ax_Hns.set_yscale('log')
  ax_heuristic.hist(heuristics, __compute_bins(heuristics), histtype="step", orientation="vertical", log=True)

  xmin, xmax = __compute_bounds(heuristics)
  ax.set_xlim(xmin=xmin, xmax=xmax)
  ax_heuristic.set_xlim(xmin=xmin, xmax=xmax)

  ymin, ymax = __compute_bounds(mH, s=sH)
  ymin = max(ymin, 0)
  ax.set_ylim(ymin=ymin, ymax=ymax)
  ax_Hns.set_ylim(ax.get_ylim())


  ax.set_xlabel(r"$(\gamma_b + \gamma_c)^{2} + (\Omega_a + \omega_b +\omega_c)^{2}$")
  ax.set_ylabel(r"$|H_{abc}|$")

  ax_Hns.set_xlabel(r"count")
  plt.setp(ax_Hns.get_yticklabels(), visible=False)
  ax_heuristic.set_ylabel(r"count")
  plt.setp(ax_heuristic.get_xticklabels(), visible=False)

  return fig, ax, ax_Hns, ax_heuristic

##################################################
def Hns_coup_Ethr(Hns_coup, system):
  """ 
  generates a scatter plot of the energy assoicated with each coupling as a function of the Ethr ranking scheme
    Ethr = yb*yc/(4*k**2*wb*wc)*(1 + (Oa+wb+wc)**2/(yb+yc)**2)  <== delegated to mode_selection.compute_Ethr()

  requires Hns_coup to have the form returned by nmode_state.Hns_coup()
  """
  network = system.network

  ### find three mode frequencies
  freqs = system.compute_3mode_freqs()
#  freqs.sort(key=lambda l: l[1])
#  freqs = [w for w, modeNo in freqs]

  ### generate detunings, mean and stdv of Hns_coup
  mH = []
  sH = []
  Ethrs = []
  for (i,j,k,kappa), h_ns_coup in Hns_coup.items():
    mh = nms.sample_mean(h_ns_coup)
    mH.append( mh )
    sH.append( nms.sample_var(h_ns_coup, xo=mh)**0.5 )
    wi = network.modes[i].w
    wj = network.modes[j].w
    wk = network.modes[k].w
    for (a,b,kappa) in network.K[i]: # find correct coupling coeff
      if ((j == a) and (k == b)) or ((j == b) and (k == a)):
        break
    else:
      raise ValueError("could not locate this coupling: (%d,%d,%d) in nmode_diagnostic.Hns_coup_Ethr()" % (i,j,k))

    if (abs(wi) > abs(wj)) and (abs(wi) > abs(wk)):
      Ethrs.append( ms.compute_Ethr(freqs[i], wj, wk, network.modes[j].y, network.modes[k].y, kappa) )
    elif (abs(wj) > abs(wk)):
      Ethrs.append( ms.compute_Ethr(freqs[j], wi, wk, network.modes[i].y, network.modes[k].y, kappa) )
    else:
      Ethrs.append( ms.compute_Ethr(freqs[k], wi, wj, network.modes[i].y, network.modes[k].y, kappa) )

  mH = abs(np.array(mH))

  fig = plt.figure()
  ax = fig.add_axes([0.15, 0.15, 0.6, 0.6])
  ax_Hns = fig.add_axes([0.775, 0.15, 0.175, 0.6])
  ax_Ethr = fig.add_axes([0.15, 0.775, 0.6, 0.175])

  ax.semilogy(Ethrs, mH, marker="o", markersize=2, markerfacecolor="none", markeredgecolor="b", linestyle="none")
  for Eth, m, s in zip(Ethrs, mH, sH):
    ax.plot([Eth,Eth], [m-s, m+s], color='b', alpha=0.2)

  ax_Hns.hist(mH, __compute_bins(mH, log=True), histtype="step", orientation="horizontal", log=True)
  ax_Hns.set_yscale('log')
  ax_Ethr.hist(Ethrs, __compute_bins(Ethrs), histtype="step", orientation="vertical", log=True)

  xmin, xmax = __compute_bounds(Ethrs)
  ax.set_xlim(xmin=xmin, xmax=xmax)
  ax_Ethr.set_xlim(xmin=xmin, xmax=xmax)

  ymin, ymax = __compute_bounds(mH, s=sH)
  ymin = max(ymin, 0)
  ax.set_ylim(ymin=ymin, ymax=ymax)
  ax_Hns.set_ylim(ax.get_ylim())

  ax.set_xlabel(r"$E_{thr}$")
  ax.set_ylabel(r"$|H_{abc}|$")

  ax_Hns.set_xlabel(r"count")
  plt.setp(ax_Hns.get_yticklabels(), visible=False)
  ax_Ethr.set_ylabel(r"count")
  plt.setp(ax_Ethr.get_xticklabels(), visible=False)

  return fig, ax, ax_Hns, ax_Ethr

####################################################################################################
#
#
#                             probability distribution (sampled over time)
#
#
####################################################################################################
def E_distrib(q, minE=1e-40, maxE=1e-10, num_bins=50, log=False, log_bins=False, n_l_m=False, mode_nums=False):
  """
  computes the distributions of mode energies over time. From this distribution, you should be able to obtain mean{A**2}, etc
    essentially just a histogram of the energies, treating each sample as independent
  """
  fig = plt.figure()
  ax = plt.subplot(1,1,1)

  weights = np.ones((len(q[0]),))/len(q[0])

  if log_bins:
    bins = np.logspace(np.log10(minE), np.log10(maxE), num_bins+1)
  else:
    bins = np.linspace(minE, maxE, num_bins+1)

  for m, Q in enumerate(q):
    if mode_nums and (m not in mode_nums):
      continue
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])

#    ax.hist([ l1**2 + l2**2 for l1,l2 in a], bins, histtype="step", normed=True, log=log, label=label)
    try:
      ax.hist([ l1**2 + l2**2 for l1,l2 in Q], bins=bins, histtype="step", weights=weights, log=log, label=label)
    except ValueError:
      print "WARNING: problem with E_distrib histogram binning. Reverting to automatic bin spacing"
      ax.hist([ l1**2 + l2**2 for l1,l2 in Q], num_bins, histtype="step", weights=weights, log=log, label=label)
      

  if log_bins:
    ax.set_xscale('log')

  return fig, ax

##################################################
def H_distrib(H, minE=1e-40, maxE=1e-10, num_bins=50, log=False, log_bins=False):
  """
  computes the distributions of the Neutron star Hamiltonian over time. In the absence of damping and external forcing (coupling to tidal potential), this number is conserved.
    essentially just a histogram of the Hamiltonian, treating each sample as independent

  expects H to a list of real numbers
  """
  fig = plt.figure()
  ax = plt.subplot(1,1,1)

  if log_bins:
    bins = np.logspace(np.log10(minE), np.log10(maxE), num_bins+1)
  else:
    bins = np.linspace(minE, maxE, num_bins+1)

#  ax.hist(H, bins, histtype="step", normed=True, log=log, label=r"$H_{\mathrm{NS}}$")
  try:
    ax.hist(H, bins, histtype="step", weights=np.ones((len(H),))/len(H), log=log, label=r"$H_{\mathrm{NS}}$")
  except ValueError:
    print "WARNING: problem with H_distrib histogram binning. Reverting to automatic bin spacing"
    ax.hist(H, num_bins, histtype="step", weights=np.ones((len(H),))/len(H), log=log, label=r"$H_{\mathrm{NS}}$")

  if log_bins:
    ax.set_xscale('log')

  return fig, ax


##################################################
#
#  multi-gen distributions through time
#
##################################################
def multi_gen_E_distrib(q, gens, minE=1e-40, maxE=1e-10, num_bins=50, log=False, log_bins=False, n_l_m=False, mode_nums=False):
  """
  computes the distributions of mode energies over time. From this distribution, you should be able to obtain mean{A**2}, etc
    essentially just a histogram of the energies, treating each sample as independent

  breaks the modes into generations defined by "gens"
  """
  fig = plt.figure()
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  weights=np.ones((len(q[0]),))/len(q[0])

  if log_bins:
    bins = np.logspace(np.log10(minE), np.log10(maxE), num_bins+1)
  else:
    bins = np.linspace(minE, maxE, num_bins+1)

  axs = []
  for genNo, gen in enumerate(gens):
    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    for m in gen:
      if mode_nums and (m not in mode_nums):
        continue
      label = "mode %d" % m
      if n_l_m:
        label += ":%s" % str(n_l_m[m])
#      ax.hist([ l1**2 + l2**2 for l1,l2 in q[m]], bins, histtype="step", normed=True, log=log, label=label)
      try:
        ax.hist([ l1**2 + l2**2 for l1,l2 in q[m]], bins, histtype="step", weights=weights, log=log, label=label)
      except:
        ax.hist([ l1**2 + l2**2 for l1,l2 in q[m]], num_bins, histtype="step", weights=weights, log=log, label=label)
    if log_bins:
      ax.set_xscale("log")

    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs

####################################################################################################
#
#
#                                      conc diagram
#
#
####################################################################################################
def conc(xs, network, annotate=0, annotate_fontsize=8, fitparams=False, fitrange=False):
  """ 
  computes a conc diagram
    xs must be of the form [x0, x1, x2, x3, ...]
  xs will then be sorted and plotted as a function of number of modes

  annotate is the number of modes to explicitly label on the plot
    assumes xs is ordered by modeNo

  if fitparams: plots the fit using nmode_state.broken_PowLaw()
  """
  s = sum(xs)
  l = len(xs)
  modeNos = np.arange(0,l+1,1)
  ys = [0]
  y = 0
  sorted_xs = [(x, modeNo) for modeNo, x in enumerate(xs)]
  sorted_xs.sort(key=lambda l: l[0], reverse=True)
  ys = ys + [x for x, modeNo in sorted_xs]
  y = sum(ys)
  ys = np.array(ys)/y
  
  fig = plt.figure()
  ax_ul = plt.subplot(2,2,1) # upper left ==> lin-lin
#  ax_ur = plt.subplot(2,2,2) # upper right
  ax_ll = plt.subplot(2,2,3) # lower left ==> lin-log
  ax_lr = plt.subplot(2,2,4) # lower right => log-log

  plt.subplots_adjust(hspace=0.05, wspace=0.05)

  ax_ul.plot(modeNos, ys, '-b')  
  ax_ll.semilogy(modeNos, ys, '-b')
  ax_lr.loglog(modeNos, ys, '-b')

  ax_ul.set_xlim(xmin=0, xmax=max(modeNos))
  ax_ll.set_xlim(xmin=0, xmax=max(modeNos))
  ax_lr.set_xlim(xmin=1, xmax=max(modeNos))

  if fitparams:
    if not fitrange:
      fitrange = [1,-1]
    modeNos = modeNos[fitrange[0]:fitrange[1]]
    fit = nms.broken_PowLaw(modeNos, fitparams)
    ax_ul.plot(modeNos, fit, '-r', label='fit', alpha=0.5)
    ax_ll.semilogy(modeNos, fit, '-r', label='fit', alpha=0.5)
    ax_lr.loglog(modeNos, fit, '-r', label='fit', alpha=0.5)

  plt.setp(ax_ul.get_xticklabels(), visible=False)
  plt.setp(ax_lr.get_yticklabels(), visible=False)

  for ind in range(min(annotate, l)):
    _y = ys[ind+1]
    modeNo = sorted_xs[ind][1]
    ax_ul.text(ind+1, _y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)
    ax_ll.text(ind+1, _y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)
    ax_lr.text(ind+1, _y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)

  return fig, [ax_ul, ax_ll, ax_lr]

##################################################
def Cconc(xs, network, annotate=0, annotate_fontsize=8, fitparams=False, fitrange=False):
  """ 
  computes a conc diagram
    xs must be of the form [x0, x1, x2, x3, ...]
  xs will then be sorted and plotted as a function of number of modes

  annotate is the number of modes to explicitly label on the plot
    assumes xs is ordered by modeNo
  """
  s = sum(xs)
  l = len(xs)
  modeNos = np.arange(0,l+1,1)
  ys = [0]
  y = 0
  sorted_xs = [(x, modeNo) for modeNo, x in enumerate(xs)]
  sorted_xs.sort(key=lambda l: l[0], reverse=True)
  for x, modeNo in sorted_xs:
    y += x
    ys.append(y)
  ys = np.array(ys)/y

  fig = plt.figure()
  ax_ul = plt.subplot(2,2,1) # upper left ==> lin-lin
#  ax_ur = plt.subplot(2,2,2) # upper right
  ax_ll = plt.subplot(2,2,3) # lower left ==> lin-log
  ax_lr = plt.subplot(2,2,4) # lower right => log-log

  plt.subplots_adjust(hspace=0.05, wspace=0.05)

  ax_ul.plot(modeNos, ys, '-b')
  ax_ll.semilogy(modeNos, ys, '-b')
  ax_lr.loglog(modeNos, ys, '-b')

  ax_ul.plot(modeNos, np.linspace(0,1,l+1), 'k', alpha=0.2)
  ax_ll.plot(modeNos, np.linspace(0,1,l+1), 'k', alpha=0.2)
  ax_lr.plot(modeNos, np.linspace(0,1,l+1), 'k', alpha=0.2)

  ax_ul.set_xlim(xmin=0, xmax=max(modeNos))
  ax_ul.set_ylim(ymin=0, ymax=1)
  ax_ll.set_xlim(xmin=0, xmax=max(modeNos))
  ax_ll.set_ylim(ymax=1)
  ax_lr.set_xlim(xmin=1, xmax=max(modeNos))
  ax_lr.set_ylim(ymax=1)

  if fitparams:
    if not fitrange:
      fitrange = [1,-1]
    modeNos = modeNos[fitrange[0]:fitrange[1]]
    fit = []
    f = 0
    for x in nms.broken_PowLaw(modeNos, fitparams): # compute cumulative curve
      f += x
      fit.append(f)

    ax_ul.plot(modeNos, fit, '-r', label='fit', alpha=0.5)
    ax_ll.semilogy(modeNos, fit, '-r', label='fit', alpha=0.5)
    ax_lr.loglog(modeNos, fit, '-r', label='fit', alpha=0.5)

  plt.setp(ax_ul.get_xticklabels(), visible=False)
  plt.setp(ax_lr.get_yticklabels(), visible=False)

  for ind in range(min(annotate, l)):
    _y = ys[ind+1]
    modeNo = sorted_xs[ind][1]
    ax_ul.text(ind+1, _y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)
    ax_ll.text(ind+1, _y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)
    ax_lr.text(ind+1, _y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)

  return fig, [ax_ul, ax_ll, ax_lr]

'''
##################################################
def multi_gen_Cconc(xs, gens, network, annotate=0, annotate_fontsize=8, log=False, cumulative=True):
  """
  computes a conc diagram 
    xs must be of the form [x0, x1, x2, x3, ...]
    gens must have the form [[modeNo00, modeNo01, ...], [modeNo10, modeNo11, ...], ...]
  """
  fig = plt.figure()
  axs = []
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  for genNo, gen in enumerate(gens):
    ys = [0]
    y = 0
    this_sorted_xs = [(xs[modeNo], modeNo) for modeNo in gen]
    this_sorted_xs.sort(key=lambda l: l[0], reverse=True)
    for x, modeNo in this_sorted_xs:
      y += x
      ys.append( y )

    modeNos = np.linspace(0,1,len(gen)+1)

    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    if log:
      ax.loglog(modeNos, modeNos, 'k', alpha=0.2)
      ax.loglog(modeNos, np.array(ys)/y, '-b')
    else:
      ax.plot(modeNos, modeNos, 'k', alpha=0.2)
      ax.plot(modeNos, np.array(ys)/y, '-b')

    for ind in range(min(annotate, len(gen))):
      _y = ys[ind+1]
      modeNo = this_sorted_xs[ind][1]
      ax.text(1.0*(ind+1)/len(gen), _y/y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)

    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs

##################################################
def multi_gen_conc(xs, gens, network, annotate=0, annotate_fontsize=8, log=False, cumulative=True):
  """
  computes a conc diagram 
    xs must be of the form [x0, x1, x2, x3, ...]
    gens must have the form [[modeNo00, modeNo01, ...], [modeNo10, modeNo11, ...], ...]
  """
  fig = plt.figure()
  axs = []
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  for genNo, gen in enumerate(gens):
    ys = [0]
    y = 0
    this_sorted_xs = [(xs[modeNo], modeNo) for modeNo in gen]
    this_sorted_xs.sort(key=lambda l: l[0], reverse=True)
    ys = ys + [x for x, modeNo in this_sorted_xs]
    y = sum(ys)

    modeNos = np.linspace(0,1,len(gen)+1)

    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    if log:
      ax.loglog(modeNos, np.array(ys)/y, '-b')
    else:
      ax.plot(modeNos, np.array(ys)/y, '-b')

    for ind in range(min(annotate, len(gen))):
      _y = ys[ind+1]
      modeNo = this_sorted_xs[ind][1]
      ax.text(1.0*(ind+1)/len(gen), _y/y, r"$\leftarrow$(%d,%d,%d)" % network.modes[modeNo].get_nlm(), ha='left', va='center', fontsize=annotate_fontsize)

    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs
'''
