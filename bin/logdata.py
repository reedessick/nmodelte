#!	/home/ressick/local/bin/Python-2.7.3/python
usage = "written to generate a large table with all pertinant log data in it"

from optparse import OptionParser
import nmode_utils as nmu
from numpy import infty

####################################################################################################

parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--cachefilename", default=False, type="string", help="the filename of a cachefile containing paths to all *.log files of interest")
parser.add_option("-o", "--output-dir", default="./", type="string", help="the directory into which the resulting table will be written. DEFAULT=\"./\"")

parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if opts.tag != "":
  opts.tag = "."+opts.tag

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
  if opts.verbose: print "reading in all *.log filenames from "+os.getcwd()
  filenames = sorted( glob.glob("./*.log") )

dat = []

### read data from logs
for filename in filenames:
  if opts.verbose: print "reading system from "+filename

  try:
    sys = nmu.load_log(filename)
  except:
    print "\tERROR when reading "+filename
    continue

  net = sys.network
  
  ### pull out data
  gens, coups = net.gens()

  this_log = [ filename, (sys.Mprim, sys.Rprim, sys.Mcomp, sys.eccentricity, sys.Porb), (len(net), len(net.to_triples())) ]

  these_gens = []

  for genNo, gen in enumerate(gens):
    max_m = -infty
    min_m =  infty
    max_l = -infty
    min_l =  infty
    max_w = -infty
    min_w =  infty

    for modeNo in gen:
      _, l, m = sys.network.nlm[modeNo]
      w, _, _ = sys.network.wyU[modeNo]

      max_l = max(max_l, l)
      min_l = min(min_l, l)
      max_m = max(max_m, m)
      min_m = min(min_m, m)
      max_w = max(max_w, abs(w))
      min_w = min(min_w, abs(w))

      if genNo < len(gens)-1:
        no_coups = len(coups[genNo])
      else:
        no_coups = 0

    these_gens.append( (len(gen), no_coups, max_l, min_l, max_m, min_m, max_w/sys.Oorb, min_w/sys.Oorb) )

  this_log.append( these_gens )

  dat.append( this_log )

### write table >> file
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
tab_preamble=r"""
\begin{table}
\caption{%s}
\begin{center}
\begin{tabular}{*{%d}{c}}
"""

tab_suffix = r"""
\end{tabular}
\end{center}
\end{table}
"""

### column headings
sum_headings = r"""
\cline{1-%d}
filename & $M_{\mathrm{prim}}/M_{\odot}$ & $R_{\mathrm{prim}}/R_{\odot}$ & $M_{\mathrm{comp}}/M_{J}$ & $\epsilon$ & $P_{\mathrm{orb}}$ [sec] & No. modes & No. triples & No. gens \\ 
\cline{1-%d}
"""
no_sum_cols = 9

spc_headings = r"""
\cline{2-%d}
& \multirow{2}{*}{No. modes} & \multirow{2}{*}{No. coups} & \multicolumn{2}{c}{$l$} & \multicolumn{2}{c}{$m$} & \multicolumn{2}{c}{$\omega/\Omega_{\mathrm{orb}}$} \\ 
& \multirow{2}{*}{}          & \multirow{2}{*}{}          &        min & max        &         min & max       &                 min & max                          \\
\cline{1-%d}
"""
no_spc_cols = 9

sum_tablename = opts.output_dir+"/log-sum%s.tex" % opts.tag
spc_tablename = opts.output_dir+"/log-spc%s.tex" % opts.tag
if opts.verbose: print "writing summary tables to:\n\t%s\n\t%s" % (sum_tablename, spc_tablename)

### a table containing general info about the networks
sum_tablefile = open(sum_tablename, "w")
print >> sum_tablefile, ldoc_preamble
if opts.tag != "":
  print >> sum_tablefile, tab_preamble % ("General System Properties: "+opts.tag[1:].replace("_","\_"), no_sum_cols)
else:
  print >> sum_tablefile, tab_preamble % ("General System Properties", no_sum_cols)
print >> sum_tablefile, sum_headings % (no_sum_cols, no_sum_cols)
sum_row_No = 0

### a series of tables specifically for each network
spc_tablefile = open(spc_tablename, "w")
print >> spc_tablefile, pdoc_preamble

for filename, (Mprim, Rprim, Mcomp, ecc, Porb), (no_modes, no_triples), these_gens in dat:
  if opts.verbose: print "working on "+filename

  ### add row to sum_tablefile
  if (sum_row_No % 3) == 0:
    print >> sum_tablefile, r"\cline{1-%d}" % no_sum_cols
  ###                   filename  Mp     Rp     Mc     e      P     nM   nT   nG
  print >> sum_tablefile, r"%s & %.3f & %.3f & %.3f & %.3f & %.1f & %d & %d & %d \\" % (filename.replace("_","\_"), Mprim, Rprim, Mcomp, ecc, Porb, no_modes, no_triples, len(these_gens)) 
  sum_row_No += 1

  ### generate new spc_table
  print >> spc_tablefile, tab_preamble % (filename.replace("_","\_"), no_spc_cols)
  print >> spc_tablefile, spc_headings % (no_spc_cols, no_spc_cols)

  ### add rows to spc_table
  for genNo, (no_modes, no_coups, max_l, min_l, max_m, min_m, max_w, min_w) in enumerate(these_gens):
    if (genNo % 3) == 0:
      print >> spc_tablefile, r"\cline{1-%d}" % no_spc_cols
    ###                        genNo   nM   nT  minL maxL minM maxM  minW   maxW
    print >> spc_tablefile, r"gen-%d & %d & %d & %d & %d & %d & %d & %.6f & %.6f \\"  % (genNo, no_modes, no_coups, min_l, max_l, min_m, max_m, min_w, max_w)
  ### end spc_table
  print >> spc_tablefile, tab_suffix

  if (sum_row_No+1)%45 == 0:
    print >> sum_tablefile, tab_suffix
    if opts.tag != "":
      print >> sum_tablefile, tab_preamble % ("General System Properties: "+opts.tag[1:].replace("_","\_"), no_sum_cols)
    else:
      print >> sum_tablefile, tab_preamble % ("General System Properties", no_sum_cols)
    print >> sum_tablefile, sum_headings % (no_sum_cols, no_sum_cols)


### end sum_table
print >> sum_tablefile, tab_suffix
print >> sum_tablefile, doc_suffix
sum_tablefile.close()

### end spc_table
print >> spc_tablefile, doc_suffix
spc_tablefile.close()


