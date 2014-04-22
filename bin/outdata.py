#!python_alias
usage = "written to generate a large table with all pertinant integration data in it"

from optparse import OptionParser
import nmode_utils as nmu
from numpy import infty

####################################################################################################

parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--cachefilename", default=False, type="string", help="the filename of a cachefile containing paths to all *.out files of interest")
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
  if opts.verbose: print "reading in all *.out filenames from "+os.getcwd()
  filenames = sorted( glob.glob("./*.out") )

dat = []

### read data from outs
for filename in filenames:
  if opts.verbose: print "reading integration data from "+filename

  try:
    ### manipulate file directly to save on memory
    f = open(filename, 'r')
    # find out how many modes there are
    for line in f:
      if line[0] != '#':
        line = [ l for l in line.strip().split()]
        N_m = (len(line) - 1) /2 # number of modes
        break
    else:
      N_m = -1
    f.close()  

    # pull out data
    f = open(filename, 'r')
    line_No = 0
    start = step = stop = 0.0
    for line in f:
      if line[0] != '#':
        line = line.strip().split()
        if len(line) == 1+2*N_m:
          if line_No == 0:
            start = float(line[0])
          elif line_No == 1:
            step = float(line[0])-start
          stop = float(line[0])

          line_No += 1
    f.close()

    del f

  except:
    print "\tERROR when reading "+filename
    continue

  dat.append( (filename, N_m, start, stop, step) )

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
filename & No. modes & $t/P_{\mathrm{orb}}[0]$ & $\mathrm{d}t/P_{\mathrm{orb}}$ & $t/P_{\mathrm{orb}}[-1]$ \\ 
\cline{1-%d}
"""

no_sum_cols = 5

sum_tablename = opts.output_dir+"/out-sum%s.tex" % opts.tag
if opts.verbose: print "writing summary tables to:\n\t%s" % (sum_tablename)

### a table containing general info about the networks
sum_tablefile = open(sum_tablename, "w")
print >> sum_tablefile, ldoc_preamble
if opts.tag != "":
  print >> sum_tablefile, tab_preamble % ("General Integration Parameters: "+opts.tag[1:].replace("_","\_"), no_sum_cols)
else:
  print >> sum_tablefile, tab_preamble % ("General Integration Parameters", no_sum_cols)
print >> sum_tablefile, sum_headings % (no_sum_cols, no_sum_cols)
sum_row_No = 0

for filename, N_m, start, stop, step in dat:
  if opts.verbose: print "working on "+filename

  ### add row to sum_tablefile
  if (sum_row_No % 3) == 0:
    print >> sum_tablefile, r"\cline{1-%d}" % no_sum_cols
  ###                   filename  nM  start  step   stop
  print >> sum_tablefile, r"%s & %d & %.1f & %.1f & %.1f \\" % (filename.replace("_","\_"), N_m, start, step, stop) 
  sum_row_No += 1

  if (sum_row_No+1)%45 == 0:
    print >> sum_tablefile, tab_suffix
    if opts.tag != "":
      print >> sum_tablefile, tab_preamble % ("General Integration Parameters: "+opts.tag[1:].replace("_","\_"), no_sum_cols)
    else:
      print >> sum_tablefile, tab_preamble % ("General Integration Parameters", no_sum_cols)
    print >> sum_tablefile, sum_headings % (no_sum_cols, no_sum_cols)

### end sum_table
print >> sum_tablefile, tab_suffix
print >> sum_tablefile, doc_suffix
sum_tablefile.close()

