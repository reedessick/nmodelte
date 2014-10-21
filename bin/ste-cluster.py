#!python_alias
usage = """python cluster.py [--options] stepkl stepkl stepkl...
 an executable that clusters ste files into pickled objects"""

import numpy as np
import pickle

import fitting
import nmode_utils as nmu

from optparse import OptionParser

#=================================================
def cluster_by_key(sorted_list, key, window):

	clusters = []
	cluster = []
	old_val = -np.infty
	for filename, element in sorted_list:
		val = element["system"][key]
		if val-old_val < window:
			cluster.append( (filename, element) )
		else:
			if cluster:
				clusters.append( cluster )
			cluster = [ (filename, element) ]
		old_val = val
	if cluster:
		clusters.append( cluster )

	return clusters
		
#=================================================
parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--Porb_window", default=10000, type="float", help="the clustering window for Porb")
parser.add_option("", "--Mcomp_window", default=10, type="float", help="the clustering window for Mcomp")

parser.add_option("","--unit-system", default="SI", type="string", help="the system of units used in the plot. Currently support either \"SI\" or \"CGS\"")

parser.add_option("-o", "--output-dir", default="./", type="string")
parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

nmu.set_units(system=opts.unit_system) ### set our system of units
energy_unit = nmu.units["energy"]
time_unit = nmu.units["time"]

if opts.tag:
	opts.tag = "_%s"%opts.tag

#=================================================
### build associations
### load data from files
data = [ ]
for filename in args:
	if opts.verbose: print filename
	sdata, mdata = nmu.load_ste(filename)
	data.append( (filename, sdata) )

### cluster data by Porb_window
if opts.verbose: print "sorting by Porb"
data.sort(key=lambda l: l[1]["system"]["Porb"]) ### sort by Porb
clusters = cluster_by_key(data, "Porb", opts.Porb_window) ### cluster

### cluster data by Mcomp_window
if opts.verbose: print "sorting by Mcomp"
tmp_clusters = []
for cluster in clusters:
	cluster.sort(key=lambda l: l[1]["system"]["Mcomp/Mjup"]) ### sort by Mcomp
	tmp_clusters += cluster_by_key(cluster, "Mcomp/Mjup", opts.Mcomp_window) ### cluster
clusters = tmp_clusters

#=================================================
### build cluster objects
if opts.verbose: print "building cluster objects"
clusters = [fitting.Cluster( [c[0] for c in cluster] ) for cluster in clusters]

#=================================================
### write clusters to file
pklname = "%s/clusters%s.pkl"%(opts.output_dir, opts.tag)
if opts.verbose: print "writing cluster objects into %s"%pklname
file_obj = open(pklname, "w")
pickle.dump( clusters, file_obj )
file_obj.close()



