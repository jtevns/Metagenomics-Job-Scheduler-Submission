#given a contig name map and binlist, rename the contigs
#to the new name given by anvio

import sys
import pandas

contig_map = sys.argv[1]
bin_list = sys.argv[2]
outfile = sys.argv[3]

map_table = pandas.read_csv(contig_map, sep="\t", header=None, index_col=1)
bin_table = pandas.read_csv(bin_list,sep=",",header=None, index_col=0)

joined = bin_table.join(map_table,how="left",lsuffix="bin_list", rsuffix="contig_map")
joined[[0,1]].to_csv(outfile,sep="\t",header=False,index=False)