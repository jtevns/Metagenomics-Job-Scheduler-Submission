#given a contig name map and binlist, rename the contigs
#to the new name given by anvio

import sys
import pandas

contig_map = sys.argv[1]
bam_header = sys.argv[2]
outfile = bam_header.replace(".sam","-anvio_names.sam")

name_dict = {}
with open(contig_map) as f:
    for line in f:
        value,key,*rest = line.split()
        name_dict[key] = value

output_list = []  
with open(bam_header) as f:
    for line in f:
        old_id = line.split("\t")[1].split(":")[1]
        if old_id in name_dict.keys():
            output_list.append(line.replace(old_id,name_dict[old_id])) 
        else:
            output_list.append(line)

output = "".join(output_list)

f = open(outfile, "w")
f.write(output)
f.close()

