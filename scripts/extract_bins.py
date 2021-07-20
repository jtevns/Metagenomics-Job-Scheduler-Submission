import sys
from Bio import SeqIO
from collections import defaultdict
import os

assem = sys.argv[1]
clustering = sys.argv[2]
out_dir = sys.argv[3]
cluster_dict = defaultdict()
bins = defaultdict(list)

if not os.path.exists(out_dir): os.makedirs(out_dir)


with open(clustering) as f:
    for line in f:
        key,value = line.split(",")
        cluster_dict[key] = value

for seq in SeqIO.parse(assem, "fasta"):
    if seq.id in cluster_dict:
        bin_num = cluster_dict[seq.id].strip()
        bins[bin_num].append(seq)

for key in bins.keys():
    SeqIO.write(bins[key], out_dir + "/" + str(key) + ".fa", "fasta")