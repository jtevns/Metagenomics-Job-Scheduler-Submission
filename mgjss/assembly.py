import os
import click
from yaml import load, dump, FullLoader
from mgjss.slurm import SlurmJob
import pandas as pds

def get_fastq_paths(table, keys):
    fq1_paths = []
    fq2_paths= []
    data = pds.read_csv(table, header=0)

    for key in keys:
        path_1 = data.loc[data["sample"] == key, "qced_fq1"].item()
        path_2 = data.loc[data["sample"] == key, "qced_fq2"].item()
        fq1_paths.append(path_1)
        fq2_paths.append(path_2)

    return(fq1_paths, fq2_paths)


# cli for inputs
@click.command()
@click.argument('fastqs', type=click.Path(exists=True, readable=True))
@click.argument('assembly_scheme', type=click.Path(exists=True, readable=True))
@click.option('--overwrite' , help='force overwrite files if needed', default=False)
@click.argument('output_dir', type=click.Path(exists=False, writable=True))
@click.option('--account' , help='choose account to submit job to', required=True)
def assemble(fastqs,assembly_scheme,overwrite,output_dir,account):
    #make output directory
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    #parse binning_scheme
    with open(assembly_scheme) as f:
        scheme = load(f, Loader=FullLoader)

    # for each sample create the binning jobs 
    for assembly in scheme.keys():
       #make sample dir in output dir
        if not os.path.exists(output_dir + "/" + assembly):os.makedirs(output_dir + "/" + assembly)

        # get paths to assembly and bams
        fqs = get_fastq_paths(fastqs, scheme[assembly])
        r1 = ",".join(fqs[0])
        r2 = ",".join(fqs[1])


        #create and submit jobs
        # 1 assemble
        megahit_job = SlurmJob(
            name="megahit-" + str(assembly),
            account=account,
            cores=36,
            mem="180gb",
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/Megahit/megahit.sif megahit -1 {} -2 {} -t 36 --presets meta-sensitive -o {}/{}/Megahit_meta-sensitive_out/".format(r1,r2,output_dir,assembly)
        )
        megahit_id = megahit_job.submit()
        print(megahit_id)
        

        # 2 assembly stats
        assembly_stats_job = SlurmJob(
            name="assembly_stats-" + str(assembly),
            account=account,
            dependency=megahit_id,
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/bbtools/bbtools.sif stats.sh {0}/{1}/Megahit_meta-sensitive_out/final.contigs.fa > {0}/{1}/Megahit_meta-sensitive_out/assembly_stats.txt".format(output_dir, assembly)
        )
        assembly_stats_jobid = assembly_stats_job.submit()
        print(assembly_stats_jobid)
