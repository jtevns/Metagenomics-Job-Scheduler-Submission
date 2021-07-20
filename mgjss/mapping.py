import os
import click
from yaml import load, dump, FullLoader
from mgjss.slurm import SlurmJob
import pandas as pds

def get_fastq_paths(data, key):
    r1 = data.loc[data["sample"] == key, "qced_fq1"].item()
    r2 = data.loc[data["sample"] == key, "qced_fq2"].item()
    return(r1, r2)

def get_ref_path(data,key):
    ref_path = data.loc[data["assembly"] == key, "path"].item()
    return(ref_path)



# cli for inputs
@click.command()
@click.argument('fastqs', type=click.Path(exists=True, readable=True))
@click.argument('assemblies', type=click.Path(exists=True, readable=True))
@click.argument('mapping_scheme', type=click.Path(exists=True, readable=True))
@click.option('--overwrite' , help='force overwrite files if needed', default=False)
@click.argument('output_dir', type=click.Path(exists=False, writable=True))
@click.option('--account' , help='choose account to submit job to', required=True)
def map(fastqs,assemblies,mapping_scheme,overwrite,output_dir,account):
    #make output directory
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    #parse binning_scheme
    with open(mapping_scheme) as f:
        scheme = load(f, Loader=FullLoader)
    
    #open fastq table
    fastq_table = pds.read_csv(fastqs, header=0)

    #open assembly table
    assembly_table = pds.read_csv(assemblies, header=0)

    # for each sample create the binning jobs 
    for ref in scheme.keys():
       #make sample dir in output dir
        if not os.path.exists(output_dir + "/" + ref):os.makedirs(output_dir + "/" + ref)

        # get paths to assembly and bams
        ref_path = get_ref_path(assembly_table, ref)
        # 1 index ref
        index_ref_job = SlurmJob(
            name="index_ref-" + str(ref),
            account=account,
            mem="32gb",
            cmd="""
                singularity exec /nfs/turbo/lsa-dudelabs/containers/bwa/bwa.sif bwa index {}
                """.format(ref_path)
        )
        index_ref_id = index_ref_job.submit()
        print(index_ref_id)

        for sample in scheme[ref]:          
            fqs = get_fastq_paths(fastq_table, sample)
            r1 = fqs[0]
            r2 = fqs[1]
            outbam = sample + "_sorted.bam"
            
            # 2 map
            bwa_job = SlurmJob(
                name="map-" + str(ref) + "-" + sample,
                account=account,
                cores=20,
                mem="64gb",
                dependency=index_ref_id,
                cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/bwa/bwa.sif bwa mem {} {} {} -t 10  | singularity exec /nfs/turbo/lsa-dudelabs/containers/bwa/bwa.sif samtools sort -o {}/{}/{} -@ 10 -".format(ref_path, r1, r2, output_dir, ref, outbam)
            )
            bwa_jobid = bwa_job.submit()
            print(bwa_jobid)


            # 3 index bam
            index_bam_job = SlurmJob(
                name="index_bam-" + ref + "-" +sample,
                account=account,
                dependency=bwa_jobid,
                cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/bwa/bwa.sif samtools index {}/{}/{}".format(output_dir, ref, outbam)
            )
            index_bam_jobid = index_bam_job.submit()
            print(index_bam_jobid)