import os
import click
from yaml import load, dump, FullLoader
from mgjss.slurm import SlurmJob
import pandas as pds

def get_bam_paths(table,ref,keys):
    bam_paths = []
    data = pds.read_csv(table, header=0)

    for key in keys:
        path = data.loc[(data["ref"] == ref) & (data["sample"] == key), "path"].item()
        bam_paths.append(path)

    return(",".join(bam_paths))

def get_assem_path(table, key):
    data = pds.read_csv(table, header=0)
    ref_path = data.loc[data["assembly"] == key, "path"].item()
    return(ref_path)

# cli for inputs
@click.command()
@click.argument('assemblies', type=click.Path(exists=True, readable=True))
@click.argument('bams', type=click.Path(exists=True, readable=True))
@click.argument('binning_scheme', type=click.Path(exists=True, readable=True))
@click.option('--overwrite' , help='force overwrite files if needed', default=False)
@click.argument('output_dir', type=click.Path(exists=False, writable=True))
@click.option('--account' , help='choose account to submit job to', required=True)
def concoct(assemblies,bams,binning_scheme,overwrite,output_dir,account):
    #make output directory
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    #parse binning_scheme
    with open(binning_scheme) as f:
        scheme = load(f, Loader=FullLoader)

    # for each sample create the binning jobs 
    for assembly in scheme.keys():
       #make sample dir in output dir
        if not os.path.exists(output_dir + "/" + assembly):os.makedirs(output_dir + "/" + assembly)
       #make fasta bins dir
        if not os.path.exists(output_dir + "/" + assembly + "/bins"):os.makedirs(output_dir + "/" + assembly + "/bins")

        # get paths to assembly and bams
        assembly_fasta = get_assem_path(assemblies,assembly)
        bams = get_bam_paths(bams,assembly,scheme[assembly])

        #make output paths
        cut_fasta = output_dir + "/" + assembly + "/contigs_10k.fasta"
        cut_bed = output_dir + "/" + assembly + "/contigs_10k.bed"
        cov_table = output_dir + "/" + assembly + "/coverage_table.tsv"
        bin_dir = output_dir + "/" + assembly + "/bins"
        concoct_output = output_dir + "/" + assembly + "/concoct_output/"
        checkm_dir = output_dir + "/" + assembly + "/checkm"
        binstats = output_dir + "/" + assembly + "/binstats.tsv"

        #create and submit jobs
        # 1 cut up contigs
        cut_contigs_job = SlurmJob(
            name="cut_contigs-" + str(assembly),
            account=account,
            cores=1,
            mem="8gb",
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/concoct/concoct.sif cut_up_fasta.py {} -c 10000 -o 0 --merge_last -b {} > {}".format(assembly_fasta, cut_bed, cut_fasta)
        )
        cut_contigs_jobid = cut_contigs_job.submit()
        print(cut_contigs_jobid)
        

        # 2 gen_cov
        gen_cov_job = SlurmJob(
            name="gen_cov-" + str(assembly),
            account=account,
            cores=1,
            mem="8gb",
            dependency=cut_contigs_jobid,
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/concoct/concoct.sif concoct_coverage_table.py {} {} > {}".format(cut_bed, bams, cov_table)
        )
        gen_cov_jobid = gen_cov_job.submit()
        print(gen_cov_jobid)

        # 3 bin
        concoct_job = SlurmJob(
            name="concoct-" + str(assembly),
            account=account,
            cores=36,
            mem="180gb",
            dependency=gen_cov_jobid,
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/concoct/concoct.sif concoct --composition_file {} --coverage_file {} -b {}".format(cut_fasta, cov_table, concoct_output)
        )
        concoct_jobid = concoct_job.submit()
        print(concoct_jobid)

        # 4 merge contigs
        merge_contigs_job = SlurmJob(
            name="merge_contigs-" + str(assembly),
            account=account,
            cores=1,
            mem="8gb",
            dependency=concoct_jobid,
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/concoct/concoct.sif merge_cutup_clustering.py {}/clustering_gt1000.csv > {}/clustering_merged.csv".format(concoct_output, concoct_output)
        )
        merge_contigs_jobid = merge_contigs_job.submit()
        print(merge_contigs_jobid)

        # 5 extract fastas
        extract_fastas_job = SlurmJob(
            name="extract_fastas-" + str(assembly),
            account=account,
            cores=1,
            mem="8gb",
            dependency=merge_contigs_jobid,
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/concoct/concoct.sif extract_fasta_bins.py {} {}/clustering_merged.csv --output_path {}".format(assembly_fasta, concoct_output, bin_dir)
        )
        extract_fastas_jobid = extract_fastas_job.submit()
        print(extract_fastas_jobid)

        # 6 checkm
        checkm_job = SlurmJob(
            name="checkm-" + str(assembly),
            account=account,
            cores=36,
            mem="180gb",
            dependency=extract_fastas_jobid,
            cmd="singularity exec /nfs/turbo/lsa-dudelabs/containers/checkm/checkm.sif checkm lineage_wf {} {} -x fa --tab_table -t 36 --pplacer_threads 36 -f {}".format(bin_dir, checkm_dir, binstats)
        )
        checkm_jobid = checkm_job.submit()
        print(checkm_jobid)