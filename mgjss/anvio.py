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

    return(bam_paths)

def get_assem_path(table, key):
    data = pds.read_csv(table, header=0)
    ref_path = data.loc[data["assembly"] == key, "path"].item()
    return(ref_path)

def get_binlist_path(table, key):
    data = pds.read_csv(table, header=0)
    binlist_path = data.loc[data["assembly"] == key, "binlist"].item()
    return(binlist_path)


# cli for inputs
@click.command()
@click.argument('assemblies', type=click.Path(exists=True, readable=True))
@click.argument('bams', type=click.Path(exists=True, readable=True))
@click.argument('bin_lists', type=click.Path(exists=True, readable=True))
@click.argument('anvio_scheme', type=click.Path(exists=True, readable=True))
@click.option('--overwrite' , help='force overwrite files if needed', default=False)
@click.argument('output_dir', type=click.Path(exists=False, writable=True))
@click.option('--account' , help='choose account to submit job to', required=True)
@click.option('--rename_contigs',  help='choose to rename contigs for anvio', is_flag=True)
def anvio(assemblies,bams,bin_lists,anvio_scheme,overwrite,output_dir,account,rename_contigs):
    # make output dirs
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    if not os.path.exists(output_dir + "/contig_dbs"): os.makedirs(output_dir + "/contig_dbs")
    if not os.path.exists(output_dir + "/profile_dbs"): os.makedirs(output_dir + "/profile_dbs")
    if not os.path.exists(output_dir + "/merged_profile_dbs"): os.makedirs(output_dir + "/merged_profile_dbs")

    profile_dbs = output_dir + "/profile_dbs"
    contig_dbs = output_dir + "/contig_dbs"
    merged_profile_dbs = output_dir + "/merged_profile_dbs"
    #parse anvio scheme
    with open(anvio_scheme) as f:
        scheme = load(f, Loader=FullLoader)

    # for each sample create the binning jobs 
    for assembly in scheme.keys():
        #init intermediate files
        bin_list = get_binlist_path(bin_lists,assembly)
        assembly_path = get_assem_path(assemblies,assembly)
        bam_paths = get_bam_paths(bams,assembly,scheme[assembly])

        #init db names
        contigs_db = assembly + "-CONTIGS.db"

        # rename contigs and files to use new names
        if rename_contigs:
            bam_rename_job_list = []
            profile_bams_job_list = []
            renamed_bam_paths = []
            #make contig map dir 
            if not os.path.exists(output_dir + "/contig_maps"): os.makedirs(output_dir + "/contig_maps")
            if not os.path.exists(output_dir + "/renamed_assemblies"): os.makedirs(output_dir + "/renamed_assemblies")
            if not os.path.exists(output_dir + "/renamed_bams"): os.makedirs(output_dir + "/renamed_bams")
            if not os.path.exists(output_dir + "/renamed_binlists"): os.makedirs(output_dir + "/renamed_binlists")

            contig_maps = output_dir + "/contig_maps"
            renamed_assemblies = output_dir + "/renamed_assemblies"
            renamed_bams = output_dir + "/renamed_bams"
            renamed_assembly = renamed_assemblies + "/" + assembly + "-renamed.fa"
            contig_map = contig_maps + "/" + assembly + "-contigs_map.tsv"
            renamed_binlist = output_dir + "/renamed_binlists/" + assembly + "-renamed_binlist.tsv"

            # rename contigs
            rename_contigs_job = SlurmJob(
                name="rename_contigs-" + str(assembly),
                account=account,
                cores=1,
                mem="8gb",
                cmd="anvi-script-reformat-fasta {} -o {} -l 0 --simplify-names -r {}/{}-contigs_map.tsv".format(assembly_path,renamed_assembly,contig_maps,assembly)
            )
            rename_contigs_jobid = rename_contigs_job.submit()
            assembly_path = renamed_assembly
            print(rename_contigs_jobid)

            #rename bam headers
            for bam_path in bam_paths:
                bam = bam_path.split("/")[-1].split(".")[0]
                rename_bam_header_job = SlurmJob(
                name="rename_bam_header-" + str(assembly + "_" + bam),
                account=account,
                cores=20,
                mem="32gb",
                dependency=rename_contigs_jobid,
                cmd="""
                    samtools view -H {0} > {1}/{2}-tmp_header.sam
                    rename_contigs_in_bam_header.py {3} {1}/{2}-tmp_header.sam
                    samtools reheader {1}/{2}-tmp_header-anvio_names.sam {0} > {1}/{4}-anvio_names.bam
                    anvi-init-bam -o {1}/{4}-anvio_names_init.bam -T 20 {1}/{4}-anvio_names.bam
                    """.format(bam_path, renamed_bams, bam, contig_map, bam)
                )
                rename_bam_header_jobid = rename_bam_header_job.submit()
                print(rename_bam_header_jobid)
                bam_rename_job_list.append(rename_bam_header_jobid)
                renamed_bam_paths.append(renamed_bams + "/" + bam + "-anvio_names_init.bam")
            bam_paths = renamed_bam_paths

            # rename bin list
            rename_contigs_job = SlurmJob(
                name="rename_binlist-" + str(assembly),
                account=account,
                cores=1,
                mem="8gb",
                dependency=",".join(bam_rename_job_list),
                cmd="rename_contigs_in_binlists.py {} {} {}".format(contig_map, bin_list, renamed_binlist)
            )
            rename_contigs_jobid = rename_contigs_job.submit()
            bin_list = renamed_binlist
            print(rename_contigs_jobid)

        # make contigs db
        gen_contigs_db_job = SlurmJob(
            name="gen_contigs_db-" + str(assembly),
            account=account,
            cores=20,
            mem="32gb",
            dependency=rename_contigs_jobid,
            cmd="anvi-gen-contigs-database -f {0} -o {1}/{2} -n '{2}' -T 20".format(renamed_assembly, contig_dbs, contigs_db)
        )
        gen_contigs_db_jobid = gen_contigs_db_job.submit()
        print(gen_contigs_db_jobid)
        
        #run hmms for completion estimates
        run_hmms_job = SlurmJob(
            name="run_hmms-" + str(assembly),
            account=account,
            cores=20,
            mem="32gb",
            dependency=gen_contigs_db_jobid,
            cmd="anvi-run-hmms -c {0}/{1} -T 20".format(contig_dbs, contigs_db)
        )
        run_hmms_jobid = run_hmms_job.submit()
        print(run_hmms_jobid)

        # make profiles
        for bam_path in bam_paths:
            bam = bam_path.split("/")[-1].split(".")[0]
            profile_bam_job = SlurmJob(
                name="profile_bam-" + str(assembly + "_" + bam),
                account=account,
                cores=36,
                mem="180gb",
                dependency=",".join([gen_contigs_db_jobid] + bam_rename_job_list),
                cmd="anvi-profile -i {} -c {}/{} -o {}/{}/{} --skip-hierarchical-clustering -T 36".format(bam_path, contig_dbs, contigs_db, profile_dbs, assembly, bam)
            )
            profile_bam_jobid = profile_bam_job.submit()
            print(profile_bam_jobid)
            profile_bams_job_list.append(profile_bam_jobid)
        
        # merge profiles
        merge_profiles_job = SlurmJob(
            name="merge_profiles-" + str(assembly),
            account=account,
            cores=1,
            mem="32gb",
            dependency=",".join([gen_contigs_db_jobid] + profile_bams_job_list),
            cmd="anvi-merge {0}/{1}/*/PROFILE.db -o {4}/{1} -c {2}/{3}".format(profile_dbs,assembly,contig_dbs,contigs_db,merged_profile_dbs)
        )
        merge_profiles_jobid = merge_profiles_job.submit()
        print(merge_profiles_jobid)

        # import bin list
        import_bins_job = SlurmJob(
            name="import_bins-" + str(assembly),
            account=account,
            cores=1,
            mem="8gb",
            dependency=merge_profiles_jobid,
            cmd="anvi-import-collection -c {}/{} -p ${}/{}/PROFILE.db -C concoct --contigs-mode {}".format(contig_dbs,contigs_db,merged_profile_dbs,assembly,bin_list)
        )
        import_bins_jobid = import_bins_job.submit()
        print(import_bins_jobid)