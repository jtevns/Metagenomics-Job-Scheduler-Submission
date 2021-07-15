from setuptools import setup, find_packages

setup(
    author='Jacob Evans',
    author_email='jtevans@umich.edu',
    description='cli for submitting batch slurm jobs',

    name='mgjss',
    version='0.1',
    packages=find_packages(),
    scripts=['bin/mgjss','scripts/rename_contigs_in_binlists.py','scripts/rename_contigs_in_bam_header.py'],

    install_requires=[],
)