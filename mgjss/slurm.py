# functions and classes for creating and submitting slurm 
# jobs

import subprocess

class SlurmJob:
    def __init__(self,name,account,cmd,email="$USER@umich.edu",time="72:00:00",cores="1",mem="8gb",dependency=""):
        self.name = name
        self.email = email
        self.account = account
        self.time = time
        self.cores = cores
        self.mem = mem
        self.dependency = dependency
        self.cmd = cmd

    def submit(self):
        dep = ""
        if self.dependency != '':
            dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(self.dependency)
 
        sbatch_command = "sbatch --job-name {} --output={}_slurm.log --mail-user={} --mail-type=FAIL --partition=standard --account={} --time {} --mem={} --cpus-per-task={} --wrap='{}' {}".format(self.name, self.name, self.email, self.account, self.time, self.mem, self.cores, self.cmd, dep)
        sbatch_response = subprocess.getoutput(sbatch_command)
        job_id = sbatch_response.split(' ')[-1].strip()
        return job_id