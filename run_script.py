from glob import glob
import subprocess as subproc
from shlex import shlex
import os
import time

for pwscf_infile in glob("*.pwscf.in"):
    if pwscf_infile.strip("in") + "out" not in glob("*.pwscf.out"):
        jobname = pwscf_infile.strip(".pwscf.in")
        with open("run_script.sh", "w") as rs:
            rs.write("#!/bin/bash\n")
            rs.write("#PBS -l walltime=04:00:00\n")
            rs.write("#PBS -l nodes=1:ppn=16\n")
            rs.write("#PBS -N QE\n")
            rs.write("#PBS -q standby\n")
            rs.write("#PBS -N %s\n" %(jobname))
#            rs.write("#PBS -m ea\n")
            rs.write("set echo\n")
            rs.write("cd $PBS_O_WORKDIR\n")
            rs.write("module load intel/16.0.1.150 impi/5.1.2.150 espresso/6.0\n")
            rs.write("mpirun -np 16 pw.x < %s > %s" %(pwscf_infile, pwscf_infile.strip(".in")+".out"))
    
        command = "qsub " + "run_script.sh"
        tmp_proc = subproc.Popen(["qsub", "run_script.sh"])
        tmp_proc.wait()
