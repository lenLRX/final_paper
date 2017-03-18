import os
import shutil

for i in range(30,0,-1):
    dir_name = "./R1_" + str(i)
    os.mkdir(dir_name)
    shutil.copyfile("./a.out", dir_name + "/a.out")
    os.chdir(dir_name)
    f = open("job.pbs","w")
    f.write("#!/bin/sh\n#PBS -q blades\n")
    f.write("#PBS -N lentask%d\n"%i)
    f.write("#PBS -l nodes=1:ppn=1\n")
    f.write("cd $PBS_O_WORKDIR\n")
    f.write("date\n")
	f.write("source /public/software/profile.d/mathlib_mkl-2017.u1.sh");
    f.write("chmod +x ./a.out\n")
    f.write("./a.out %d\n"%i)
    f.write("date\n")
    f.close()
    os.system("qsub job.pbs")
    os.chdir("..")