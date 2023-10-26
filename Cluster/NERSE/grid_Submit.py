#!/usr/bin/env python
#author: Wenbin Zhao
#email: wenbinzhao237@gmail.com

from subprocess import call
import sys
import random
import time
def submit(fold_id_start=0, nfold = 1, nevent = 10, random_number = 0):
    jobs = '''#!/bin/bash

hostname
date
 

# First link the file

fold_id_start2=$(({fold_id_start} + $1 * {nfold}))
fold_id_end=$(({fold_id_start} + $1 * {nfold} + {nfold}))
# Then run the framework
for (( ii=$fold_id_start2; ii<$fold_id_end; ii++ ))
do

cp -r event0 Playground/job-$ii
cd Playground/job-$ii

# Generate the pythia parton
#sleep 1s
cd pythia_parton
./mymain06 {nevent} $(({random_number} + $ii * 12345 + $1 * 11))
cd ../

# ZPC for parton cascade
cd  ZPC
mkdir ana
ln -s ../pythia_parton/parton_info.dat ./
./exec
rm -r ana/parton-collisionsHistory.dat
rm -r ana/zpc.res
cd ../
rm -rf pythia_parton/parton_info.dat

# fragmentation and urqmd
cd hadronization_urqmd
cd fragmentation
ln -s ../../ZPC/ana/zpc.dat ./
./main_string_fragmentation {nevent}
rm -r ../../ZPC/ana/*

cd ../urqmd_code
    # script to run urqmd
    cd osc2u
    ln -s ../../fragmentation/hadrons_frag1.dat ./
    ./osc2u.e < hadrons_frag1.dat > run.log
    rm -r ../../fragmentation/hadrons_frag1.dat
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh > run.log
    rm -fr OSCAR.input
    rm -rf run.log
    cd ..
cd ../
cd ../
# jet finding of final hadrons
cd fastjet_hadron
ln -s ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat ./
./fastjet_hadron {nevent}
rm -r ../hadronization_urqmd/urqmd_code/urqmd/particle_list.dat
cd ../

# Save the final results into folder
mkdir results
mv fastjet_hadron/final_state_hard_hadrons.bin ./results/$ii.bin
rm -r fastjet_hadron
rm -r hadronization_urqmd 
rm -r pythia_parton 
rm -r ZPC 
cd ../../
done
'''.format(fold_id_start=fold_id_start, nfold = nfold, nevent = nevent, random_number = random_number)
    job_name = "NSC3_%s.sh"%(fold_id_start)
    with open(job_name, 'w') as fout:
        fout.write(jobs)

    jobs='''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --account=m3148
#SBATCH --error=/pscratch/sd/w/wenbinz/V2inJet/log/%J.err
#SBATCH --output=/pscratch/sd/w/wenbinz/V2inJet/log/%J.out
#SBATCH --array=0-{nfold}
chmod -R 777 {job_name}
./{job_name} $SLURM_ARRAY_TASK_ID
    '''.format(job_name = job_name, nfold = nfold)
    job_name2 = "Submit_%s.sh"%(fold_id_start)
    with open(job_name2, 'w') as fout:
        fout.write(jobs)
    call(['sbatch', job_name2])
    #call(['mv', job_name, 'jobs/'])
    #call(['mv', job_name2, 'jobs/'])
# comparing the grid size dependence of the hypersf cube
# etaos_ymin = 0.08 for both ampt runs, not 0.2

if __name__=='__main__':
    import sys
    fold_id_start = int(int(sys.argv[1]))
    nfold = int(int(sys.argv[2]))
    nevent = int(int(sys.argv[3]))
    # Get the current system time as a seed
    seed = int(time.time())

    # Seed the random number generator
    random.seed(seed)

    # Generate a random number between 0 and 10,000,000
    random_number = random.randint(0, 10**6)
    #for n in range(0,nods):
    #mmid = fold_id_start * tot_ev + n 
    submit(fold_id_start * nfold * nfold, nfold, nevent, random_number + fold_id_start)
    #start_num += jobs_per_cpu 

