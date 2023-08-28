#!/bin/bash
#SBATCH --job-name provadduu
#SBATCH -N1 --ntasks-per-node=1 --cpus-per-task=48 --mem=80GB
#SBATCH --time=24:00:00
#SBATCH --account=IscrC_CSoVCS
#SBATCH --partition=g100_usr_prod
# go to submission dir
cd $SLURM_SUBMIT_DIR
module load profile/chem-phys
module load autoload g16

. $g16root/g16/bsd/g16.profile       # for bash script

#Functions
trimfchk () {
    local chkname=$1
    local fchkname=${chkname%.chk}.fchk
    formchk $chkname 2>&1 > /dev/null
    awk '
        NR==1{flag1=1} /Atom fragment info /{flag1=0} flag1
        /Energy/
        /Cartesian Gradient/{flag2=1}/Nonadiabatic coupling/{flag2=0} flag2
        /Dipole Moment/{flag3=1} /NMR/{flag3=0} flag3
    ' $fchkname > ${fchkname%.fchk}_trim.fchk
    rm $fchkname
    rm $chkname
}

run_gaussian () {
        local jxyz=$1
        g16 ${jxyz}
        trimfchk ${jxyz%gjf}chk
}



export GAUSS_SCRDIR=$CINECA_SCRATCH/g16/test1  # def. tmp folder in $CINECA_SCRATCH
mkdir  -p $GAUSS_SCRDIR                      # the dir must exist

#put you input data in file test1.com, for example taking from g16 directory>
cp *.gjf $GAUSS_SCRDIR

pjob=4
cd $GAUSS_SCRDIR

name=()
#for i in ${prefix}*.gjf
for i in *.gjf
do
        name+=($i)
done
njob=${#name[@]}
let cycles=($njob+$pjob-1)/$pjob


for (( i=0; i<cycles; i++ ));
do
        for (( j=0; j<pjob; j++));
        do
                jobid=$(( i*pjob+j ))
                if [ $jobid -lt $njob ]; then
                        run_gaussian ${name[jobid]} &
                fi
        done
        wait
done
wait

cp *.fchk $SLURM_SUBMIT_DIR/
cd $SLURM_SUBMIT_DIR
rm -r $GAUSS_SCRDIR

