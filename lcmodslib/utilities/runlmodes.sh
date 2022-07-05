#!/bin/bash
# author m.fuse
# exe


### Functions ###
Help()
{
   # Display Help
   echo "Runs the jobs required for lmodes analysis"
   echo "and save the reuired quantities in trimmed fchk"
   echo "to reduce the disk usage."
   echo "Syntax: runlmodes.sh prefix ncpu"
   echo 
   echo "prefix          The prefix of Gaussian input files"
   echo "njob            The number of jobs to run in parallel"
   echo "options:"
   echo "h     Print this Help."
   echo
}

check_integer (){
        if [ "$1" -eq "$1" ] 2>/dev/null
        then
                echo 0
        else
                echo 1
        fi
}

trimfchk () {
    local chkname=$1
    local fchkname=${chkname%.chk}.fchk
    ${GAUSS_EXEDIR}/formchk $chkname 2>&1 > /dev/null
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
        $GAUSS_EXEDIR/g16 ${jxyz}
        trimfchk ${jxyz%gjf}chk
}

### MAIN ###
while getopts ":h" option; do
    case $option in
        h) Help
           exit;;
        \?) echo "Error: Invalid option"
           exit;;
    esac
done

argv=("$@")
nargv=${#argv[@]}
if [ $nargv -lt 2 ]
then
        echo "not enough argument"
        exit 1
fi
prefix=$argv[1]
pjob=$argv[2]
#Checks the integers
isint=$( check_integer $pjob )
if ![[ $isint -eq 0 ]]
then
    echo $pjob "Not a number."
    Help
    exit 1
fi

name=()
for i in ${prefix}*.gjf
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
                        run_gaussian ${name[jobid]} $PBS_O_WORKDIR &
                fi
        done
        wait
done
wait
