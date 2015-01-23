#************************************************************
#
#                     Run Job on PBS Systems
#
# Methodology: Instead of submitting the program directly, 
# submit a script as the initial job. The script launches
# your program and many other operations and schedules that
# you wanted. This approach will provide you the most 
# flexibility for running your case.
#
# This script can be executed with arguments, these arguments
# will all be passed to the executable program in the form that
# these arguments will show after the program name.
#
#************************************************************

#! /bin/bash
set -e

#******************** Set Job Parameters ********************
#
# Delete characters between quotes to omit an unwanted item.
#

#* the command line for running your program
Executable_Pragram="mpirun starccmplus"

#* job type: [mpi] or [serial] or [threaded]
Job_Type="mpi"

#* extra memory request
Extra_Memory_Request="6g"

#* number of processors need to use
Number_Of_Processors="2"

#* how long need to run: [h] for hours, [d] for days
Time_To_Run="1h"

#* output file for standard output of the program
Output_File="output"

#* error file for standard error stream of the program
Error_File="errfile"

#* special pool to link
Link_Pool="NRAP_1213"

#******************* Extra Information *********************

# show your jobs: [sqjobs] options: [-a], [-u someuser]

# to kill, suspend or resume your jobs, use 
# sqkill/suspend/resume with the job ID as shown by sqjobs. 

# run program interactively by submitting a screen bash 
# command (screen -D -fn -m) as a job. set:
# Executable_Program="-o /dev/null screen -D -fn -m bash"

# show executing node: qstat -f -l <JOBID> | egrep exec_host

# ssh to node and attach the running screen session:
# ssh -t <NODE> screen -r

# run commands on all allocated nodes: pbsdsh -o <COMMAND>

# start MPI programs on nodes: mpirun <COMMAND>

# check system default MPI library: sqsub -vd ...

# /home and /work directories on SHARCNET are remote to 
# clusters, thus, to obtain the best file system throughput
# should use the /scratch file system, but also notice that
# /scratch does not have unified access and expiries for
# each two months. Use /tmp if possible as it's a filesystem
# that is mounted from the node's local disk.

# set the permissions on the base directory to only be 
# accessible to yourself: chmod 700 ~/

# How to archive my data? 
# cp /scratch/$USER/$SIMULATION /archive/$USER/$SIMULATION

#############################################################
#    generally, following configures do not need to edit
#############################################################

#********** obtain the correct form of parameters **********
#
# set the job type, use serial as default. 
# using the --nompirun flag when submitting parallel jobs
# will not auto trigger the mpirun of the job, therefore,
# it enables to encapsulate MPI jobs in a shell script.
#

JobType=""
if [[ $Job_Type == "thread" ]]; then
    JobType="-q thread"
elif [[ $Job_Type == "mpi" ]]; then
    JobType="-q mpi --nompirun"
else
    JobType="-q serial"
fi

ExtraMemoryRequest=""
if [[ -n $Extra_Memory_Request ]]; then
    ExtraMemoryRequest="--mpp=$Extra_Memory_Request"
fi

NumberOfProcessors=""
if [[ -n `echo $Number_Of_Processors | grep "^[0-9]*$"` ]]; then
    NumberOfProcessors="-n $Number_Of_Processors"
else
    echo "illegal number of processors, exit..."
    exit
fi

TimeToRun=""
if [[ -n $Time_To_Run ]]; then
    TimeToRun="-r $Time_To_Run"
else
    echo "illegal time to run, exit..."
    exit
fi

OutputFile=""
if [[ -n $Output_File ]]; then
    OutputFile="-o $Output_File"
fi

ErrorFile=""
if [[ -n $Error_File ]]; then
    ErrorFile="-e $Error_File"
fi

if [[ -n $Link_Pool ]]; then
    JobType="-q "$Link_Pool" -f mpi --nompirun"
fi

#***************** generate the scripts file ****************
ScriptDir="."
QsubFile="execution.sh"
echo "******************************************************"
if [[ -f $ScriptDir/$QsubFile ]]; then
    echo "use the existing scripts file: $QsubFile..."
else
    echo "generate and use the scripts file: $QsubFile..."
cat > $ScriptDir/$QsubFile <<EOF
#************************************************************
#
#                Current Job Schedule
#
# Techniquely, you can do whatever you want with the allocated
# resources as part of a job - multiple MPI subjobs, serial 
# sections, etc. However,  the non-MPI portions of this script
# will run serially. This wastes cycles on all but one of the
# processors - a serious concern for long serial sections and
# /or jobs with many cpus.
#
#************************************************************

#! /bin/bash

echo "start at  \`date +'%F %k:%M:%S'\`"

# this part indicates any preprocessing.
echo "preprocessing..."

# now run programs
echo "running..."
$Executable_Pragram $@

# this part indicases any postprocessing.
echo "postprocessing..."

echo "finish at  \`date +'%F %k:%M:%S'\`"
EOF
    chmod +x "$ScriptDir/$QsubFile"
fi
echo "sqsub $JobType $NumberOfProcessors $ExtraMemoryRequest $TimeToRun $OutputFile $ErrorFile ./$QsubFile"
echo "sqsub $JobType $NumberOfProcessors $ExtraMemoryRequest $TimeToRun $OutputFile $ErrorFile ./$QsubFile" > "$Output_File"
sqsub $JobType $NumberOfProcessors $ExtraMemoryRequest $TimeToRun $OutputFile $ErrorFile ./$QsubFile
echo "Job submitted!"
echo "******************************************************"

