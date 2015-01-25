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
# they show after the program name.
#
#************************************************************

#! /bin/bash
set -e

#******************** Set Job Parameters ********************
#
# Delete characters between quotes to omit an unwanted item.
#

#* the command line for running your program
Executable_Pragram="mpirun artracfd"

#* queue name: [serial], [threaded], [mpi], [nonstandard queues]
Queue_Name="mpi"

#* specify certain flags to modify behavior. [mpi], [gpu],
# [threaded] (these flags will be ignored for standard queues)
Specify_Flags="mpi"

#* run time: [h] hours, [d] days, [runtime] unlimited
Time_To_Run="1h"

#* number of CPU processors (it will be ignored if [serial])
Number_Of_CPU_Processors="2"

#* number of GPU processors (works only if [gpu])
Number_Of_GPU_Processors=""

#* output file for standard output of the program
Output_File="output"

#* error file for standard error stream of the program
Error_File="errfile"

#* extra memory request, default is [2g]
Extra_Memory_Request=""

#* test mode: [false], [true]
Test_Mode="false"

#* wait for a list of jobs to complete. [JobID[,JobID...]]
Wait_For=""

#* record JobID of current job to file
JobID_File="JobID"

#* provides a name for the current job.
Job_Name=`pwd`

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
# set the queue name and related flags 
#
# using the --nompirun flag when submitting parallel jobs
# will not auto trigger the mpirun of the job, therefore,
# it enables to encapsulate MPI jobs in a shell script.
#

CPUProcessors=""
if [[ -n `echo $Number_Of_CPU_Processors | grep "^[0-9]*$"` ]]; then
    CPUProcessors="$Number_Of_CPU_Processors"
else
    if [[ $Specify_Flags != "gpu" ]]; then
        echo "illegal number of CPU processors, exit..."
        exit
    fi
fi

GPUProcessors=""
if [[ $Specify_Flags == "gpu" ]]; then
    if [[ -n `echo $Number_Of_GPU_Processors | grep "^[0-9]*$"` ]]; then
        GPUProcessors="$Number_Of_GPU_Processors"
    else
        echo "illegal number of GPU processors, exit..."
        exit
    fi
fi

QueueName=""
if [[ $Queue_Name == "serial" ]]; then
    QueueName="-q serial"
elif [[ $Queue_Name == "threaded" ]]; then
    QueueName="-q threaded -n $CPUProcessors"
elif [[ $Queue_Name == "mpi" ]]; then
    QueueName="-q mpi -n $CPUProcessors --nompirun"
else
    if [[ $Specify_Flags != "gpu" ]]; then
        QueueName="-q $Queue_Name -f $Specify_Flags -n $CPUProcessors --nompirun"
    else
        QueueName="-q $Queue_Name -f $Specify_Flags --gpp=$GPUProcessors"
    fi
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

ExtraMemoryRequest=""
if [[ -n $Extra_Memory_Request ]]; then
    ExtraMemoryRequest="--mpp=$Extra_Memory_Request"
fi

WaitFor=""
if [[ -n $Wait_For ]]; then
    WaitFor="-w $Wait_For"
fi

TestMode=""
if [[ $Test_Mode == "true" ]]; then
    TestMode="--test"
fi

JobIDFile=""
if [[ -n $JobID_File ]]; then
    JobIDFile="--idfile=$JobID_File"
fi

JobName=""
if [[ -n $Job_Name ]]; then
    JobName="-j $Job_Name"
fi
#***************** generate the scripts file ****************
ScriptDir="."
QsubFile="execution.sh"
echo "******************************************************"
if [[ -f $ScriptDir/$QsubFile ]]; then
    echo "Submit the existing script file: $QsubFile..."
else
    echo "Generate and submit the script file: $QsubFile..."
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
echo "sqsub $TestMode $QueueName $ExtraMemoryRequest $TimeToRun $OutputFile $ErrorFile $WaitFor $JobIDFile $JobName ./$QsubFile"
echo "sqsub $TestMode $QueueName $ExtraMemoryRequest $TimeToRun $OutputFile $ErrorFile $WaitFor $JobIDFile $JobName ./$QsubFile" > "$Output_File"
sqsub $TestMode $QueueName $ExtraMemoryRequest $TimeToRun $OutputFile $ErrorFile $WaitFor $JobIDFile $JobName ./$QsubFile
echo "Job submitted!"
echo "******************************************************"

