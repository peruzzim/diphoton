#!/bin/bash
#################################
# PSI Tier-3 example batch Job  #
#################################

##### CONFIGURATION ##############################################
# Output files to be copied back to the User Interface
# (the file path must be given relative to the working directory)
#OUTFILES="myout.txt myerr.txt"

# Output files to be copied to the SE after your job has run
# (as above the file path must be given relative to the working directory)
# the file will be named with the same base name on the SE
#SEOUTFILES="mybigfile"
#
# By default, the files will be copied to $USER_SRM_HOME/$username/$JOBDIR,
# but here you can define the subdirectory under your SE storage path
# to which files will be copied (uncomment option)
#SEUSERSUBDIR="mytestsubdir/somedir"
#
# User's CMS hypernews name (needed for user's SE storage home path
# USER_SRM_HOME below)
HN_NAME=YOUR_CMS_HN_NAME

# set DBG=1 for additional debug output in the job report files
# DBG=2 will also give detailed output on SRM operations
DBG=0

#### The following configurations you should not need to change
# The SE's user home area (SRMv2 URL)
USER_SRM_HOME="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/"

# Top working directory on worker node's local disk. The batch
# job working directory will be created below this
TOPWORKDIR=/scratch/`whoami`

# Basename of job sandbox (job workdir will be $TOPWORKDIR/$JOBDIR)
JOBDIR=sgejob-$JOB_ID
##################################################################


############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# # Job name (defines name seen in monitoring by qstat and the
# #     job script's stderr/stdout names)
# #$ -N example_job
# 
# # Run time soft and hard limits hh:mm:ss
# # soft=CPU time; hard=Wallclock time
# #$ -l s_rt=7:00:00,h_rt=07:05:00
# 
# ### Specify the queue on which to run
# #$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#  #$ -o /shome/username/mydir/
#  #$ -e /shome/username/mydir/

##################################################################



##### MONITORING/DEBUG INFORMATION ###############################
DATE_START=`date +%s`
echo "Job started at " `date`
cat <<EOF
################################################################
## QUEUEING SYSTEM SETTINGS:
HOME=$HOME
USER=$USER
JOB_ID=$JOB_ID
JOB_NAME=$JOB_NAME
HOSTNAME=$HOSTNAME
TASK_ID=$TASK_ID
QUEUE=$QUEUE

EOF

echo "################################################################"

if test 0"$DBG" -gt 0; then
   echo "######## Environment Variables ##########"
   env
   echo "################################################################"
fi


##### SET UP WORKDIR AND ENVIRONMENT ######################################
STARTDIR=`pwd`
WORKDIR=$TOPWORKDIR/$JOBDIR
RESULTDIR=$STARTDIR/$JOBDIR
if test -e "$WORKDIR"; then
   echo "ERROR: WORKDIR ($WORKDIR) already exists! Aborting..." >&2
   exit 1
fi
mkdir -p $WORKDIR
if test ! -d "$WORKDIR"; then
   echo "ERROR: Failed to create workdir ($WORKDIR)! Aborting..." >&2
   exit 1
fi

cd $WORKDIR
cat <<EOF
################################################################
## JOB SETTINGS:
STARTDIR=$STARTDIR
WORKDIR=$WORKDIR
RESULTDIR=$RESULTDIR
EOF

###########################################################################
## YOUR FUNCTIONALITY CODE GOES HERE


source /swshare/ROOT/thisroot.sh
cd /shome/peruzzi/shape_studies/theory_marco || exit 1
root -q -b -l ~/shape_studies/template_production.C+O\(\"$1\",\"effunf\",0,\"$2\"\)

echo "Cleaning up $WORKDIR"
rm -rf $WORKDIR

###########################################################################
DATE_END=`date +%s`
RUNTIME=$((DATE_END-DATE_START))
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0
