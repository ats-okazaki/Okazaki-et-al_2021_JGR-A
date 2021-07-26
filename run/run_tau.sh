#!/bin/bash
#PBS -j oe
#PBS -l nodes=1:ppn=30
#PBS -N tau_all
cd $PBS_O_WORKDIR

#set -ex

#-------------------------------------------------------------------------------------------------

#export INTEL_LICENSE_FILE=28518@172.19.6.1
export DIR=$PWD
export time_start=19000101
export NE=30    # ensemble size
export INFLTYPE="RTPS"  # choose from MULT, RTPP, RTPS

for LOCAL in 3000 2000 4000 1000 5000 6000 ; do
for INFL in 0.0 0.6 0.9 ; do
for RES in 1dy 3dy 10dy 1mo 3mo 6mo 12mo ; do

export RES=${RES}

#-------------------------------------------------------------------------------------------------

##
## Simulation length & Time-stepping
##
SFX=${RES:$((${#RES}-2)):$((${#RES}-1))}
TRES=$(basename $RES $SFX)
NSTEP=$(($TRES*300))
#NSTEP=$(($TRES*100))
if [ ${SFX} == "mo" ] ; then 
    tunit="month"
elif [ ${SFX} == "dy" ] ; then
    tunit="day"
fi
time_end=$(date --date "${time_start} ${NSTEP} ${tunit}" +%Y%m%d)

##
## Inputs 
##
export EXP="${RES}/L${LOCAL}.${INFLTYPE}${INFL}"
export WDIR=$DIR/$EXP
export SRCDIR="~/speedy_ver41.5/work/source"
export OBSDIR="~/speedy_ver41.5/work/nature_icsea-2/obs/snr2/${TRES}${SFX}"
export INITDIR="~/speedy_ver41.5/work/ens_spinup/gues"
export DATADIR="~/speedy_ver41.5/work/data"
export LOCAL=${LOCAL}
export INFL=${INFL}
export EDIM=${NE}

##
## check status
##
if [ $(ls ${EXP}/anal/mean | wc -l) -eq 101 -a $(ls ${EXP}/gues/mean | wc -l) -eq 101 ] ; then
    echo ${EXP} ' --skipped'
    continue
fi

##
## Initialization
##
cd $DIR
sh init_tau.sh
ln -s ~/speedy_ver41.5/work/obsnet/obslist_all.t30.bin $WDIR/obsnet

#-------------------------------------------------------------------------------------------------
##
## Run
##
time1=${time_start}
while (($(date --date $time1 +%Y%m%d) <= $(date --date ${time_end} +%Y%m%d))) ; do

    year=$(date --date $time1 +%Y)
    time_next=$(date --date "$time1 ${TRES} ${tunit}" +%Y%m%d)

    #--- distribute ---#
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        echo $time1 > $WDIR/gues/$MEM/fort.120
        echo "N"    > $WDIR/gues/$MEM/fort.121
    done

    #--- run speedy ---#
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        cd $WDIR/gues/$MEM
        mpirun -n 1 ./imp.exe  &
        #mpirun -np 1 dplace -s1 ./imp.exe  &
    done

    wait

    #--- obsope. ---#
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        cd $WDIR/gues/$MEM
        if [ ${SFX} = "mo" -a ${TRES} -gt 1 ]; then
            ./mon2ann.exe attm000_${year}.grd attm${TRES}${SFX}_${time1}.grd ${TRES}
        elif [ -f attm000_${year}.grd -a -s attm000_$((year+1)).grd ] ; then
            mv            attm000_$((year+1)).grd attm${TRES}${SFX}_${time1}.grd
        else
            mv            attm000_${year}.grd attm${TRES}${SFX}_${time1}.grd
        fi
    done

    #--- letkf ---#
    cd $WDIR
    #export omp_num_threads=${NE}
    time ./letkf.exe $time1 obsnet ${OBSDIR}/obs${TRES}${SFX}_${time1}.grd ${OBSDIR}/std${TRES}${SFX}.grd attm${TRES}${SFX}_${time1}.grd attm${TRES}${SFX}_${time1}.grd

    #--- time-mean update ---#
    cd $WDIR
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        ./tau.exe $WDIR/anal/${MEM}/attm${RES}_${time1}.grd $WDIR/gues/${MEM}/attm${RES}_${time1}.grd $WDIR/gues/${MEM}/${time_next}00.grd $WDIR/anal/${MEM}/${time_next}00.grd
    done

    #--- restart ---#
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        cd $WDIR/gues/$MEM
        cp $WDIR/anal/$MEM/${time_next}00.grd fort.3
        rm attm000_*.grd atva000_*.grd atdf000_*.grd
    done

    time1=${time_next}
    
done # time


##
## clean directory
##
for MEM in $(seq -f %03g 1 ${NE}) ; do
    cd $WDIR/gues/$MEM
    rm fort.*
    cd $WDIR/anal/$MEM
    rm fort.*
done

#-------------------------------------------------------------------------------------------------

done # res
done # infl
done # local


#--- report RMSE & SPRD ---#
grads -blc calc_rmse-sprd_brief.gs
