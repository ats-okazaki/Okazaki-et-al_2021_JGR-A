#!/bin/sh
#PBS -j oe
#PBS -l nodes=1:ppn=30
#PBS -N errg-E
cd $PBS_O_WORKDIR

#set -ex

# setting
ens=30
time_start_0=19000101

TOP=$PWD
DATADIR="/data/user/speedy_ver41.5/work/data"


#########################################################################
###
### compile
###


#########################################################################
###
### run
###
for exp in long ; do
for res in 3mo ; do
#for res in 1dy 3dy 10dy 1mo 3mo 6mo 12mo ; do

if [ ${exp} = 'ctrl' ] ; then
    INITDIR="/data/user/speedy_ver41.5/work/tau_all_snr2/"
elif [ ${exp} = 'long' ] ; then
    INITDIR="/data/user/speedy_ver41.5/work/tau_all_long_snr2"
fi

sfx=${res:$((${#res}-2)):$((${#res}-1))}
dt=$(basename $res $sfx)
if [ ${sfx} == 'mo' ] ; then
    tunit='month'
elif [ ${sfx} == 'dy' ] ; then
    tunit='day'
fi

# params
if [ ${exp} == 'ctrl' ] ; then
    if [ ${res} == '1dy' ] ; then PARAMS='L2000.RTPS0.9' ; 
    elif [ ${res} == '3dy' ]  ; then PARAMS='L2000.RTPS0.0' ;
    elif [ ${res} == '10dy' ] ; then PARAMS='L2000.RTPS0.9' ;
    elif [ ${res} == '1mo' ]  ; then PARAMS='L3000.RTPS0.6' ;
    elif [ ${res} == '3mo' ]  ; then PARAMS='L3000.RTPS0.0' ;
    elif [ ${res} == '6mo' ]  ; then PARAMS='L3000.RTPS0.9' ;
    elif [ ${res} == '12mo' ] ; then PARAMS='L3000.RTPS0.6' ; fi
elif [ ${exp} == 'long' ] ; then
    if [ ${res} == '1dy' ] ; then PARAMS='L2000.RTPS0.6' ; 
    elif [ ${res} == '3dy' ]  ; then PARAMS='L2000.RTPS0.0' ;
    elif [ ${res} == '10dy' ] ; then PARAMS='L3000.RTPS0.6' ;
    elif [ ${res} == '1mo' ]  ; then PARAMS='L3000.RTPS0.0' ;
    elif [ ${res} == '3mo' ]  ; then PARAMS='L4000.RTPS0.6' ;
    elif [ ${res} == '6mo' ]  ; then PARAMS='L3000.RTPS0.0' ;
    elif [ ${res} == '12mo' ] ; then PARAMS='L3000.RTPS0.9' ; fi
fi

#time_start=$(date --date "${time_start_0} $((${dt}*5)) ${tunit}" +%Y%m%d)
#time_end=$(date --date "${time_start_0} $((${dt}*65)) ${tunit}" +%Y%m%d)
time_start=19081001
time_end=19081001
echo $time_start $time_end
while (($(date --date ${time_start} +%Y%m%d) <= $(date --date ${time_end} +%Y%m%d))) ; do

    cd $TOP
    DIR=${TOP}/${exp}/${res}/${time_start}

    #--- Init ---#
    if [ ${sfx} = "mo" ] ; then
        sh init.sh $((dt*5))${sfx}
    elif [ ${sfx} = "dy" ] ; then
        sh init.sh $((dt*5))${sfx} "true"
    fi
    for stat in gues ; do
        for MEM in $(seq -f %03g 1 $ens) ; do
            mkdir -p $DIR/$stat/$MEM || true
            cp imp.exe $DIR/$stat/$MEM/
            # ic
            cp $INITDIR/$res/$PARAMS/anal/$MEM/${time_start}00.grd $DIR/$stat/$MEM/fort.3
            # bc
            cp $DATADIR/* $DIR/$stat/$MEM/
            # 
            echo 1    > $DIR/$stat/$MEM/fort.2
            echo 000 >> $DIR/$stat/$MEM/fort.2
            echo $time_start > $DIR/$stat/$MEM/fort.120
            echo "N"  > $DIR/$stat/$MEM/fort.121
        done
    done
    
    #--- Run SPEEDY ---#
    for MEM in $(seq -f %03g 1 $ens) ; do
        cd $DIR/gues/$MEM
        mpirun -np 1 dplace -s1 ./imp.exe &
    done

    wait # DO NOT REMOVE; necessary for mpirun
    
    #--- clean directory ---#
    for stat in gues ; do
        for MEM in $(seq -f %03g 1 $ens) ; do
            rm $DIR/$stat/$MEM/imp.exe
            rm $DIR/$stat/$MEM/mon2ann.exe
            rm $DIR/$stat/$MEM/fort.*
            rm $DIR/$stat/$MEM/atdf*
            rm $DIR/$stat/$MEM/atva*
            #rm $DIR/$stat/$MEM/??????????.grd
        done
    done

    time_start=$(date --date "${time_start} ${dt} ${tunit}" +%Y%m%d)

done

done # res
done # exp


