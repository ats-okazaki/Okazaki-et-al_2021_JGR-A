#!/bin/sh
#PBS -j oe
#PBS -l nodes=1:ppn=30
#PBS -N errg-noint
cd $PBS_O_WORKDIR

export INTEL_LICENSE_FILE=28518@172.19.6.1

set -ex

# setting
ens=30
time_start=18700101
time_end=18750101

TOP=$PWD
DATADIR="/data/user/speedy_ver41.5/work/data"


#########################################################################
###
### compile
###

ifort -O3 -assume byterecl -convert big_endian usr_tau.f90 -o tau.exe
ifort -O3 -assume byterecl -convert big_endian conv_hist_mon2ann.f -o mon2ann.exe

#########################################################################
###
### run
###
for exp in long ; do
for res in 1dy 3dy 10dy 1mo 3mo 6mo 12mo ; do

if [ ${exp} = 'ctrl' ] ; then
    INITDIR="/data/user/speedy_ver41.5/work/ens_spinup/gues"
elif [ ${exp} = 'long' -o ${exp} == 'short' ] ; then
    INITDIR="/data/user/speedy_ver41.5/work/ens_spinup_${exp}/gues"
fi

sfx=${res:$((${#res}-2)):$((${#res}-1))}
dt=$(basename $res $sfx)
if [ ${sfx} == 'mo' ] ; then
    tunit='month'
elif [ ${sfx} == 'dy' ] ; then
    tunit='day'
fi


time_start=18700101
while (($(date --date ${time_start} +%Y%m%d) <= $(date --date ${time_end} +%Y%m%d))) ; do

    year=$(date --date $time_start +%Y)
    mon=$(date --date $time_start +%m)

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
            # bc
            cp $DATADIR/* $DIR/$stat/$MEM/
            # 
            echo 1    > $DIR/$stat/$MEM/fort.2
            echo 000 >> $DIR/$stat/$MEM/fort.2
            echo $time_start > $DIR/$stat/$MEM/fort.120
            echo "N"  > $DIR/$stat/$MEM/fort.121
            # ic
            cp $INITDIR/${MEM}/${year}${mon}0100.grd $DIR/$stat/${MEM}/fort.3
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
            #rm $DIR/$stat/$MEM/mon2ann.exe
            rm $DIR/$stat/$MEM/fort.*
            rm $DIR/$stat/$MEM/atdf*
            rm $DIR/$stat/$MEM/atva*
            #rm $DIR/$stat/$MEM/??????????.grd
        done
    done

    time_start=$(date --date "$time_start 1 month" +%Y%m%d)

done

done # res
done # exp


##--- clean directory ---#
#time_start=18600101
#while (($(date --date ${time_start} +%Y%m%d) <= $(date --date ${time_end} +%Y%m%d))) ; do
#    DIR=${TOP}/${time_start}
#    for stat in gues ; do
#        for MEM in $(seq -f %03g 1 $ens) ; do
#            rm $DIR/$stat/$MEM/imp.exe
#            rm $DIR/$stat/$MEM/fort.*
#            rm $DIR/$stat/$MEM/atdf*
#            rm $DIR/$stat/$MEM/atva*
#            rm $DIR/$stat/$MEM/??????????.grd
#        done
#    done
#    time_start=$(date --date "$time_start 1 month" +%Y%m%d)
#done



##--- Ensemble Mean & Spread---#
#while (($(date --date $time1 +%Y%m%d) <= $(date --date $tend +%Y%m%d))) ; do
#
#    year=$(date --date $time1 +%Y)
#    time_next=$(date --date "$time1 12 month" +%Y%m%d)
#
#    #--- time-mean  ---#
#    for MEM in $(seq -f %03g 1 $ens) ; do
#        cd $DIR/gues/$MEM
#        ./mon2ann.exe attm000_${year}.grd attm1yr_${time1}.grd 12
#    done
#
#    wait
#
#    #--- Ensemble mean and spread ---#
#    cd $DIR
#    ./stats.exe gues attm1yr_${time1}.grd
#
#    wait
#
#    time1=${time_next}
#
#done
#
