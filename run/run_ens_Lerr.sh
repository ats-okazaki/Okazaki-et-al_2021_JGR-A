#!/bin/sh
#PBS -j oe
#PBS -l nodes=1:ppn=30
#PBS -N errg-L
cd $PBS_O_WORKDIR

export INTEL_LICENSE_FILE=28518@172.19.6.1

set -ex

# setting
ens=30
#time_start=18700101
#time_end=18700201
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
for exp in ctrl ; do
#for res in 3mo; do
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


#time_start=18700201
time_start=18700101
while (($(date --date ${time_start} +%Y%m%d) <= $(date --date ${time_end} +%Y%m%d))) ; do

    time_start_init=$(date --date "$time_start ${dt} ${tunit} ago" +%Y%m%d)
    year=$(date --date $time_start_init +%Y)
    mon=$(date --date $time_start_init +%m)

    cd $TOP
    DIR=${TOP}/${exp}/${res}/${time_start}

    #--- Init (1st) ---#
    sh init.sh ${res}
    for stat in gues ; do
        for MEM in $(seq -f %03g 1 $ens) ; do
            mkdir -p $DIR/$stat/$MEM || true
            cp imp.exe $DIR/$stat/$MEM/
            cp mon2ann.exe $DIR/$stat/$MEM/
            # bc
            cp $DATADIR/* $DIR/$stat/$MEM/
            # 
            echo 1    > $DIR/$stat/$MEM/fort.2
            echo 000 >> $DIR/$stat/$MEM/fort.2
            echo $time_start_init  > $DIR/$stat/$MEM/fort.120
            echo "N"  > $DIR/$stat/$MEM/fort.121
            # ic
            cp $INITDIR/${MEM}/${year}${mon}0100.grd $DIR/$stat/${MEM}/fort.3
        done
    done
    
    #--- Run SPEEDY (1st) ---#
    for MEM in $(seq -f %03g 1 $ens) ; do
        cd $DIR/gues/$MEM
        mpirun -np 1 dplace -s1 ./imp.exe &
    done
    
    wait # DO NOT REMOVE; necessary for mpirun


    #--- time-mean ---#
    for MEM in $(seq -f %03g 1 ${ens}) ; do
        cd $DIR/gues/$MEM
        if [ ${sfx} = "mo" -a ${dt} -gt 1 ]; then
            if [ $(ls -l attm000_${year}.grd | awk '{print $5}') -lt $((1917760*${dt})) \
                -a -s attm000_$((year+1)).grd ] ; then 
                cat attm000_$((year+1)).grd >> attm000_${year}.grd
            fi
                ./mon2ann.exe attm000_${year}.grd attm${dt}${sfx}_${time_start}.grd ${dt}
        elif [ -f attm000_${year}.grd -a -s attm000_$((year+1)).grd ] ; then
            mv            attm000_$((year+1)).grd attm${dt}${sfx}_${time_start}.grd
        else
            mv            attm000_${year}.grd attm${dt}${sfx}_${time_start}.grd
            if [ ! -s attm${dt}${sfx}_${time_start}.grd ] ; then
                echo $DIR/$stat/$MEM/attm${dt}${sfx}_${time_start}.grd
                stop
            fi
        fi
    done


    #--- Init (2nd) ---#
    cd $TOP
    if [ ${sfx} = "mo" ] ; then
        sh init.sh $((dt*5))${sfx}
    elif [ ${sfx} = "dy" ] ; then
        sh init.sh $((dt*5))${sfx} "true"
    fi
    for stat in gues ; do
        for MEM in $(seq -f %03g 1 $ens) ; do
            cp imp.exe $DIR/$stat/$MEM/
            echo $time_start > $DIR/$stat/$MEM/fort.120
            # ic
            ./tau.exe $DIR/$stat/030/attm${dt}${sfx}_${time_start}.grd \
                      $DIR/$stat/$MEM/attm${dt}${sfx}_${time_start}.grd \
                      $DIR/$stat/$MEM/${time_start}00.grd \
                      $DIR/$stat/$MEM/a${time_start}00.grd
            cp $DIR/$stat/$MEM/a${time_start}00.grd $DIR/$stat/$MEM/fort.3
            (cd ${DIR}/${stat}/${MEM} ; rm attm000_*.grd atva* atdf* || true)
        done
    done


    #--- Run SPEEDY (2nd) ---#
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
