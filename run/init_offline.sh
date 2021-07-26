#!/bin/sh
set -x

##
## Default
##
LOCAL=${LOCAL:-3000.d0}
UPDATE_T=${UPDATE_T:-".false."}
UPDATE_Q=${UPDATE_Q:-".false."}
UPDATE_PS=${UPDATE_PS:-".false."}
UPDATE_UV=${UPDATE_UV:-".false."}
UPDATE_SST=${UPDATE_SST:-".false."}
UPDATE_LST=${UPDATE_LST:-".false."}

##
## clean
##
rm $(cd $DIR ; rm *.o *.exe)
rm $WDIR/gues/*/*.grd
rm $WDIR/anal/*/*.grd

##
## Simulation length & Time-stepping
##
SFX=${RES:$((${#RES}-2)):$((${#RES}-1))}
TRES=$(basename $RES $SFX)
if [ ${SFX} == "mo" ] ; then 
    NMONTS=${TRES}
    NSTOUT=-1
    NDAYS="ndays"
    NDAY=$(($TRES*30))
    CPP_OPT="OUT_FREQ_MON"
elif [ ${SFX} == "dy" ] ; then
    NMONTS=1
    NSTOUT=$(($TRES*48)) # x NSTEPS
    NDAYS=${TRES}
    NDAY=${TRES}
    CPP_OPT="OUT_FREQ_DAY"
fi

##
## Pre-process (cls_instep.h)
##
sed -e "s/<NMONTS>/${NMONTS}/" cls_instep.H |\
sed -e "s/<NSTOUT>/${NSTOUT}/" \
> cls_instep.h

##
## Pre-process (at_gcm.f)
##
sed -e "s/<NDAYS>/${NDAYS}/" at_gcm.F \
> at_gcm.f

##
## Pre-process (ppo_iogrid.f)
##
cpp -D${CPP_OPT} ppo_iogrid.F |\
awk 'NR>6{print $0}' \
> ppo_iogrid.f

##
## Pre-process (usr_iau.f)
##
sed -e "s/<UPDATE_T>/${UPDATE_T}/" usr_iau.F |\
sed -e "s/<UPDATE_Q>/${UPDATE_Q}/" |\
sed -e "s/<UPDATE_PS>/${UPDATE_PS}/" |\
sed -e "s/<UPDATE_UV>/${UPDATE_UV}/" |\
sed -e "s/<UPDATE_SST>/${UPDATE_SST}/" |\
sed -e "s/<UPDATE_LST>/${UPDATE_LST}/" |\
sed -e "s/<NDAY>/${NDAY}/" \
> usr_iau.f

##
## make speedy
##
make imp.exe

##
## compile
## 
cp $SRCDIR/letkf.f90 .
cp $SRCDIR/conv_hist_mon2ann.f . 
sed -e "s/<EDIM>/${EDIM}/" letkf.f90 |\
sed -e "s/<LOCAL>/${LOCAL}/" |\
sed -e "s/<INFLTYPE>/${INFLTYPE}/" |\
sed -e "s/<INFL>/${INFL}/" \
> tmp.f90
ifort -O3 -assume byterecl -convert big_endian -mkl tmp.f90 -o letkf.exe
ifort -O3 -assume byterecl -convert big_endian conv_hist_mon2ann.f -o mon2ann.exe

##
## Initial and External Conditions
## 
for stat in anal gues ; do
    mkdir -p $WDIR/$stat/mean
    mkdir -p $WDIR/$stat/sprd
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        mkdir -p $WDIR/$stat/$MEM
    done
done
cp letkf.exe mon2ann.exe $WDIR/

# exe
for stat in gues ; do
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        cp $DIR/imp.exe $WDIR/$stat/$MEM/
        cp $DIR/mon2ann.exe $WDIR/$stat/$MEM/
    done
done

# BC
for stat in gues ; do
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        cp $DATADIR/fort.* $WDIR/$stat/$MEM/
        echo 1    > $WDIR/$stat/$MEM/fort.2
        echo 000 >> $WDIR/$stat/$MEM/fort.2
    done
done

# IC
for stat in gues ; do
    for MEM in $(seq -f %03g 1 ${NE}) ; do
        cp $INITDIR/$MEM/1900010100.grd $WDIR/$stat/$MEM/fort.3
    done
done

