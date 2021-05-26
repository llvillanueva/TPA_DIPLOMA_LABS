#!/bin/bash
set -e # kill the script on errors

# change the main variable and units
setname(){
    local file="$1"
    local oname="$2"
    local nname="$3"
    local lname="$4"
    local units="$5"

    ncrename -v ${oname},${nname} -h $file
    ncatted -O -a units,${nname},o,c,"${units}" $file
    ncatted -O -a long_name,${nname},o,c,"${lname}" $file
}


# functions for bias and RMSE
rmse(){
# this calculates the RMSE
local file1=$1
local file2=$2
outfile=RMSE_${file1%???}_${file2%???}.nc
cdo -sqrt -fldmean -sqr -sub $file1 $file2 $outfile
}

bias(){
# bias calc
local file1=$1
local file2=$2
outfile=BIAS_${file1%???}_${file2%???}.nc
cdo -fldmean -sub $file1 $file2 $outfile
}

# programme to parametrize Sundqvist scheme from ECMWF input
# Tompkins@ictp.it

rhcrit_list=`seq 0.65 0.02 0.95`
cwc_thresh=1.e-7
cdo="cdo -s"

# extract the fields
$cdo selvar,r ecmwf_data.nc rh.nc
$cdo selvar,cc ecmwf_data.nc tiedtkecc.nc
$cdo selvar,ciwc ecmwf_data.nc ciwc.nc
$cdo selvar,clwc ecmwf_data.nc clwc.nc

# calculate total cloud water
$cdo add ciwc.nc clwc.nc cwc.nc
echo setname cwc.nc clwc cwc "specific cloud water content" "kg kg**-1"

# scale RH to 0-1
$cdo divc,100 rh.nc rhf.nc

for rhcrit in $rhcrit_list ; do 

    echo $rhcrit 

    sfile1=sundqvistcc_mask0_rhc${rhcrit}.nc
    sfile2=sundqvistcc_mask1_rhc${rhcrit}.nc

    # 1-RH/1-RHcrit
    rhcritv=`echo 1 - $rhcrit | bc` 
    $cdo divc,${rhcritv} -addc,1 -mulc,-1 rhf.nc t1.nc

    # max of ratio and zero
    $cdo mul -gec,0 t1.nc t1.nc t2.nc

    # min of ratio and 1.0
    $cdo addc,1 -mul -lec,0 -subc,1 t2.nc -subc,1 t2.nc t3.nc

    # 1 - SQRT of this - SUNDQVIST
    $cdo addc,1 -mulc,-1 -sqrt t3.nc $sfile1
   
    # add mask for total cloud water - MASK SUNDQVIST
    $cdo mul -gec,$cwc_thresh cwc.nc $sfile1 $sfile2

    # change the meta data
    echo setname $sfile1 r cc "cloud cover" "fraction"
    echo setname $sfile2 r cc "cloud cover" "fraction"

    # diagnostics
    bias $sfile1 tiedtkecc.nc
    bias $sfile2 tiedtkecc.nc
    rmse $sfile1 tiedtkecc.nc
    rmse $sfile2 tiedtkecc.nc

    # clean up
    rm -f t?.nc 
done # loop over rhcrit
