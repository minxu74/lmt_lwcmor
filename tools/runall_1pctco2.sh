#!/usr/bin/env bash



export HDF5_USE_FILE_LOCKING=FALSE

#topdir=/global/homes/m/minxu/scratch/C4MIP_CMOR/20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/
#topdir=/global/homes/m/minxu/scratch/C4MIP_CMOR/20190320_gswpv2_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/
#topdir="/global/homes/m/minxu/scratch/C4MIP_CMOR/20190224_princeton_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/"
#topdir="/global/homes/m/minxu/scratch/C4MIP_CMOR/20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC"
#expnam=''

prefix=/global/cfs/cdirs/m2467/prj_minxu/1pctCO2/

#-topdirs=(\
#-20180910_cruncepv8_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
#-20180910_princeton_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
#-20190224_princeton_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC \
#-20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
#-20190320_gswpv2_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC \
#-20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC) 

#topdirs=(\
#20180910_cruncepv8_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
#20180910_princeton_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
#20190224_princeton_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC \
#20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC) 

#topdirs=(20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC) 
#topdirs=(\
#20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
#20190320_gswpv2_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC)

topdirs=(\
20200509.CO21PCTBGC_RUBISCO_CNCTC20TR_OIBGC.ne30_oECv3.compy  \
20200802.CO21PCTFUL_RUBISCO_CNCTC20TR_OIBGC.ne30_oECv3.compy)

wkdir=$(pwd)

for topdir in "${topdirs[@]}"; do
    echo $topdir
    
    trange=000101-015012
    expnam=''
    if [[ "$topdir" == *"PCTBGC"* ]]; then
       expnam='BGC'
    fi

    if [[ "$topdir" == *"PCTFUL"* ]]; then
       expnam=''
    fi

    echo $expnam $trange

    if [[ ! -d "$prefix/$topdir/org/lnd/Amon" ]]; then
       mkdir -p $prefix/$topdir/org/lnd/Amon
    fi

    if [[ ! -d "$prefix/$topdir/org/lnd/Lmon" ]]; then
       mkdir -p $prefix/$topdir/org/lnd/Lmon
    fi

    if [[ ! -d "$prefix/$topdir/org/lnd/Emon" ]]; then
       mkdir -p $prefix/$topdir/org/lnd/Emon
    fi

    if [[ ! -d "$prefix/$topdir/org/lnd/LImon" ]]; then
       mkdir -p $prefix/$topdir/org/lnd/LImon
    fi

    cd $prefix/$topdir/org/lnd/Amon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/Lmon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/Emon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/LImon && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir


    #-python c4mip_cmor.py  --datadir $prefix/$topdir/rgr/ \
    #-                      --cmordir $prefix/$topdir/cmor/ \
    #-                      --mipname C4MIP \
    #-                      --modname E3SM \
    #-                      --expname 1pctCO2${expnam} \
    #-                      --insname RUBISCO\
    #-                      --tabname fx\
    #-                      --relname atmos\
    #-                      --tmrange $trange \
    #-                      --sourceid E3SM-1-1


    python c4mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname Lmon\
                          --relname land\
                          --tmrange $trange \
                          --sourceid E3SM-1-1


    python c4mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname Amon\
                          --relname atmos\
                          --tmrange $trange \
			  --sourceid E3SM-1-1


    python c4mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname Emon\
                          --relname land\
                          --tmrange $trange \
			  --sourceid E3SM-1-1
 
    #LImon
    python c4mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname LImon\
                          --relname land\
                          --tmrange $trange \
                          --sourceid E3SM-1-1


done
