#!/usr/bin/env bash


export HDF5_USE_FILE_LOCKING=FALSE

prefix=/global/homes/m/minxu/scratch/LS3MIP_PUBLISH_TEST/

topdirs=(\
20180910_cruncepv8_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20180910_princeton_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC)

wkdir=$(pwd)

for topdir in "${topdirs[@]}"; do
    echo $topdir
    
    trange=185001-201412
    expnam=''
    if [[ "$topdir" == *"_cruncep"* ]]; then
       expnam='-cruNcep'
       trange=185001-201612
    fi

    if [[ "$topdir" == *"_princeton"* ]]; then
       expnam='-princeton'
       trange=185001-201212
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

    cd $prefix/$topdir/org/lnd/Amon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/Lmon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/Emon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/LImon && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir


    python trendy_cmor.py --datadir $prefix/$topdir/org/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname TRENDY \
                          --modname E3SM \
                          --expname land-hist${expnam} \
                          --insname E3SM-Project\
                          --tabname fx\
                          --relname atmos\
                          --tmrange $trange \
                          --sourceid E3SM-2-0


    python ls3mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname TRENDY \
                          --modname E3SM \
                          --expname land-hist${expnam} \
                          --insname E3SM-Project\
                          --tabname Lmon\
                          --relname land\
                          --tmrange $trange \
                          --sourceid E3SM-2-0


    python ls3mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname TRENDY \
                          --modname E3SM \
                          --expname land-hist${expnam} \
                          --insname E3SM-Project\
                          --tabname Amon\
                          --relname atmos\
                          --tmrange $trange \
			  --sourceid E3SM-2-0


    python ls3mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname TRENDY \
                          --modname E3SM \
                          --expname land-hist${expnam} \
                          --insname E3SM-Project\
                          --tabname Emon\
                          --relname land\
                          --tmrange $trange \
			  --sourceid E3SM-2-0
 
    #LImon
    python ls3mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname TRENDY \
                          --modname E3SM \
                          --expname land-hist${expnam} \
                          --insname E3SM-Project\
                          --tabname LImon\
                          --relname land\
                          --tmrange $trange \
                          --sourceid E3SM-2-0


done
