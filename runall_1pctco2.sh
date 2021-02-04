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

#topdirs=(\
#20200509.CO21PCTBGC_RUBISCO_CNCTC20TR_OIBGC.ne30_oECv3.compy  \
#20200802.CO21PCTFUL_RUBISCO_CNCTC20TR_OIBGC.ne30_oECv3.compy)
#topdirs=(20200509.CO21PCTCTL_RUBISCO_CNCTC1850_OIBGC.ne30_oECv3.compy)
topdirs=(20200622.CO21PCTRAD_RUBISCO_CNCTC20TR_OIBGC.ne30_oECv3.compy)

gridtype='rgr'

wkdir=$(pwd)

for topdir in "${topdirs[@]}"; do
    echo $topdir
    
    trange=000101-015012
    expnam=''
    if [[ "$topdir" == *"PCTBGC"* ]]; then
       expnam='-bgc'
    fi

    if [[ "$topdir" == *"PCTFUL"* ]]; then
       expnam=''
    fi
    if [[ "$topdir" == *"PCTRAD"* ]]; then
       expnam='-rad'
    fi

    echo $expnam $trange

    #if [[ ! -d "$prefix/$topdir/${gridtype}/lnd/Amon" ]]; then
    #   mkdir -p $prefix/$topdir/${gridtype}/lnd/Amon
    #fi

    if [[ ! -d "$prefix/$topdir/${gridtype}/atm/Amon" ]]; then
       mkdir -p $prefix/$topdir/${gridtype}/atm/Amon
    fi
    if [[ ! -d "$prefix/$topdir/${gridtype}/atm/Omon" ]]; then
       mkdir -p $prefix/$topdir/${gridtype}/atm/Omon
    fi

    if [[ ! -d "$prefix/$topdir/${gridtype}/lnd/Lmon" ]]; then
       mkdir -p $prefix/$topdir/${gridtype}/lnd/Lmon
    fi

    if [[ ! -d "$prefix/$topdir/${gridtype}/lnd/Emon" ]]; then
       mkdir -p $prefix/$topdir/${gridtype}/lnd/Emon
    fi

    if [[ ! -d "$prefix/$topdir/${gridtype}/lnd/LImon" ]]; then
       mkdir -p $prefix/$topdir/${gridtype}/lnd/LImon
    fi

    cd $prefix/$topdir/${gridtype}/atm/Amon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/${gridtype}/atm/Omon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/${gridtype}/lnd/Lmon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/${gridtype}/lnd/Emon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/${gridtype}/lnd/LImon && ln -sf ../*.nc . && ln -sf $wkdir/bounds/1pctco2/bounds.nc . && cd $wkdir

    mkdir -p $prefix/$topdir/cmor


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


    python c4mip_cmor.py --datadir $prefix/$topdir/$gridtype/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname Lmon\
                          --relname land\
                          --tmrange $trange \
                          --sourceid E3SM-1-1


    python c4mip_cmor.py --datadir $prefix/$topdir/$gridtype/atm/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname Amon\
                          --relname atmos\
                          --tmrange $trange \
			  --sourceid E3SM-1-1


    python c4mip_cmor.py --datadir $prefix/$topdir/$gridtype/lnd/ \
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
    python c4mip_cmor.py --datadir $prefix/$topdir/$gridtype/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname LImon\
                          --relname land\
                          --tmrange $trange \
                          --sourceid E3SM-1-1
    #ocean

    python c4mip_cmor.py --datadir $prefix/$topdir/$gridtype/atm/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname C4MIP \
                          --modname E3SM \
                          --expname 1pctCO2${expnam} \
                          --insname RUBISCO\
                          --tabname Omon\
                          --relname ocean\
                          --tmrange $trange \
                          --sourceid E3SM-1-1

done
