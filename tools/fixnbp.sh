#!/usr/bin/env bash

#topdir=/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/
#topdir=/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190320_gswpv2_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/
#topdir="/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190224_princeton_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/"
#topdir="/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC"
#expnam=''

prefix=/global/homes/m/minxu/scratch/LS3MIP_CMOR/

topdirs=(\
20180910_cruncepv8_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20180910_princeton_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20190224_princeton_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC \
20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20190320_gswpv2_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC \
20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC) 

#topdirs=(20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC) 

for topdir in "${topdirs[@]}"; do
    echo $topdir
    
    trange=185001-201412
    expnam=''
    if [[ "$topdir" == *"_cruncep"* ]]; then
       expnam='-cruNcep'
    fi

    if [[ "$topdir" == *"_princeton"* ]]; then
       expnam='-princeton'
    fi

    if [[ "$topdir" == *"_princeton_cn"* ]]; then
       trange=185001-201212
    fi
    echo $expnam $trange


    python ls3mip_cmor.py --datadir $prefix/$topdir/org/lnd/ \
                          --cmordir $prefix/$topdir/cmor/ \
                          --mipname LS3MIP \
                          --modname E3SM \
                          --expname land-hist${expnam} \
                          --insname RUBISCO\
                          --tabname Lmon\
                          --relname land\
                          --tmrange $trange \
			  --sourceid E3SM-1-1
done
