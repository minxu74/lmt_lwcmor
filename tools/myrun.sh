#!/usr/bin/env bash

#topdir=/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/
#topdir=/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190320_gswpv2_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/
#topdir="/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190224_princeton_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC/"
topdir="/global/homes/m/minxu/scratch/LS3MIP_CMOR/20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC"
expnam=''

python ls3mip_cmor.py --datadir $topdir/org/ \
	              --cmordir $topdir/cmor_compress/ \
		      --mipname LS3MIP \
		      --modname E3SM \
		      --expname land-hist${expnam} \
		      --insname RUBISCO\
		      --tabname fx\
		      --relname atmos\
                      --tmrange 185001-201412\
		      --sourceid E3SM-1-1

#-exit

python ls3mip_cmor.py --datadir $topdir/org/lnd/ \
	              --cmordir $topdir/cmor_compress/ \
		      --mipname LS3MIP \
		      --modname E3SM \
		      --expname land-hist${expnam} \
		      --insname RUBISCO\
		      --tabname Amon\
		      --relname atmos\
                      --tmrange 185001-201412\
		      --sourceid E3SM-1-1

python ls3mip_cmor.py --datadir $topdir/org/lnd/ \
	              --cmordir $topdir/cmor_compress/ \
		      --mipname LS3MIP \
		      --modname E3SM \
		      --expname land-hist${expnam} \
		      --insname RUBISCO\
		      --tabname Lmon\
		      --relname land\
                      --tmrange 185001-201412\
		      --sourceid E3SM-1-1 

python ls3mip_cmor.py --datadir $topdir/org/lnd/ \
	              --cmordir $topdir/cmor_compress/ \
		      --mipname LS3MIP \
		      --modname E3SM \
		      --expname land-hist${expnam} \
		      --insname RUBISCO\
		      --tabname Emon\
		      --relname land\
                      --tmrange 185001-201412\
		      --sourceid E3SM-1-1
