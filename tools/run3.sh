#1/usr/bin/env bash

python ls3mip_cmor.py --datadir /global/homes/m/minxu/scratch/LS3MIP_CMOR/20180910_princeton_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC/org/lnd/ \
	              --cmordir /global/homes/m/minxu/scratch/LS3MIP_CMOR/20180910_princeton_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC/cmor/ \
		      --mipname LS3MIP \
		      --modname E3SM \
		      --expname land-hist-princeton \
		      --insname RUBISCO\
		      --tabname Emon\
		      --relname land\
                      --tmrange 185001-201412\
		      --sourceid E3SM-1-1

