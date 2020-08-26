#!/usr/bin/env bash


echo "PrePare checking ..."

prefix=/global/homes/m/minxu/scratch/LS3MIP_CMOR_FX

topdirs=(\
20180910_cruncepv8_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20180910_princeton_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20190224_princeton_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC \
20190319_gspwv2_hcru_hcru_S6_hcru_hcru_ICB20TRCNPRDCTCBC \
20190320_gswpv2_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC \
20190328_cruncepv8_cn_hcru_hcru_S6_hcru_hcru_ICB20TRCNRDCTCBC) 

TablePath=/global/homes/m/minxu/MyGit/MySrc/CMOR/CMIP6_work/cmip6-cmor-tables/Tables/

for tdir in "${topdirs[@]}"; do
    echo "Checking the case $tdir"

    /bin/rm -f $prefix/$tdir/cmor/log_cmorchecker_$tdir.txt
    /bin/rm -f $prefix/$tdir/cmor/checksum_$tdir.md5
    for ncf in $(find $prefix/$tdir/cmor -name '*.nc' ); do
        NumErrors=`PrePARE --table-path $TablePath $ncf | grep -i "Number of file with error(s)" |cut -d ':' -f 2`
	CfcResult=`cfchecks $ncf | tail -3 | tr '\n' ' '`
	printf "%-90s| %s\n" $(basename -- $ncf) "PrePARE (ERRORS: $NumErrors) CFcheck ($CfcResult)" >> $prefix/$tdir/cmor/log_cmorchecker_$tdir.txt
	md5sum $ncf >> $prefix/$tdir/cmor/checksum_$tdir.md5
    done
done



