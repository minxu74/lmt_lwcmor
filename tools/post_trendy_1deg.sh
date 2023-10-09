#!/usr/bin/env bash


export HDF5_USE_FILE_LOCKING=FALSE

prefix=/qfs/people/xumi699/compyfs/scratch/forxj/trendy_postprocess/new_fv09

topdirs=(\
20230926_S3_f09_f09_ICB20TRCNPRDCTCBC)
expnams=(S3)

#-topdirs=(\
#-20230830_S3_f09_f09_ICB20TRCNPRDCTCBC \
#-20230830_S2_f09_f09_ICB20TRCNPRDCTCBC \
#-20230830_S1_f09_f09_ICB20TRCNPRDCTCBC \
#-20230830_S0_f09_f09_ICB20TRCNPRDCTCBC)

wkdir=$(pwd)

k=0
for topdir in "${topdirs[@]}"; do
    echo $topdir
    
    expnam="-${expnams[$k]}"
    trange=170001-202212

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

    if [[ ! -d "$prefix/$topdir/cmor" ]]; then
       mkdir -p $prefix/$topdir/cmor
    fi

    cd $prefix/$topdir/org/lnd/Amon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/Lmon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/Emon  && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir
    cd $prefix/$topdir/org/lnd/LImon && ln -sf ../*.nc . && ln -sf $wkdir/bounds.nc . && cd $wkdir


    if [[ $expnam == "-S3" ]]; then
        python trendy_cmor.py --datadir $prefix/$topdir/org/ \
                              --cmordir $prefix/$topdir/cmor/ \
                              --mipname TRENDY \
                              --modname E3SM \
                              --expname trendy${expnam} \
                              --insname E3SM-Project\
                              --tabname fx\
                              --relname land\
                              --tmrange $trange \
                              --sourceid E3SM-2-0
    fi


    #-python trendy_cmor.py --datadir $prefix/$topdir/org/lnd/ \
    #-                      --cmordir $prefix/$topdir/cmor/ \
    #-                      --mipname TRENDY \
    #-                      --modname E3SM \
    #-                      --expname trendy${expnam} \
    #-                      --insname E3SM-Project\
    #-                      --tabname Lmon\
    #-                      --relname land\
    #-                      --tmrange $trange \
    #-                      --sourceid E3SM-2-0



    #-python trendy_cmor.py --datadir $prefix/$topdir/org/lnd/ \
    #-                      --cmordir $prefix/$topdir/cmor/ \
    #-                      --mipname TRENDY \
    #-                      --modname E3SM \
    #-                      --expname trendy${expnam} \
    #-                      --insname E3SM-Project\
    #-                      --tabname Amon\
    #-                      --relname land\
    #-                      --tmrange $trange \
    #-    		  --sourceid E3SM-2-0

    #-python trendy_cmor.py --datadir $prefix/$topdir/org/lnd/ \
    #-                      --cmordir $prefix/$topdir/cmor/ \
    #-                      --mipname TRENDY \
    #-                      --modname E3SM \
    #-                      --expname trendy${expnam} \
    #-                      --insname E3SM-Project\
    #-                      --tabname Emon\
    #-                      --relname land\
    #-                      --tmrange $trange \
    #-                      --sourceid E3SM-2-0

    k=$((k+1))
done
