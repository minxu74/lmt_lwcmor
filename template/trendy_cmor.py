#!/usr/bin/env python

import subprocess
import datetime
import uuid
import numpy as np
import json
import glob, sys
from lwcmorlib import lwcmor
import argparse


"""
    This is for cmorization LS3MIP
    The hard-coded part can be a config file in futute
"""

parser = argparse.ArgumentParser()

parser.add_argument("--datadir",  help="the directory contains the original data", required=True)
parser.add_argument("--cmordir",  help="the directory contains the processed cmorized data", required=True)
parser.add_argument("--mipname",  help="MIP name (C4MIP, LS3MIP ...)", required=True)
parser.add_argument("--modname",  help="model name (E3SM, CESM2 ...)")
parser.add_argument("--expname",  help="experiment name (historical, land-hist ...)", required=True)
parser.add_argument("--insname",  help="institution name (RUBISCO, E3SM ...)")
parser.add_argument("--tabname",  help="table name (Lmon, Amon, Emon ...)", required=True)
parser.add_argument("--relname",  help="real name (land, atmos ...)", required=True)
parser.add_argument("--tmrange",  help="time range (185001-201412)", required=True)
parser.add_argument("--sourceid", help="source id (E3SM-1-1)", required=True)

args = parser.parse_args()

data_dir = args.datadir
cmor_dir = args.cmordir

mip_name = args.mipname
if not args.modname:
   mod_name = args.modname
else:
   mod_name = 'E3SM' 

exp_name = args.expname

if not args.insname:
   ins_name = args.insname
else:
   ins_name = 'RUBISCO'

tab_name = args.tabname
rel_name = args.relname

tm_range = args.tmrange

source_id = args.sourceid


timeshift= 0

# areacella_fx_CESM2_esm-hist_r1i1p1f1_gn.nc
UserInput = {}
#UserInput["mip_name"             ] = "CMIP"
#UserInput["tab_name"             ] = "fx"
#UserInput["rel_name"             ] = "atmos"
UserInput["mip_name"             ] = "trendy2023"
UserInput["tab_name"             ] = tab_name
UserInput["rel_name"             ] = rel_name
UserInput["mod_name"             ] = mod_name
UserInput["time_range"           ] = tm_range

UserInput["ModelRltDir"          ] = data_dir
UserInput["CmorRltDir"           ] = cmor_dir

UserInput["timeshift"            ] = timeshift

UserInput["experiment_id"        ] = exp_name
UserInput["institution_id"       ] = ins_name
UserInput["grid_label"           ] = "gn"
UserInput["nominal_resolution"   ] = "50 km"
UserInput["further_info_url"     ] = "https://e3sm.org"
UserInput["creation_date"        ] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
UserInput["source_component"     ] = ["atmos", "atmosChem", "land", "ocean", "ocnBgchem", "seaIce"]
UserInput["source_id"            ] = source_id
UserInput["forcing_index"        ] = np.int32(1)
UserInput["initialization_index" ] = np.int32(1)
UserInput["physics_index"        ] = np.int32(1)
UserInput["realization_index"    ] = np.int32(1)
UserInput["variant_info"         ] = "S1 experiment simulation"
UserInput["branch_method"        ] = "no parent" 
UserInput["branch_time_in_child" ] = 674885.
UserInput["branch_time_in_parent"] = 0.


UserInput["parent_activity_id"   ] = "no parent"
UserInput["parent_experiment_id" ] = "no parent"
UserInput["parent_mip_era"       ] = "no parent"
UserInput["parent_source_id"     ] = "no parent"
UserInput["parent_source_id"     ] = "no parent"
UserInput["parent_time_units"    ] = "no parent"
UserInput["parent_variant_label" ] = "no parent"


UserInput["data_version"         ] = "v20230831"



# user the shell to get
myuuid = str(uuid.uuid4())
UserInput["tracking_id"          ] = "hdl:21.14100/" + myuuid
UserInput["calendar"             ] = "noleap"
UserInput["timeunits"            ] = "days since 0001-01-01 00:00:00"

# get the cmorjson version
jsonversion = subprocess.check_output("cd cmorjson && git describe --always", shell=True)

if tab_name == 'fx':
   with open("cmorjson/jsonfiles/cmor/{}_{}_{}.json".format("CMIP", UserInput["tab_name"], UserInput["mod_name"])) as rdf:
       cmor_json = json.load(rdf)
else:
   with open("cmorjson/jsonfiles/cmor/{}_{}_{}.json".format(UserInput["mip_name"], UserInput["tab_name"], UserInput["mod_name"])) as rdf:
       cmor_json = json.load(rdf)
cmordicts = [ cm for cm in cmor_json["variables"] if int(cm["confidence"]) >= 90.0 ]

if mip_name == "LS3MIP" or mip_name == "trendy2023":
    xrel_name = "land"
else:
    xrel_name = rel_name

# get all variable definitons
for jsnf in glob.glob("cmorjson/jsonfiles/{}/*".format(mod_name.lower())):
    if "mon" in jsnf and xrel_name in jsnf:
        print (jsnf)
        with open (jsnf) as rdf:
             modeljson = json.load(rdf)
        compdicts = modeljson["variables"]



fmiss = open("missingvars_"+tab_name+"_"+exp_name+".txt", "w")


for i, cmordict in enumerate(cmordicts):

    #if i > 5:
    #    break

       #incVars = ["gpp", "tas", "areacella", "mrtws"]
       #incVars = ["rsds", "rlds", "rsus", "rlus", "hfss", "hfls", "evspsbl", "mrro", "mrso", "mrtws", "tran", "evspsblsoi", "evspsblveg", "mrsos", "nbp", 
       #           "gpp", "ra", "rh", "lai", "tas", "ps", "pr", "prsn", "snc", "snw", "snd", "tsl", "sftlf"]
       #incVars = ["mrsos" , "snw", "snd"]
       incVars = ["tas", "pr", "rsds", "mrso", "mrro", "evapotrans",
                   "cVeg", "cLitter", "cSoil", "cProduct", 
                   "gpp", "ra", "npp", "rh", "fFire", "fLuc", "nbp", 
                   "landCoverFrac", "oceanCoverFrac", "burntArea", "lai",
                   "cLeaf", "cWood", "cRoot", "cCwd", "cSoilpools", 
                   "fVegLitter", "fLeafLitter", "fWoodLitter", 
"fRootLittera", "fLitterSoil", "fVegSoil", "rhpool", "fAllocLeaf", 
"fAllocWood", "fAllocRoot", "fFireCveg", "fFireLittera", "fFireCsoil", 
"tsl", "msl", "evspsblveg"] 


       #if cmordict["cmvar"] != "gpp" and cmordict["cmvar"] != "tas" and cmordict["cmvar"] != "areacella" and cmordict["cmvar"] != "mrtws":
       #    continue

       if cmordict["cmvar"] not in incVars:
           continue

    try:
       lwc = lwcmor(cmordict, compdicts, UserInput)
       lwc.UnitsConversion()
       lwc.CmorvarCompute()


       if len(lwc.missingvars) > 0:
          fmiss.write(','.join(lwc.missingvars)+",")

       variant_label = "r{}i{}p{}f{}".format(UserInput["realization_index"], UserInput["initialization_index"], UserInput["physics_index"], UserInput["forcing_index"])
       optgblattrs={}
       optgblattrs.update({"further_info_url":"https://furtherinfo.es-doc.org/CMIP6.{}.{}.{}.{}.{}".format(\
                          UserInput["institution_id"],UserInput["source_id"],UserInput["experiment_id"],'none', variant_label)})
       #CMIP6.NCAR.CESM2.1pctCO2-bgc.none.r1i1p1f1"
       optgblattrs.update({"model_doi_url":"https://doi.org/10.11578/E3SM/dc.20180418.36"})
       optgblattrs.update({"relationship":lwc.cmrvar["relationship"]})
       optgblattrs.update({"sympyinput": lwc.cmrela})
       optgblattrs.update({"CmorJsonVersion": jsonversion.strip()})
       optgblattrs.update({"contact": "Forrest M. Hoffman <forrest@ornl.gov>"})

       lwc.GetAttributes(optgblattrs, {})
       lwc.Cmorization()
       #print (lwc)
       #print (lwc.cmfnam, lwc.cmdnam)
    except:
       print ('Error in', cmordict["cmvar"], sys.exc_info()[0])


