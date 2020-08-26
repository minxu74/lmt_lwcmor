#!/usr/bin/env python

'''
   This is the land model testbed (LMT) light-weighted cmorization toolkit
'''

import sympy
import numpy as np
import cf_units as cf
import netCDF4 as nc4
import os, sys, json, re
from shutil import copyfile

from mpi4py import MPI

import glob

def parse_rela(relastr):
    '''
       parse the relationship string and save the dimensions of variables
       to a dictionary and the relationship string only containing the 
       variables and operators
    '''

    import re
    
    #teststr="aaa(1)*bbb(0:10)"
    
    if '(' in relastr:
    
       p = re.compile(r"\w+\(")
       vnms=p.findall(relastr)

       print (vnms)
    
       p = re.compile(r"\(([0-9,:]+)\)")
       dims=p.findall(relastr)
    
       reladims={}
    
       for v, d in zip(vnms, dims):
           reladims[v.replace('(','')] = d
       print (dims)
    
       relaeqns=re.sub("\(*\)*", '', re.sub('\([0-9,:]+\)', '', relastr))

       return reladims, relaeqns
    
    else:
       return {}, relastr


def assign_varattrs(ncdf_var, VarMap):

    
    # hard-coded

    ncdf_var.setncattr('time'          ,  VarMap["attributes"]["time"])
    ncdf_var.setncattr('time_label'    ,  VarMap["attributes"]["time_label"])
    ncdf_var.setncattr('time_title'    ,  VarMap["attributes"]["time_title"])
    ncdf_var.setncattr('cell_measures' ,  VarMap["attributes"]["cell_measures"]) 
    ncdf_var.setncattr('cell_methods'  ,  VarMap["attributes"]["cell_methods"]) 



    # from data requestion python API
    ncdf_var.setncattr('description'   ,  VarMap["attributes"]["description"]) 
    ncdf_var.setncattr('comment'       ,  VarMap["attributes"]["comment"])
    ncdf_var.setncattr('standard_name' ,  VarMap["attributes"]["standard_name"])
    ncdf_var.setncattr('long_name'     ,  VarMap["attributes"]["long_name"])
    ncdf_var.setncattr('title'         ,  VarMap["attributes"]["long_name"])
    ncdf_var.setncattr('id'            ,  VarMap["attributes"]["id"])
    ncdf_var.setncattr('variable_id'   ,  VarMap["attributes"]["variable_id"])
    ncdf_var.setncattr('frequency'     ,  VarMap["attributes"]["frequency"])
    ncdf_var.setncattr('positive'      ,  VarMap["attributes"]["positive"])
    ncdf_var.setncattr('modeling_realm',  VarMap["attributes"]["modeling_realm"])
    ncdf_var.setncattr('type'          ,  VarMap["attributes"]["type"])
    ncdf_var.setncattr('prov'          ,  VarMap["attributes"]["prov"])
    ncdf_var.setncattr('mipTable'      ,  VarMap["attributes"]["miptable"])
    ncdf_var.setncattr('positive'      ,  VarMap["attributes"]["positive"])
    ncdf_var.setncattr('units'         ,  VarMap["attributes"]["units"])

    if ncdf_var.dtype == np.single: 
        ncdf_var.setncattr('_FillValue'    , np.float32(1.e20))
        ncdf_var.setncattr('missing_value' , np.float32(1.e20))
    elif ncdf_var.dtype == np.double:
        ncdf_var.setncattr('_FillValue'    , np.float64(1.e20))
        ncdf_var.setncattr('missing_value' , np.float64(1.e20))

    elif ncdf_var.dtype == np.int16:
        ncdf_var.setncattr('_FillValue'    , np.int16(-9999))
        ncdf_var.setncattr('missing_value' , np.int16(-9999))



    return


#mipname #realname

def assign_gblattrs(ncdf_gbl, VarMap):
    """
    #x-Conventions xum
    #x-activity_id xum
    #-creation_date xum
    #-data_specs_version xum
    #x-experiment xum
    #x-experiment_id xum
    #x-forcing_index xum
    #x-frequency xum
    #x-further_info_url xum
    #-grid xum
    #-grid_label xum
    #x-initialization_index xum
    #x-institution xum
    #x-institution_id xum
    #x-license xum
    #x-mip_era xum
    #x-nominal_resolution xum
    #x-physics_index xum
    #-product xum
    #x-realization_index xum
    #-realm xum
    #-source xum
    #-source_id xum
    #-source_type xum
    #x-sub_experiment xum
    #x-sub_experiment_id xum
    #x-table_id xum
    #x-tracking_id xum
    #-variable_id xum
    #-variant_label xum
    """

    import load_CMIP6CV as cv
    import uuid

    # input: VarMap

    mipcv = cv.cmip6cv()

    UserInput["experiment_id"        ] = '1pctCO2-bgc'
    UserInput["institution_id"       ] = 'RUBISCO'
    UserInput["grid_label"           ] = "gr"
    UserInput["nominal_resolution"   ] = "100 km"
    UserInput["further_info_url"     ] = "https://bgc-feedbacks.org"
    UserInput["creation_date"        ] = "2020-06-11"
    UserInput["source_component"     ] = ["atmos", "atmosChem", "land", "ocean", "ocnBgchem", "seaIce"]
    UserInput["source_id"            ] = "E3SM-1-1"
    UserInput["forcing_index"        ] = 1
    UserInput["initialization_index" ] = 1
    UserInput["physics_index"        ] = 1
    UserInput["realization_index"    ] = 1
    UserInput["variant_info"         ] = "p = 1, for CTC with CNP"
    UserInput["branch_method"        ] = "50 yrs spinup on DOE Compy from the 969 yrs piControl spinup on NERSC Edison " \
                                         "to account for the machine differences"
    UserInput["branch_time_in_child" ] = 0.0
    UserInput["branch_time_in_parent"] = 371935.0
    UserInput["parent_time_units"    ] = "days since 0001-01-01"
    # user the shell to get
    myuuid = str(uuid.uuid4())
    UserInput["tracking_id"          ] = "hdl:21.14100/" + myuuid

    # from dreq 
    DreqInput["frequency"            ] = VarMap["attributes"]["frequency"]
    DreqInput["realm"                ] = VarMap["attributes"]["modeling_realm"]
    DreqInput["table_id"             ] = VarMap["attributes"]["miptable"]
    DreqInput["variable_id"          ] = VarMap["cmvar"]

    mipcv.GetGlobalAttributes(UserInput, DreqInput) 

    # sort ???
    for gkey in mipcv.ReqGblAttrs.keys():
        setattr(ncdf_gbl, gkey, ReqGblAttrs[gkey]) 

    for gkey in mipcv.OptGblAttrs.keys():
        setattr(ncdf_gbl, gkey, OptGblAttrs[gkey])





    #-ncdf_gbl.Conventions = "CF-1.7 CMIP-6.2"
    #-ncdf_gbl.mip_era = "CMIP6"
    #-ncdf_gbl.product = "model-output"
    #-ncdf_gbl.further_info_url = "https://www.bgc-feedbacks.org"

    #-ncdf_gbl.activity_id    = mip_filt
    #-ncdf_gbl.table_id       = tab_filt
    #-ncdf_gbl.experiment_id  = exp_name
    #-ncdf_gbl.institution_id = ins_name
    #-ncdf_gbl.variant_label  = ens_name

    #-ncdf_gbl.source_id      = mod_name

    #-ncdf_gbl.experiment     = cv.dict_exp['experiment_id'][exp_name]['experiment']
    #-ncdf_gbl.institution    = cv.dict_ins['institution_id'][ins_name]
    #-ncdf_gbl.license        = cv.dict_lic['license'][0].replace('<Your Centre Name>', ins_name).\
    #-                          replace('<some URL maintained by modeling group>', 'https://www.bgc-feedbacks.org')


    #-for frq in cv.dict_frq['frequency'].keys():
    #-    if frq in tab_id:
    #-       ncdf_gbl.frequency = frq
    #-       break

    #-nomgrid=[]
    #-for nom in cv.dict_nom['nominal_resolution']:
    #-    if 'degree' in nom.split(' ')[1]:
    #-       nomgrid.append(110)
    #-    else:
    #-       nomgrid.append(float(nom.split(' ')[0]))
    #-ncdf_gbl.nominal_resolution = cv.dict_nom['nominal_resolution'][nomgrid.index(min(nomgrid, key=lambda x:abs(x-gridsize)))]


    #-ncdf_gbl.realization_index     = ens_name.split('r')[1].split('i')[0]
    #-ncdf_gbl.initialization_index  = ens_name.split('i')[1].split('p')[0]
    #-ncdf_gbl.physics_index         = ens_name.split('p')[1].split('f')[0]
    #-ncdf_gbl.forcing_index         = ens_name.split('p')[1].split('f')[1]

    #-ncdf_gbl.sub_experiment = "none"
    #-ncdf_gbl.sub_experiment_id = "none"


    #-# user the shell to get
    #-myuuid = str(uuid.uuid1())
    #-ncdf_gbl.tracking_id = "hdl:21.14100/" + myuuid

    #-ncdf_gbl.grid = '0.5x0.5 degree'
    #-ncdf_gbl.grid_label = 'gn'
    #-ncdf_gbl.realm = "land"



    #-src_str = cv.dict_sid['source_id'][mod_name]['source_id']+': (' +str(dict_sid['source_id'][mod_name]['release_year']) + '): '
    #-mykeys = dict_sid['source_id'][mod_name]['model_component'].keys()
    #-mykeys.sort()

    #-for mykey in keys:
    #-    if dict_sid['source_id'][mod_name]['model_component'][mykey]['description'] != 'none':
    #-        src_str += mykey + ': ' + dict_sid['source_id'][mod_name]['model_component'][mykey]['description'] + '; '


    #-ncdf_gbl.source = src_str
    #-ncdf_gbl.source_id = mod_name
    #-ncdf_gbl.source_type = "LAND"

    #-ncdf_gbl.variable_id = var_name

    #-ncdf_gbl.variant_label = "Same land model configuration, including representation of land cover, land use, \
    #-                          and land management, as used in coupled CMIP6 historical simulation with all applicable land-use features active. \
    #-                          Shared with LS3MIP and LUMIP"
    #-

    #-#options
    #-ncdf_gbl.model_doi_url = "http://dx.doi.org/10.11578/E3SM/dc.20180418.36"
    #-ncdf_gbl.contact = "forrest@climatemodeling.org"



    #-ncdf_gbl.parent_activity_id = "no parent" ;
    #-ncdf_gbl.parent_experiment_id = "no parent" ;
    #-ncdf_gbl.parent_mip_era = "no parent" ;
    #-ncdf_gbl.parent_source_id = "no parent" ;
    #-ncdf_gbl.parent_time_units = "no parent" ;
    #-ncdf_gbl.parent_variant_label = "no parent" ;


    return


def search_cmorlist(mip_filt, tab_filt):
    """
       Get a cmorlist from a MIP name
    """
    from dreqPy import dreq
    from dreqPy.extensions import collect


    dq = dreq.loadDreq()
    collect.add( dq )
    mip = dq.coll['mip'].items[0]._labelDict[mip_filt]

    mvars = mip._get__CMORvar()
    expts = mip._get__expt()

    MIP_Expts = set([ expt for expt in dq.coll['experiment'].items
                    if expt.uid in expts ])
    MIP_Expts = list(MIP_Expts)

    # Find the 'CMORvar'-items, associated with these requestVars:

    cmorlist = {}
    mvarlist = {}
    for cmvar in dq.coll['CMORvar'].items:
        if cmvar.label == 'lwsnl':
           print (cmvar.mipTable, tab_filt)
        if cmvar.uid in mvars and tab_filt in cmvar.mipTable:
           cmorlist[cmvar.label] = cmvar
           mvarlist[cmvar.label] = dq.inx.uid[cmvar.vid] 
           

    return cmorlist, mvarlist
            

def define_filename(var_name, rel_name, mod_name, exp_name, ens_grid, ens_time, ens_name='r1i1p1f1'):
    # mrro_Lmon_CESM2_historical_r1i1p1f1_gn_185001-201412.n
    ens_file= var_name + '_' + rel_name + '_' + mod_name + '_' + exp_name + '_' + ens_name + '_' + ens_grid + '_' + ens_time + '.nc'

    return ens_file

def define_dirlevel(var_name, institid, exp_name, ens_grid, date_ver):

    var_dirn = root_dir + '/' + mip_name + '/' + institid + '/' + mod_name + '/' + \
               exp_name + '/' + ens_name + '/' + miptable + '/' + var_name + '/' + ens_grid + '/' + data_ver + '/'
    # ~/CMIP6/LS3MIP/NCAR/CESM2/land-hist/r1i1p1f1/Lmon/var/gn/ver

    return var_dirn





comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


mip_filt='CMIP'
tab_filt='Lmon'

cmorlist={}
mvarlist={}

cmorlist, mvarlist = search_cmorlist(mip_filt, tab_filt)


print (cmorlist.keys())


ens_grid = 'gn'
ens_time = '000101-015012'



mod_name = 'E3SM'
exp_name = '1pctCO2-bgc'
mip_name = 'C4MIP'
ins_name = 'RUBISCO'
tab_name = 'Amon'
rel_name = 'atmos'

datadir="/global/cscratch1/sd/minxu/data/20191123.CO21PCTRAD_RUBISCO_CNPCTC20TR_OIBGC.I1900.ne30_oECv3.compy/rgr/atm/"
cmordir="/global/homes/m/minxu/scratch/data/mytest/"

timeshift=1899


if rank == 0:
    with open("cmorjson/jsonfiles/cmor/{}_{}_{}.json".format(mip_name, tab_name, mod_name)) as rdf:
        cmor_json = json.load(rdf)
    cmordicts = [ cm for cm in cmor_json["variables"] if int(cm["confidence"]) >= 90.0 ]

    
    # get all variable definitons

    for jsnf in glob.glob("cmorjson/jsonfiles/{}/*".format(mod_name.lower())):
        if "mon" in jsnf and rel_name in jsnf:
            with open ("cmorjson/jsonfiles/{}/".format(mod_name.lower()) + jsnf) as rdf:
                 modeljson = json.load(rdf)
            compdicts = modeljson["variables"]
else:
    cmordicts = None
    compdicts = None
cmordicts = comm.bcast(cmordicts, root=0)
compdicts = comm.bcast(cmordicts, root=0)

sys.exit()


#-
#-
#-varsneeded=[]
#-for c in cmor_json:
#-    print (c['cmvar'], c['relationship'])
#-    if c['relationship'].strip() != '':
#-       temp=re.sub('[+\-\*\/]',',',re.sub("\(*\)*", '', re.sub('\([0-9,:]+\)', '', c['relationship'])))
#-
#-       print ('xxx', re.sub("\(*\)*", '', re.sub('\([0-9,:]+\)', '', c['relationship'])))
#-       if ',' in temp:
#-          print (temp, temp.split(','))
#-          for t in temp.split(','):
#-              if not t.replace('.','',1).isdigit():
#-                 varsneeded.append(t.strip())
#-       else:
#-          varsneeded.append(temp)
#-#this list should be used by the workflow
#-print (list(set(varsneeded)))


# parallelism 
nseg = max(int(len(cmordicts)/size), 1)
nbgn = rank * nseg
nend = min(nbgn + nseg, len(cmordicts))


print (rank, nbgn, nend, size, nseg)
sys.exit()

for cmordict in cmordicts[nbgn:nend]:

#for cmordict in cmordicts:
    if int(cmordict['confidence']) < 80 or not cmordict['cmvar'] in cmorlist.keys():
       continue

    cmortask = ['copy', 'rename']
    tag_units = cf.Unit(cmordict['units'])

    print ('xxx', cmordict['units'], cmordict['cmvar'])

    if len(cmordict['_children']) > 1:

       var_units=[c['units'] for c in cmordict['_children']]

       print (var_units)

       if '+' in cmordict['relationship'] or '-' in cmordict['relationship']:

          unitset = list(set(var_units))
          print (cmordict['relationship'])
          print (unitset[0])

          if len(unitset) == 1:

             try:
                print( unitset[0] == 'mm/s' )

                if unitset[0] == 'mm/s':
                   src_units = cf.Unit('kg m-2 s-1')
                else:
                   src_units = cf.Unit(unitset[0])
             except:
                if 'C' in unitset[0]: 
                   src_units = cf.Unit(unitset[0].replace('C',''))
                elif 'N' in unitset[0]:
                   src_units = cf.Unit(unitset[0].replace('N',''))
                else:
                   print ("cannot convert", unitset)
               
          else:
             print ("unit error, cannot add and substract")

     
       elif '*' in cmordict['relationship']:
           for i, u in enumerate(var_units):
               if i == 0: 
                  src_units = cf.Unit(u)
               else:
                  src_units *= cf.Unit(u)
       
       else:

            print ("unit error other operators")

           
    elif len(cmordict['_children']) == 1:
       print(cmordict['_children'][0]["units"])
       if cmordict['_children'][0]["units"] == 'proportion' or cmordict['_children'][0]["units"] == 'unitless'  or \
          cmordict['_children'][0]["units"] == 'none':
          src_units = cf.Unit('1')

       elif cmordict['_children'][0]["units"] == 'mm/s':
          src_units = cf.Unit('kg m-2 s-1')

       else:
          try:
             src_units = cf.Unit(cmordict['_children'][0]['units'])
          except:
             if 'C' in cmordict['_children'][0]['units']: 
                src_units = cf.Unit(cmordict['_children'][0]['units'].replace('C',''))
             elif 'N' in cmordict['_children'][0]['units']:
                src_units = cf.Unit(cmordict['_children'][0]['units'].replace('N',''))
             else:
                print ("cannot convert")
               
       #check if the units is convertible after calculations

    print (tag_units, '---', src_units)
    if tag_units.is_convertible(src_units):
       if tag_units.convert(1., src_units) != 1.:
          cmortask.append("units")
       print(tag_units.convert(1., src_units))
    
    else:
       print("units are not convertible\n")
       sys.exit()



    print (src_units.convert(1., tag_units))

    print (cmortask)

    if cmordict['relationship'] != '':
       #preprocessing
       cmoreq = cmordict['relationship']

       cmdims, cmrela = parse_rela(cmoreq)

       print (cmoreq)
       print (cmdims)
       print (cmrela)

       if cmdims == {}: 

          rel = sympy.sympify(cmoreq)
          xxx = tuple(rel.free_symbols)
          fuc = sympy.lambdify(xxx, rel, 'numpy')

          print(type(rel.free_symbols))
          print(rel.free_symbols)
          print(rel)

          xlist=[]

          #QFLX_000101_015012.nc
          for mvar in xxx:
              with nc4.Dataset(datadir + '/' + str(mvar)+"_000101_015012.nc", "r") as f:
                   if len(f.variables[str(mvar)].dimensions) >=4: 
                      print ("huge array")
                      xlist.append(f.variables[str(mvar)])
                   else:
                      xlist.append(f.variables[str(mvar)][...])
                   vardimnm = f.variables[str(mvar)].dimensions
                   vardimno = f.variables[str(mvar)].shape

          #--newvar = fuc(*xlist)

       else:
          rel = sympy.sympify(cmrela)
          xxx = tuple(rel.free_symbols)
          fuc = sympy.lambdify(xxx, rel, 'numpy')
          xlist=[]
          klist=[]

          print (type(rel.free_symbols))
          for mvar in xxx:

              cvar = str(mvar)
              print (cvar)

              with nc4.Dataset(datadir + cvar + "_185001_201412.nc", "r") as f:
                  if cvar in cmdims.keys():
                     if ':' not in cmdims[cvar]: 
                        xlist.append(f.variables[cvar][:,int(cmdims[cvar]),:,:])
                     else:
                        kbgn = int(cmdims[cvar].split(':')[0])
                        kend = int(cmdims[cvar].split(':')[1])

                        varsum = np.sum(f.variables[cvar][:,kbgn:kend+1,:,:], axis=1)
                        xlist.append(varsum)

                        #-for k in range(kbgn, kend+1):
                        #-    klist.append(f.variables[cvar][:,k,:,:])

                     vardimnm = (f.variables[cvar].dimensions[0],) 
                     for vn in f.variables[cvar].dimensions[2:]:
                         vardimnm = vardimnm + (vn,)
                     vardimno = (f.variables[cvar].shape[0],)
                     for vn in f.variables[cvar].shape[2:]:
                         vardimno = vardimno + (vn,)

                            
                  else:

                     if len(f.variables[str(mvar)].dimensions) >=4: 
                        print ("huge array")
                        xlist.append(f.variables[str(mvar)])
                     else:
                        xlist.append(f.variables[str(mvar)][...])
                     vardimnm = f.variables[str(mvar)].dimensions 
                     vardimno = f.variables[str(mvar)].shape
                       


          print (vardimnm)
          print (vardimno)

          #--for x in xlist:
          #--    print (x.shape)
          #--for k in klist:
          #--    print (k.shape)

          #--if len(vardimnm) >=4:
          #--   klist = []
          #--   for it in range(vardimno[0]):
          #--       for xv in xlist:
          #--           klist.append(xv[it])
          #--       newvar = fuc(*klist)

          #--else:
          #--   newvar = fuc(*xlist)
          #-if not klist:

          #-   print ("no complex")
          #-   newvar = fuc(*xlist)

          #-else:

          #-   for k, ks in enumerate(klist):
          #-       tlist = []
          #-       tlist = xlist.copy()
          #-       print ('k=', k)
          #-       print (len(tlist))
          #-       print (len(xlist))
          #-       tlist.append(ks)
          #-       print (len(tlist))
          #-       if k == 0:
          #-          newvar = fuc(*tlist)
          #-       else:
          #-          newvar = newvar + fuc(*tlist)
                
       if 'units'  in cmortask:
          idx = cmortask.index('units')
          cmortask.insert(idx, 'compute')
       else:
          cmortask.append('compute')

                      
       #-print (newvar.shape)

       # the relationship if level




    var_name = cmordict['cmvar']
    rel_name = cmorlist[var_name].mipTable

    cmorfn = cmordir + define_filename(var_name, rel_name, mod_name, exp_name, ens_grid, ens_time, ens_name='r1i1p1f1')

    origfn = datadir + cmordict['_children'][0]['cmvar'] + "_000101_015012.nc"
    origvn = cmordict['_children'][0]['cmvar']

    for task in cmortask:
         if type(task).__name__ == 'list':
            for subtask in task:
                 print (subtask)
         else:
            print (cmordict['cmvar'], task)

         if task == 'copy':

            with nc4.Dataset(origfn, "r") as forig, nc4.Dataset(cmorfn, "w") as fcmor:

                # copy global attributes all at once via dictionary
                #fcmor.setncatts(forig.__dict__)
                # copy dimensions
                for name, dimension in forig.dimensions.items():
                    if name in vardimnm or name in ['nb', 'nbnd']:
                       fcmor.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))


            #copyfile(origfn, cmorfn)

         #with nc4.Dataset(cmorfn, "r+") as fcmor:
         #   if task == 'rename':
         #        fcmor.renameVariable(origvn,  cmordict['cmvar']) 
         #        fcmor[cmordict['cmvar']].setncattr('long_name',  cmordict['longname'])

                #mxu implemented cmip6 requirements

                fcmor.createVariable(cmordict['cmvar'], forig.variables[origvn].datatype, vardimnm, complevel=2, shuffle=True, fill_value=1.e20)

                # copy all file data except for the excluded
                for name, variable in forig.variables.items():
                    if name != origvn:

                        extdim=True
                        for vd in variable.dimensions:
                            if vd not in vardimnm:
                               extdim = False
                               break
                        if name == "time_bnds" or name == "time_bounds":
                            extdim = True
                        if extdim:
                           print (name, variable.datatype,  variable.dimensions)
                           x = fcmor.createVariable(name, variable.datatype, variable.dimensions)
                           fcmor[name].setncatts(forig[name].__dict__)
                           fcmor[name][:] = forig[name][:]
                           # copy variable attributes all at once via dictionary
                if timeshift != 0:
                    #time_bounds=time_bounds+$(($difyear*(-365))); time=time_bounds(:,1);

                    print ("TIME SHIFTING !!!!!!")

                    if 'time_bounds' in fcmor.variables.keys(): 
                       fcmor['time_bounds'][...] = fcmor['time_bounds'][...] - timeshift * 365
                       fcmor['time'][:] = fcmor['time_bounds'][:,1]
                    elif 'time_bnds' in fcmor.variables.keys():
                       fcmor['time_bnds'][...] = fcmor['time_bnds'][...] - timeshift * 365
                       fcmor['time'][:] = fcmor['time_bnds'][:,1]

                    else:
                       print ('cannot time shifting due to no time variables')



         if task == 'compute':
            with nc4.Dataset(origfn, "r") as forig, nc4.Dataset(cmorfn, "r+") as fcmor:
              print ('computing {}'.format(cmordict['cmvar']))
              #-fcmor.variables[cmordict['cmvar']][...] = newvar[...]

              if len(vardimnm) >=4:

                 print ("huge data", len(xlist))
                 for it in range(vardimno[0]):
                     klist = []
                     for xv in xlist:
                         klist.append(xv[it])
                     newvar = fuc(*klist)
                     fcmor.variables[cmordict['cmvar']][it] = newvar

              else:
                 newvar = fuc(*xlist)
                 fcmor.variables[cmordict['cmvar']][...]= newvar

            

         if task == 'units':
            with nc4.Dataset(cmorfn, "r+") as fcmor:
               print ("change units", src_units.convert(1., tag_units))
               fcmor.variables[cmordict['cmvar']][...] = src_units.convert(fcmor.variables[cmordict['cmvar']][...], tag_units)
               #fcmor.variables[cmordict['cmvar']].setncattr('units',      cmordict['units'])
              #-fcmor[cmordict['cmvar']].setncattr('_FillValue', 1.e33)


         # writing meta data
         with nc4.Dataset(cmorfn, "r+") as fcmor:
            fcmor.variables[cmordict['cmvar']].setncattr('units',      cmordict['units'])
            fcmor.variables[cmordict['cmvar']].setncattr('missingvalue', 1.0e20)
            assign_varattrs(fcmor.variables[cmordict['cmvar']], cmordict['cmvar'], cmorlist, mvarlist)
            #-assign_gblattrs(fcmor)
            

sys.exit()




print(tag_units.is_convertible(src_units))

print(tag_units.convert(1., src_units))

print (tag_units, src_units)




f = nc4.Dataset("/global/homes/m/minxu/scratch/ILAMB_WCYCLE_20190319/cbgc_eca_prc/20190308.BDRD_CNPECACNT_20TR/rgr/atm/atm/QFLX_190001_200612.nc", "r")


oldvar = f.variables[cmordict['_children'][0]['cmvar']]

newvar = tag_units.convert(f.variables[cmordict['_children'][0]['cmvar']][...], src_units)

print (oldvar[...] == newvar)

newvar.longname = cmordict['longname']
newvar.units = cmordict['units']



print (oldvar.__dict__)
print (newvar.__dict__)
print (newvar.shape)



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='Light-weighted cmorize tool')


    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                               help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                           const=sum, default=max,
                           help='sum the integers (default: find the max)')

    parser.add_argument('--sum', dest='accumulate', action='store_const',
                           const=sum, default=max,
                           help='sum the integers (default: find the max)')

    parser.add_argument('--sum', dest='accumulate', action='store_const',
                           const=sum, default=max,
                           help='sum the integers (default: find the max)')
#-toexclude = ['ExcludeVar1', 'ExcludeVar2']
#-
#-with netCDF4.Dataset("in.nc") as src, netCDF4.Dataset("out.nc", "w") as dst:
#-    # copy global attributes all at once via dictionary
#-    dst.setncatts(src.__dict__)
#-    # copy dimensions
#-    for name, dimension in src.dimensions.items():
#-        dst.createDimension(
#-            name, (len(dimension) if not dimension.isunlimited() else None))
#-    # copy all file data except for the excluded
#-    for name, variable in src.variables.items():
#-        if name not in toexclude:
#-            x = dst.createVariable(name, variable.datatype, variable.dimensions)
#-            dst[name][:] = src[name][:]
#-            # copy variable attributes all at once via dictionary
#-            dst[name].setncatts(src[name].__dict__)
#-
#-with netCDF4.Dataset(file1) as src, netCDF4.Dataset(file2) as dst:
#-
#-  for name, dimension in src.dimensions.iteritems():
#-    dst.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
#-
#-  for name, variable in src.variables.iteritems():
#-
#-    # take out the variable you don't want
#-    if name == 'some_variable': 
#-      continue
#-
#-    x = dst.createVariable(name, variable.datatype, variable.dimensions)
#-    dst.variables[x][:] = src.variables[x][:]
