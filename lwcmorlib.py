#!/usr/bin/env python

import sympy, json, glob, re, os
import cf_units as cf
import numpy as np
import netCDF4 as nc4
from sympy.parsing.sympy_parser import parse_expr


class lwcmor():

    def __init__(self, cmorvmap, compdict, userinput):
       #validate the relationship

       cmoreq = cmorvmap['relationship']
       #cmoreq = "sum@levgrnd=SOILLIQ(1:15)#{'levgrnd'}"
       #cmoreq = "SOILLIQ[0:14]"
       #cmoreq = "SOILLIQ*levgrnd/sum(levgrnd)#{sum@levgrnd:0-14;levgrnd:inline=SOILLIQ}"
       #cmoreq = "hybi*PS+hyai*100000#{shape:time,ilev,lat,lon}"
       #cmoreq = "hybm*PS+hyam*float(100000)#{np.float32@;shape:time,lev,lat,lon}"

       print ('cmoreq', cmoreq)

       if '#' in cmoreq:
          cmtemp = cmoreq.split('#')[0].strip()
          cmdmns = cmoreq.split('#')[1].strip()
       else:
          cmtemp = cmoreq
          cmdmns = "{}"

       if "=" in cmtemp:
          cmrcnd = cmtemp.split('=')[0].strip()
          cmoreq = cmtemp.split('=')[1].strip()
       else:
          cmrcnd = ''
          cmoreq = cmtemp.split('=')[0].strip()

       if '[' in cmoreq:
          p = re.compile(r"\w+\[")
          vnms=p.findall(cmoreq)

          p = re.compile(r"\[([0-9,:]+)\]")
          dims=p.findall(cmoreq)

          cmdims={}

          for v, d in zip(vnms, dims):
              cmdims[v.replace('[','')] = d

          cmrela = re.sub("\[*\]*", '', re.sub('\[[0-9,:]+\]', '', cmoreq))

       else:
          cmdims = {}
          cmrela = cmoreq

       if cmrcnd != '':
          cmoper = cmrcnd.split('@')[0].strip()
          if len(cmrcnd.split('@')) == 2:
             cmopdm = cmrcnd.split('@')[1]
          else:
             cmopdm = None
       else:
          cmoper = None
          cmopdm = None
           
       if cmoper:
          cmrela = cmoper + "(" + cmrela + ")"
       print (cmorvmap["cmvar"], 'relationship1', cmrela)
       #cmrela = "dot(SOILLIQ, SOILICE)"

       print (cmdims, cmrela)

       syrela = sympy.sympify(cmrela)
       syvars = tuple(syrela.free_symbols)
       syfunc = sympy.lambdify(syvars, syrela, "numpy")

       #-xa = sympy.IndexedBase('SOILLIQ')
       #-n, m, j, i = sympy.symbols('n m j i', integer = True)

       #-xs = sympy.Sum(xa[:,m,:,:], (m,0,14))
       #-fn = sympy.lambdify(xa, xs)

       #-t1 = np.arange(30).reshape((5,3,2))
       #-t2 = np.arange(3).reshape((1,3,1))

       #-print (t1)
       #-print (t2)
       #-print (syfunc(*[t1, t2]))

       #-print (syvars, cmrela)

       # check if the variables are in the model outputs
       self.modvar = {}
       for vn in syvars:
           if vn.name in compdict.keys():
               self.modvar[vn.name] = compdict[vn.name]
       print (len(syvars))
       print (syvars)
       print (self.modvar.keys())
       if len(syvars) != len(self.modvar.keys()):
           raise ValueError

       self.cmrvar = cmorvmap

       self.syrela = syrela
       self.syvars = syvars
       self.syfunc = syfunc
       self.cmrela = cmrela
       self.cmdims = cmdims   # apply to individual variables
       self.cmdmns = cmdmns   # apply to all variables
       self.cmopdm = cmopdm   # reserved for future use

       #if len(cmdmns.split(',')) != len(cmdims.keys()):
       #   raise ValueError
       self.userinput = userinput

       self.cmfnam = None
       self.cmdnam = None
       self.cmvval = None
       self.cmvshp = None
       self.cmvdim = None

       self.missingvars = []



       #overwite
       print (self.cmrvar["dimensions"])
       print (self.cmrvar["cmvar"])
       if "time2" in self.cmrvar["dimensions"] and "time" not in self.cmrvar["dimensions"]:
           self.cmrvar["dimensions"]["time"] = self.cmrvar["dimensions"]["time2"]

       print (self.cmrvar["dimensions"]["shape"])
       if "time" in ','.join(self.cmrvar["dimensions"]["shape"]):
          self.cmrvar["dimensions"]["time"]["units"] = userinput["timeunits"]
          self.cmrvar["dimensions"]["time"]["calendar"] = userinput["calendar"]

       return


    def __FindModVarFile(self, sv):
        """
            Find the model varaible file
        """

        if sv.name+':inline' in self.cmdmns:
            p = re.compile(sv.name+r":inline=(\w+)")
            fname = p.findall(self.cmdmns)[0]
        elif 'plev19' in self.cmrvar["dimensions"]:   # standard pressure level
            fname = sv.name + "_plev19"
        elif 'fx' in self.cmrvar["attributes"]["miptable"]:
            fname = 'atm_area_landfrac'
        else:
            fname = sv.name
        fnabs = "{}/{}/{}_{}.nc".format(self.userinput["ModelRltDir"], self.cmrvar["attributes"]["miptable"], 
                                 fname, self.userinput["time_range"].replace('-','_'))

        return fnabs


    def __Map2npType(self, CmorType, ifv): 
        """
            Mapping the data types defined in the CMIP6 data request (CmorType) 
            to numpy types
        """

        if CmorType == 'real':
            DimVarType = np.float32
            Fill_Value = np.float32(1.0e20)
        elif CmorType == 'double':
            DimVarType = np.float64
            Fill_Value = np.float64(1.0e20)
        elif CmorType == 'int':
            DimVarType = np.int32
            Fill_Value = np.int32(-9999)
        else:
            DimVarType = np.float32
            Fill_Value = np.float32(1.0e20)

        return (DimVarType, Fill_Value) if ifv else DimVarType


    def UnitsConversion(self):

       # make the units be paresed by udunits
       src_units = []
       for vn in self.syvars:
           vr = self.modvar[vn.name]["attributes"]
           print ('xxxxx', vr, vn.name, self.syvars)

           if "units" not in vr.keys():
               print ("Warning: units are not in the record")
               if "volume mixing ratio" in vr["long_name"]:
                   vr["units"] = "ppm"
               elif 'hybrid B coefficient' in vr["long_name"]:
                   vr["units"] = "1"
               elif 'hybrid A coefficient' in vr["long_name"]:
                   vr["units"] = "1.e-5Pa"
               elif 'gll grid areas' in vr["long_name"]:
                   vr["units"] = "steradian"
               elif 'land fraction' in vr["long_name"]:
                   vr["units"] = "1"
           print (vr, vr["units"])
           if vr["units"] == 'proportion' or vr["units"] == 'unitless' or vr["units"] == 'none' or vr["units"] == 'fraction':
              vr["units"] = '1'
           elif vr["units"] == 'mm/s':
              vr["units"] = 'kg m-2 s-1'
           elif vr["units"] == 'm/s' and ("precipitation" in vr["long_name"] or "snow" in vr["long_name"]):
              vr["units"] = '1000kg m-2 s-1'

           if vn.name == 'TWS':
              self.modvar[vn.name]["attributes"]["units"] ="kg m-2"

           try:
              src_units.append(cf.Unit(vr["units"]))
           except:
              if 'C' in vr["units"]:
                 vr["units"] = vr["units"].replace('C','')
                 src_units.append(cf.Unit(vr["units"]))
              elif 'N' in cmordict['_children'][0]['units']:
                 vr["units"] = vr["units"].replace('N','')
                 src_units.append(cf.Unit(vr["units"]))
              else:
                 print ("cannot convert")

       print ("relationship", self.cmrela)

       if ('+' in self.cmrela or '-' in self.cmrela):
           subeq = re.split(r'\+|\-', self.cmrela)
           sunit = []
           print ('xumdeb', subeq, self.cmrela)
           for eq in subeq:

               if not eq.isnumeric() and eq != '':
                  symeq = sympy.sympify(eq)
                  symvr = tuple(symeq.free_symbols)
                  symfn = sympy.lambdify(symvr, symeq)

                  symut = []
                  for sv in symvr:
                      symut.append(cf.Unit(self.modvar[sv.name]["attributes"]["units"]))

                  print (eq, 'xxx', symut)
                  sunit.append(symfn(*symut))
           print ('yyy', sunit)
           if len(set(sunit)) == 1:
               rlt_units = sunit[0]
           else:
               print (sunit)
               raise ValueError
       else:
           rlt_units = self.syfunc(*src_units)


       if vr["units"] == "steradian":
           rlt_units = cf.Unit("m2")   # conversion is done by relationship
       #check if the  
       tag_units = cf.Unit(self.cmrvar["attributes"]["units"])

       print (tag_units, '---', rlt_units)
       if tag_units.is_convertible(rlt_units):
          #if tag_units.convert(1., rlt_units) != 1.:
          #   cmortask.append("units")
          #print(tag_units.convert(1., rlt_units))
          self.src_units = rlt_units
          self.tag_units = tag_units
       else:
          print("units are not convertible\n")
          sys.exit()

    def GetAttributes(self, OptGblAttrs={}, OptVarAttrs={}):
       #Gbl and Var attrbutes
       import load_CMIP6CV as cv
       import uuid

       mipcv = cv.cmip6cv()
       # from dreq 
       DreqInput = {}
       DreqInput["frequency"            ] = self.cmrvar["attributes"]["frequency"]

       #"land landIce" for SOILICE
       #-if len(self.cmrvar["attributes"]["modeling_realm"].split(' ')) > 1:
       #-    if self.cmrvar["attributes"]["miptable"][0:2] == 'Lm':
       #-        DreqInput["realm"                ] = "land"
       #-    elif self.cmrvar["attributes"]["miptable"][0:2] == 'Li':
       #-        DreqInput["realm"                ] = "landIce"
       #-    elif self.cmrvar["attributes"]["miptable"][0:2] == 'Am':
       #-        DreqInput["realm"                ] = "atmos"
       #-    else:
       #-        print (self.cmrvar["attributes"]["modeling_realm"])
       #-        raise ValueError
       #-else:
       #-     DreqInput["realm"                ] = self.cmrvar["attributes"]["modeling_realm"]
       DreqInput["realm"                ] = self.cmrvar["attributes"]["modeling_realm"]
       DreqInput["table_id"             ] = self.cmrvar["attributes"]["miptable"]
       DreqInput["variable_id"          ] = self.cmrvar["cmvar"]

       mipcv.GetGlobalAttributes(self.userinput, DreqInput)
       #print (mipcv.ReqGblAttrs)
       self.GblAttrs = {**mipcv.ReqGblAttrs, **OptGblAttrs}


       # valattribues

       # from data requestion python API
       ReqVarAttrs = {}
       ReqVarAttrs.update({'description'   :  self.cmrvar["attributes"]["description"]})
       ReqVarAttrs.update({'comment'       :  self.cmrvar["attributes"]["comment"]})
       ReqVarAttrs.update({'standard_name' :  self.cmrvar["attributes"]["standard_name"]})
       ReqVarAttrs.update({'long_name'     :  self.cmrvar["attributes"]["longname"]})
       ReqVarAttrs.update({'title'         :  self.cmrvar["attributes"]["longname"]})
       ReqVarAttrs.update({'id'            :  self.cmrvar["attributes"]["id"]})
       ReqVarAttrs.update({'variable_id'   :  self.cmrvar["attributes"]["variable_id"]})
       ReqVarAttrs.update({'frequency'     :  self.cmrvar["attributes"]["frequency"]})
       ReqVarAttrs.update({'positive'      :  self.cmrvar["attributes"]["positive"]})
       ReqVarAttrs.update({'modeling_realm':  self.cmrvar["attributes"]["modeling_realm"]})
       ReqVarAttrs.update({'type'          :  self.cmrvar["attributes"]["type"]})
       ReqVarAttrs.update({'prov'          :  self.cmrvar["attributes"]["prov"]})
       ReqVarAttrs.update({'mipTable'      :  self.cmrvar["attributes"]["miptable"]})
       ReqVarAttrs.update({'positive'      :  self.cmrvar["attributes"]["positive"]})
       ReqVarAttrs.update({'units'         :  self.cmrvar["attributes"]["units"]})


       for ak in self.cmrvar["attributes"].keys():
           if ak == 'miptable':
              ReqVarAttrs.update({'mipTable'      :  self.cmrvar["attributes"]["miptable"]})
           elif ak != "_FillValue" and ak != "missing_value" and ak != "priority" and ak != "mip":
              ReqVarAttrs.update({ak              :  self.cmrvar["attributes"][ak]})

       print ('attrbute length', set(self.cmrvar["attributes"])-set(ReqVarAttrs.keys()))
       print ('attrbute length', set(ReqVarAttrs.keys()) - set(self.cmrvar["attributes"]))


       #ReqVarAttrs.update({'_FillValue'    :  'inline'})
       #ReqVarAttrs.update({'missing_value' :  self.userinput["missing_value"]})


       self.VarAttrs = {**ReqVarAttrs, **OptVarAttrs}
       print (ReqVarAttrs)
   

       #if ncdf_var.dtype == np.single:
       #    ReqVarAttrs.setncattr('_FillValue'    , np.float32(1.e20))
       #    ReqVarAttrs.setncattr('missing_value' , np.float32(1.e20))
       #elif ncdf_var.dtype == np.double:
       #    ncdf_var.setncattr('_FillValue'    , np.float64(1.e20))
       #    ncdf_var.setncattr('missing_value' , np.float64(1.e20))
       #elif ncdf_var.dtype == np.int16:
       #    ncdf_var.setncattr('_FillValue'    , np.int16(-9999))
       #    ncdf_var.setncattr('missing_value' , np.int16(-9999))


       # Now we can get the filename and directory name

       #<variable_id>_<table_id>_<experiment_id >_<source_id>_<member_id>_<grid_label>[_<time_range>].nc

       if ReqVarAttrs["mipTable"] == "fx":
          self.cmfnam = "{}_".format(ReqVarAttrs["variable_id"]) +\
                        "{}_".format(ReqVarAttrs["mipTable"]) +\
                        "{}_".format(mipcv.ReqGblAttrs["source_id"]) +\
                        "{}_".format(mipcv.ReqGblAttrs["experiment_id"]) +\
                        "{}_".format(mipcv.ReqGblAttrs["variant_label"]) +\
                        "{}.nc".format(mipcv.ReqGblAttrs["grid_label"]) 

       else:
          self.cmfnam = "{}_".format(ReqVarAttrs["variable_id"]) +\
                        "{}_".format(ReqVarAttrs["mipTable"]) +\
                        "{}_".format(mipcv.ReqGblAttrs["source_id"]) +\
                        "{}_".format(mipcv.ReqGblAttrs["experiment_id"]) +\
                        "{}_".format(mipcv.ReqGblAttrs["variant_label"]) +\
                        "{}_".format(mipcv.ReqGblAttrs["grid_label"]) +\
                        "{}.nc".format(self.userinput["time_range"])


       #Directory structure = <mip_era>/
       #                        <activity_id>/
       #                          <institution_id>/
       #                            <source_id>/
       #                              <experiment_id>/
       #                                <member_id>/         <== variant_label
       #                                  <table_id>/
       #                                    <variable_id>/
       #                                      <grid_label>/
       #                                        <version>/

       #"{}/".format(mipcv.ReqGblAttrs["activity_id"]) +\
       self.cmdnam = "{}/".format("CMIP6") +\
                     "{}/".format(self.userinput["mip_name"]) +\
                     "{}/".format(mipcv.ReqGblAttrs["institution_id"]) +\
                     "{}/".format(mipcv.ReqGblAttrs["source_id"]) +\
                     "{}/".format(mipcv.ReqGblAttrs["experiment_id"]) +\
                     "{}/".format(mipcv.ReqGblAttrs["variant_label"]) +\
                     "{}/".format(mipcv.ReqGblAttrs["table_id"]) +\
                     "{}/".format(mipcv.ReqGblAttrs["variable_id"]) +\
                     "{}/".format(mipcv.ReqGblAttrs["grid_label"]) +\
                     "{}/".format("v20200611")

       return


    def CmorvarCompute(self):

       #check the consistency with dimensions
       ReqDimLong = [self.cmrvar['dimensions'][dk]["out_name"] for dk in self.cmrvar['dimensions'].keys() if dk != 'shape']
       ReqDimShrt = [self.cmrvar['dimensions'][dk]["title"] for dk in self.cmrvar['dimensions'].keys() if dk != 'shape']

       if "E3SM" in self.userinput["source_id"] or "CESM" in self.userinput["source_id"]:
           ReqDimLong.append("levgrnd")
           ReqDimShrt.append("levgrnd")

           ReqDimLong.append("ilev")
           ReqDimShrt.append("ilev")

       # operation is dependent, the operator and operating dimension is set
       if self.cmdims == {}:
          xlist = []
          self.vardimnm = ()
          self.vardimno = ()

          print ('xxx1', self.cmdmns)
          if self.cmdmns != '{}':
              if 'shape' in self.cmdmns:
                  p = re.compile(r"shape:(.*)}")
                  commondimnms = p.findall(self.cmdmns)[0].split(',')
                  self.vardimnm = commondimnms

              else:
                  commondimnms = []
                  for sv in self.syvars:
                      fnabs = self.__FindModVarFile(sv)
                      #dimension conform
                      if os.path.exists(fnabs):
                          with nc4.Dataset(fnabs, "r") as f:
                              if len(commondimnms) < len(f.variables[sv.name].dimensions):
                                 commondimnms =  list(f.variables[sv.name].dimensions)
                                 self.vardimnm = f.variables[sv.name].dimensions
                                 self.vardimno = f.variables[sv.name].shape

                      else:
                          print ("Add {} into the variable list".format(sv.name))
                          print (fnabs, self.cmdmns)
                          self.missingvars.append(sv.name)

                  print (commondimnms)

              actshape = [1]*len(commondimnms)
              for sv in self.syvars:
                  fnabs = self.__FindModVarFile(sv)
                  if os.path.exists(fnabs):
                      with nc4.Dataset(fnabs, "r") as f:
                           if len(f.variables[sv.name].dimensions) >=4:
                              print ("huge array skip")
                           if f.variables[sv.name].dimensions == tuple(commondimnms):
                               xlist.append(np.float32(f.variables[sv.name][...]))
                               self.vardimnm = commondimnms
                               self.vardimno = f.variables[sv.name][...].shape
                           else:
                               conshape = [1]*len(commondimnms)
                               for dm, ds in zip(f.variables[sv.name].dimensions, f.variables[sv.name].shape):
                                   if dm in tuple(commondimnms):
                                       k = commondimnms.index(dm)
                                       conshape[k] = ds
                                       actshape[k] = ds
                               xlist.append(np.float32(f.variables[sv.name][...].reshape(tuple(conshape))))
              self.vardimnm = commondimnms
              self.vardimno = tuple(actshape)

              print ('actshape', actshape)
                  
          else:   # no speicifc law

              print ("no special law")
              for sv in self.syvars:
                  fnabs = self.__FindModVarFile(sv)
                  print('fnabs=', fnabs, os.path.exists(fnabs))
                  if os.path.exists(fnabs):
                      print ('aaaa')
                      with nc4.Dataset(fnabs, "r") as f:
                           print('xumdeb1', f.variables[sv.name].dimensions)
                           if len(f.variables[sv.name].dimensions) >=4:
                              print ("huge array skip")
                           xlist.append(f.variables[sv.name][...])
                           self.vardimnm = f.variables[sv.name].dimensions
                           self.vardimno = f.variables[sv.name].shape

          for dn in self.vardimnm:
              print ('dn', dn)
              if dn == 'plev':
                  dn = 'plev19'
              if dn not in ReqDimLong and dn not in ReqDimShrt:
                  print ("error in dimesion")
                  sys.exit()

          if len(xlist) == len(self.syvars):
             if '@' in self.cmdmns and self.cmdmns.split('@')[1].split(';')[0] !='':
                 print (len(self.cmdmns.split('@')))
                 print (self.cmdmns.split('@'))
                 operator = self.cmdmns.split('@')[0]



                 temp = self.cmdmns.split('@')[1].split(';')[0].split(':')
                 opdimnam = temp[0]
                 opdimval = temp[1]

                 if '-' not in opdimval:
                    kbgn = int(opdimval)
                    kend = int(opdimval)
                 else:
                    kbgn = int(opdimval.split('-')[0])
                    kend = int(opdimval.split('-')[1])
                 # remove the operated dimensions
                 rmindx = self.vardimnm.index(opdimnam)
                 self.vardimnm = tuple([x for x in self.vardimnm if x != opdimnam])
                 self.vardimno = self.vardimno[:rmindx] + self.vardimno[rmindx+1:]

                 print (opdimnam, kbgn, kend)
                 self.cmvval = np.sum(self.syfunc(*xlist), axis=1)
             else:

                 for x in xlist:
                     print ('dtype', x.dtype)
                 self.cmvval = self.syfunc(*xlist)
                 print (self.cmvval.dtype)

             print ('xumdeb', self.cmrela)
             print ('shape', self.cmvval.shape)


       else:
          # the variables are independent and the operator is sum in default:
          xlist = []
          self.vardimnm = ()
          self.vardimno = ()
          print ("in complex computing")
          for sv in self.syvars:
              fnabs = "{}/{}/{}_{}.nc".format(self.userinput["ModelRltDir"], self.cmrvar["attributes"]["miptable"], 
                      sv.name, self.userinput["time_range"].replace('-','_'))
              if os.path.exists(fnabs):
                  with nc4.Dataset(fnabs, "r") as f:
                      if len(f.variables[sv.name].dimensions) >=4:
                         print ("huge array skip")
                         print (sv.name, self.cmdims)
                      if sv.name in self.cmdims.keys():
                         if ':' not in self.cmdims[sv.name]:
                            xlist.append(f.variables[sv.name][:,int(self.cmdims[sv.name]),:,:])
                         else:
                            kbgn = int(self.cmdims[sv.name].split(':')[0])
                            kend = int(self.cmdims[sv.name].split(':')[1])

                            # assume it is sum, but we can define 
                            varsum = np.sum(f.variables[sv.name][:,kbgn:kend+1,:,:], axis=1)
                            xlist.append(varsum)
                         #remove the dimensions
                         vardimnm = (f.variables[sv.name].dimensions[0],)
                         for vn in f.variables[sv.name].dimensions[2:]:
                             vardimnm = vardimnm + (vn,)
                         vardimno = (f.variables[sv.name].shape[0],)
                         for vn in f.variables[sv.name].shape[2:]:
                             vardimno = vardimno + (vn,)
                         print(vardimnm, vardimno)
                         if len(vardimnm) > len(self.vardimnm):
                               self.vardimnm = vardimnm
                               self.vardimno = vardimno
                      else:
                         xlist.append(f.variables[sv.name][...])
                         vardimnm = f.variables[sv.name].dimensions
                         vardimno = f.variables[sv.name].shape
                         if len(vardimnm) > len(self.vardimnm):
                               self.vardimnm = vardimnm
                               self.vardimno = vardimno
              else:
                  print ("Add {} into the variable list".format(sv.name))
                  print (fnabs, 'else')
                  self.missingvars.append(sv.name)

          for dn in self.vardimnm:
              if dn not in ReqDimLong and dn not in RegDimShrt:
                  print ("error in dimesion")
                  sys.exit()

          if len(xlist) == len(self.syvars):
              self.cmvval = self.syfunc(*xlist)
              print ('shape', self.cmvval.shape)

          return

    def Cmorization(self):

        Cmip6CompressionOpt = {"zlib": True, "complevel":1, "shuffle":False}
        Cmip6NetcdfOpt = {"format": "NETCDF4_CLASSIC"}

        if self.cmvval is None:
            return

        for sv in self.syvars:
            if sv.name+':inline' not in self.cmdmns:
                break
        if 'plev19' in self.cmrvar["dimensions"]:   # standard pressure level
             fname = sv.name + '_plev19'
        elif self.cmrvar["cmvar"] == 'pfull' or self.cmrvar["cmvar"] == 'phalf':
             fname = "tmp_vert"
        elif 'fx' in self.cmrvar["attributes"]["miptable"]:
             fname = 'atm_area_landfrac'
        else:
             fname = sv.name
        fnabs = "{}/{}/{}_{}.nc".format(self.userinput["ModelRltDir"], self.cmrvar["attributes"]["miptable"], 
                 fname, self.userinput["time_range"].replace('-','_'))
        origvn = sv.name

        if not os.path.isdir(self.userinput["CmorRltDir"] + "/" + self.cmdnam):
            os.system("mkdir -p " + self.userinput["CmorRltDir"] + "/" + self.cmdnam)
        fncmr = self.userinput["CmorRltDir"] + "/" + self.cmdnam + "/" + self.cmfnam

        bnd_name = None

        fnbnd = "{}/{}/{}".format(self.userinput["ModelRltDir"], self.cmrvar["attributes"]["miptable"],"bounds.nc")

        with nc4.Dataset(fnabs, "r") as forig, nc4.Dataset(fncmr, "w", **Cmip6NetcdfOpt) as fcmor, nc4.Dataset(fnbnd, "r") as fbnds:
            #create dimensions
            for name, dimension in forig.dimensions.items():
                print ('nam=', name, forig.dimensions, self.vardimnm)
                if name in ['nb', 'nbnd']:
                    bnd_name = name
                #if name in self.vardimnm or name in ['nb', 'nbnd', 'hist_interval']:
                if name in self.vardimnm or name in ['nb', 'nbnd']:
                    if name == "levgrnd":
                        dim_name = "depth"
                    elif name == "ilev":
                        dim_name = "lev"
                    else:
                        dim_name = name
                    fcmor.createDimension(dim_name, (len(dimension) if not dimension.isunlimited() else None))
            if not bnd_name:
               bnd_name = "nbnd"  #default
               fcmor.createDimension(bnd_name, 2)

            # copy all file data except for the excluded
            # extended dimensions
            DimBnd = {}
            DimKey = {}
            for name, variable in forig.variables.items():
                if name != origvn:
                   extdim=True
                   for vd in variable.dimensions:
                        if vd not in self.vardimnm:
                           extdim = False
                           break
                   if "bnd" in name or "bound" in name or "hist_interval" in name:
                       extdim = True
                   if extdim:
                      print (name, variable.datatype,  variable.dimensions)
                      cmname = None
                      dimkey = None
                      dmnlst = None

                      print("name=", name)
                      #check the dim variable name is in the cmor standard
                      for dk in self.cmrvar['dimensions'].keys():
                          print ('xxx', dk, self.cmrvar['dimensions'][dk])
                          if name in dk:
                             cmname = self.cmrvar['dimensions'][dk]["out_name"]
                             dmnlst = variable.dimensions
                             dimkey = dk
                             break

                      if name == "levgrnd":
                          cmname = 'depth'
                          print ('indepth', name, variable.datatype, variable.dimensions)
                          dmnlst = ("depth")
                          dimkey = 'sdepth'

                      if name == "ilev":
                          cmname = self.cmrvar['dimensions']["alevhalf"]["out_name"]
                          dmnlst = ("lev")
                          dimkey = "alevhalf"

                      if  name == "lev":
                          cmname = self.cmrvar['dimensions']["alevel"]["out_name"]
                          dmnlst = ("lev")
                          dimkey = "alevel"

                      print ('dimkey', dimkey)

                      if dimkey:
                          if "type" in self.cmrvar['dimensions'][dimkey].keys():
                              DimVarType = self.__Map2npType(self.cmrvar['dimensions'][dimkey]['type'], False)
                              #if self.cmrvar['dimensions'][dimkey]['type'] == 'real':
                              #    DimVarType = np.float32
                              #elif self.cmrvar['dimensions'][dimkey]['type'] == 'double':
                              #    DimVarType = np.float64
                              #elif self.cmrvar['dimensions'][dimkey]['type'] == 'int':
                              #    DimVarType = np.int32
                          else:
                              DimVarType = np.float32

                          x = fcmor.createVariable(cmname, DimVarType, dmnlst)

                      print ('xxx2', name, dimkey)
                      if cmname:  
                         #fcmor[name].setncatts(forig[name].__dict__)

                         if 'time' in cmname:

                             print (forig[name].units)
                             oldtimeunit = cf.Unit(forig[name].units, calendar=forig[name].calendar)
                             newtimeunit = cf.Unit(self.cmrvar['dimensions']['time']["units"], forig[name].calendar)

                             fcmor[cmname][:] = oldtimeunit.convert(forig[name][:], newtimeunit)

                         else:

                             print (cmname, name)
                             fcmor[cmname][:] = forig[name][:]

                         # fix errors caused by JS
                         for k, v in self.cmrvar['dimensions'][dimkey].items():
                             if not isinstance(v, str):
                                 try:
                                     DimVarType = self.__Map2npType(self.cmrvar['dimensions'][dimkey]["type"], False)

                                     if isinstance(v, np.ndarray):
                                        self.cmrvar['dimensions'][dimkey][k] = v.astype(DimVarType)
                                     else:

                                        if isinstance(v, list):
                                           self.cmrvar['dimensions'][dimkey][k] = np.array([v]).astype(DimVarType) 
                                        else:
                                           self.cmrvar['dimensions'][dimkey][k] = np.array([v]).astype(DimVarType)[0]

                                     print (np.array([v]).astype(DimVarType)[0], v)
                                     
                                     #if self.cmrvar['dimensions'][dimkey]["type"] == "real":
                                     #   self.cmrvar['dimensions'][dimkey][k] = np.float32(v)
                                     #elif self.cmrvar['dimensions'][dimkey]["type"] == "double":
                                     #   self.cmrvar['dimensions'][dimkey][k] = np.float64(v)
                                     #elif self.cmrvar['dimensions'][dimkey]["type"] == "int":
                                     #   self.cmrvar['dimensions'][dimkey][k] = np.int32(v)
                                 except:
                                     raise ValueError
                         #+mxu in order to pass cfchecker, need to remove type:double attributes
                         # otherwise there will be warning, 
                         #fcmor[cmname].setncatts(self.cmrvar['dimensions'][dimkey])
                         fcmor[cmname].setncatts({xmkey:xmval for xmkey, xmval in self.cmrvar['dimensions'][dimkey].items() if xmkey != "title"})


                         print ('xxx3', dimkey, self.cmrvar['dimensions'][dimkey].keys())
                         if "bounds" in self.cmrvar['dimensions'][dimkey].keys():
                             DimBnd[name] = self.cmrvar['dimensions'][dimkey]["bounds"]
                             DimKey[name] = dimkey
                             print ('xxx4', dimkey, name, DimBnd[name])

            # write the bnd data
            for name, variable in forig.variables.items():
                if 'nb' in variable.dimensions or 'nbnd' in variable.dimensions or 'nbound' in variable.dimensions or 'hist_interval' in variable.dimensions:
                    dimvnm = name.split('_')[0]
                    if dimvnm in DimBnd.keys():
                       #x = fcmor.createVariable(DimBnd[dimvnm], variable.datatype, variable.dimensions)
                       if "type" in self.cmrvar["dimensions"][DimKey[dimvnm]].keys():

                           DimVarType = self.__Map2npType(self.cmrvar["dimensions"][DimKey[dimvnm]]["type"], False)
                       else:
                           DimVarType = np.float32

                       DimBndDims = tuple([bnd_name if x == 'hist_interval' else x for x in variable.dimensions])
                       print ('xumdeb', variable.dimensions, DimBndDims)
                       #x = fcmor.createVariable(DimBnd[dimvnm], DimVarType, variable.dimensions)
                       x = fcmor.createVariable(DimBnd[dimvnm], DimVarType, DimBndDims, **Cmip6CompressionOpt)

                       fcmor[DimBnd[dimvnm]][...] = forig[name][...].astype(DimVarType)

                       #+mxu there should be no attrributes for bnd variables based on CF-1.7
                       #-fcmor[DimBnd[dimvnm]].setncatts({"units": self.cmrvar['dimensions'][DimKey[dimvnm]]["units"]})
                       if 'time' in DimBnd[dimvnm]:
                           fcmor[DimBnd[dimvnm]][...] = oldtimeunit.convert(forig[name][...], newtimeunit)
                           #-mxu
                           #-fcmor[DimBnd[dimvnm]].setncatts({"calendar": self.cmrvar['dimensions'][DimKey[dimvnm]]["calendar"]})
                           #make sure the time is in the middle of bnds
                           fcmor["time"][:] = 0.5 * (fcmor[DimBnd[dimvnm]][:,0] + fcmor[DimBnd[dimvnm]][:,1])
                # special handler
                if name == "levgrnd" and name in DimBnd.keys():

                   dimvnm = name
                   if "type" in self.cmrvar["dimensions"][DimKey[dimvnm]].keys():

                       DimVarType = self.__Map2npType(self.cmrvar["dimensions"][DimKey[dimvnm]]["type"], False)
                   else:
                       DimVarType = np.float32


                   x = fcmor.createVariable(DimBnd["levgrnd"], DimVarType, ("depth", bnd_name), **Cmip6CompressionOpt)

                   #-mxu
                   #-fcmor[DimBnd["levgrnd"]].setncatts({"units": self.cmrvar['dimensions'][DimKey["levgrnd"]]["units"]})

                   fcmor[DimBnd["levgrnd"]][0, 0] = 0.0  
                   fcmor[DimBnd["levgrnd"]][-1,1] = 42.1023  

                   for k, d in enumerate(fcmor["depth"][:]): 
                       if k < fcmor["depth"][:].shape[0] - 1:
                          fcmor[DimBnd["levgrnd"]][k+1,0] = 0.5 * (d + fcmor["depth"][k+1])  
                          fcmor[DimBnd["levgrnd"]][k  ,1] = 0.5 * (d + fcmor["depth"][k+1])
                   
                if (name == "lev" or name =="ilev") and name in DimBnd.keys():

                   dimvnm = name
                   if "type" in self.cmrvar["dimensions"][DimKey[dimvnm]].keys():
                       DimVarType = self.__Map2npType(self.cmrvar["dimensions"][DimKey[dimvnm]]["type"], False)
                   else:
                       DimVarType = np.float32

                   x = fcmor.createVariable(DimBnd[name], DimVarType, ("lev", bnd_name), **Cmip6CompressionOpt)

                   #-mxu
                   #-fcmor[DimBnd[name]].setncatts({"units": self.cmrvar['dimensions'][DimKey[name]]["units"]})

                   arrsize = forig["ilev"][...].shape[0]
                   if name == "ilev":
                      fcmor[DimBnd[name]][0:-1,0] = forig["ilev"][0:-1]
                      fcmor[DimBnd[name]][0:-1,1] = forig["ilev"][1:arrsize]
                      fcmor[DimBnd[name]][-1,0] = fcmor[DimBnd[name]][-2,0]
                      fcmor[DimBnd[name]][-1,1] = fcmor[DimBnd[name]][-2,1]
                   else:
                      fcmor[DimBnd[name]][:,0] = forig["ilev"][0:-1]
                      fcmor[DimBnd[name]][:,1] = forig["ilev"][1:arrsize]

                      # copy variable attributes all at once via dictionary

            if 'nb' not in forig.dimensions and 'nbound' not in forig.dimensions and 'nbnd' not in forig.dimensions:
                for dimvnm in DimBnd.keys():
                   if DimBnd[dimvnm] in fbnds.variables.keys():
                      if "type" in self.cmrvar["dimensions"][DimKey[dimvnm]].keys():
                          DimVarType = self.__Map2npType(self.cmrvar["dimensions"][DimKey[dimvnm]]["type"], False)
                      else:
                          DimVarType = np.float32
                      x = fcmor.createVariable(DimBnd[dimvnm], DimVarType, fbnds.variables[DimBnd[dimvnm]].dimensions, **Cmip6CompressionOpt)

                      fcmor[DimBnd[dimvnm]][...] = fbnds[DimBnd[dimvnm]][...].astype(DimVarType)
                      #-
                      #-fcmor[DimBnd[dimvnm]].setncatts({"units": self.cmrvar['dimensions'][DimKey[dimvnm]]["units"]})

                
            if self.userinput['timeshift'] != 0:
                timeshift = self.userinput['timeshift']
                #time_bounds=time_bounds+$(($difyear*(-365))); time=time_bounds(:,1);
                print ("TIME SHIFTING !!!!!!")
                if 'time_bounds' in fcmor.variables.keys():
                   fcmor['time_bounds'][...] = fcmor['time_bounds'][...] - timeshift * 365
                   fcmor['time'][:] = fcmor['time_bounds'][:,1]
                elif 'time_bnds' in fcmor.variables.keys():

                   print (timeshift)
                   print (fcmor['time_bnds'][0,0])
                   fcmor['time_bnds'][...] = fcmor['time_bnds'][...] - timeshift * 365
                   fcmor['time'][:] = fcmor['time_bnds'][:,1]
                   print (fcmor['time_bnds'][0,0])
                else:
                   print ('cannot time shifting due to no time variables')

            VarDataType, fillvalue = self.__Map2npType(self.cmrvar["attributes"]["type"], True)
            self.userinput["missing_value"] = fillvalue

            var_dimnms = ['depth' if dmn == 'levgrnd' else dmn for dmn in self.vardimnm]
            var_dimnms = ['lev' if dmn == 'ilev' else dmn for dmn in var_dimnms]

            print (var_dimnms, self.vardimnm, fcmor.dimensions)
            #fcmor.createVariable(self.cmrvar["cmvar"], forig.variables[origvn].datatype, var_dimnms, complevel=2, shuffle=True, fill_value=fillvalue)
            #-fcmor.createVariable(self.cmrvar["cmvar"], VarDataType, var_dimnms, complevel=2, shuffle=True, fill_value=fillvalue)
            fcmor.createVariable(self.cmrvar["cmvar"], VarDataType, var_dimnms, fill_value=fillvalue, **Cmip6CompressionOpt)

            print('xumdeb', self.cmrvar["cmvar"], forig.variables[origvn].datatype)

            if self.cmrvar["cmvar"] == "burntFractionAll":
                print('xudeb', self.src_units, self.tag_units, self.cmrvar["cmvar"])

                ndays=[0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
                for i, ti in enumerate(fcmor['time'][:]):
                    #print (fcmor['time'].getncattr('units'), fcmor['time'].getncattr('calendar')) 
                    tempdate=nc4.num2date(ti, units=fcmor['time'].getncattr('units'), calendar=fcmor['time'].getncattr('calendar'))

                    #print (tempdate.month, ndays[tempdate.month])

                    fcmor[self.cmrvar["cmvar"]][i,:,:] = self.src_units.convert(self.cmvval[i,:,:], self.tag_units) * 86400 * ndays[tempdate.month]



            else:
                fcmor[self.cmrvar["cmvar"]][...] = self.src_units.convert(self.cmvval[...], self.tag_units)

            self.VarAttrs["missing_value"] = fillvalue


            # variable attributes
            external_variables_val = "none"
            for vk in self.VarAttrs.keys():

                if vk == "cell_measures" and self.VarAttrs[vk] != "none" and self.VarAttrs[vk]:
                   print (vk, self.VarAttrs[vk], 'xxx')
                   external_variables_val = self.VarAttrs[vk].split(':')[1].strip()

                # ignore empty positive and title
                if (vk == "positive" and not self.VarAttrs[vk]) or vk == "title":
                    print ("skip " + vk)
                elif (vk == "cell_measures" and not self.VarAttrs[vk]):
                    print ("skip " + vk)
                else:
                    fcmor[self.cmrvar["cmvar"]].setncattr(vk, self.VarAttrs[vk])

            # add CF attributes
            fcmor[self.cmrvar["cmvar"]].setncattr("coordinates", " ".join(var_dimnms))
            # write the global attributes

            if external_variables_val != "none":
               setattr(fcmor, 'external_variables',external_variables_val)
            for gk in self.GblAttrs.keys():
                setattr(fcmor, gk, self.GblAttrs[gk])
                if gk == "realization_index":

                    #setattr(fcmor, gk, np.int32(self.GblAttrs[gk]))
                    print (type(fcmor.getncattr(gk)))
                    print (type(self.GblAttrs[gk]))
                if gk == "source" and self.userinput["mod_name"] == "E3SM":
                    e3sm_sourceline="E3SM-1-1 (2019): aerosol: MAM4 (MAM4, with resuspension, marine organics, and secondary organics same grid as atmos); "\
                                    "atmos: EAM (EAMv1.1, cubed sphere spectral-element grid; 5400 elements with p=3; 1 deg average grid spacing; "\
                                    "90 x 90 x 6 longitude/latitude/cubeface; 72 levels; top level 0.1 hPa);"\
                                    "atmosChem: LINOZ (LINOZ v2, Troposphere specified oxidants for aerosols, Stratosphere linearized interactive ozone, same grid as atmos); "\
                                    "land: ELM (ELMv1.1, same grid as atmos; active biogeochemistry using the Converging Trophic Cascade plant and soil carbon "\
                                    "and nutrient mechanisms to represent carbon, nitrogen and phosphorus cycles; "\
                                    "MOSART v1.1, 0.5 degree latitude/longitude grid); "\
                                    "ocean: MPAS-Ocean (MPAS-Oceanv6.0, oEC60to30 unstructured SVTs mesh with 235160 cells and 714274 edges, "\
                                    "variable resolution 60 km to 30 km; 60 levels; top grid cell 0-10 m); "\
                                    "ocnBgchem: BEC (BEC, Biogeochemical Elemental Cycling model, NPZD-type with C/N/P/Fe/Si/O; same grid as ocean); "\
                                    "seaIce: MPAS-Seaice (MPAS-Seaicev6.0, same grid as ocean);"
                    setattr(fcmor, gk, e3sm_sourceline)



# test code

if __name__ == "__main__":

   import subprocess
   import datetime
   import uuid

   mod_name = 'E3SM'
   exp_name = '1pctCO2-bgc'
   mip_name = 'C4MIP'
   ins_name = 'RUBISCO'
   tab_name = 'Lmon'
   rel_name = 'land'
   
   datadir="/global/cscratch1/sd/minxu/data/20191123.CO21PCTRAD_RUBISCO_CNPCTC20TR_OIBGC.I1900.ne30_oECv3.compy/rgr/atm/"
   cmordir="/global/homes/m/minxu/scratch/data/mytest/"
   
   timeshift=1899

   # areacella_fx_CESM2_esm-hist_r1i1p1f1_gn.nc

   UserInput = {}
   #UserInput["mip_name"             ] = "CMIP"
   #UserInput["tab_name"             ] = "fx"
   #UserInput["rel_name"             ] = "atmos"
   UserInput["mip_name"             ] = "C4MIP"
   UserInput["tab_name"             ] = "Lmon"
   UserInput["rel_name"             ] = "land"
   #UserInput["tab_name"             ] = "Amon"
   #UserInput["rel_name"             ] = "atmos"
   UserInput["mod_name"             ] = "E3SM"
   UserInput["time_range"           ] = "000101-015012"
   UserInput["ModelRltDir"          ] = "/global/cscratch1/sd/minxu/data/20191123.CO21PCTBGC_RUBISCO_CNPCTC20TR_OIBGC.I1900.ne30_oECv3.compy/rgr/"
   UserInput["CmorRltDir"           ] = "/global/cfs/cdirs/m2467/prj_minxu/cmip6/"

   UserInput["experiment_id"        ] = '1pctCO2-bgc'
   UserInput["institution_id"       ] = 'RUBISCO'
   UserInput["grid_label"           ] = "gr"
   UserInput["nominal_resolution"   ] = "100 km"
   UserInput["further_info_url"     ] = "https://bgc-feedbacks.org"
   UserInput["creation_date"        ] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
   UserInput["source_component"     ] = ["atmos", "atmosChem", "land", "ocean", "ocnBgchem", "seaIce"]
   UserInput["source_id"            ] = "E3SM-1-1"
   UserInput["forcing_index"        ] = np.int32(1)
   UserInput["initialization_index" ] = np.int32(1)
   UserInput["physics_index"        ] = np.int32(1)
   UserInput["realization_index"    ] = np.int32(1)
   UserInput["variant_info"         ] = "p = 1, for CTC with CNP"
   UserInput["branch_method"        ] = "50 yrs spinup on DOE Compy from the 969 yrs piControl spinup on NERSC Edison " \
                                        "to account for the machine differences"
   UserInput["branch_time_in_child" ] = 0.0
   UserInput["branch_time_in_parent"] = 371935.0
   UserInput["parent_time_units"    ] = "days since 0001-01-01"
   # user the shell to get
   myuuid = str(uuid.uuid4())
   UserInput["tracking_id"          ] = "hdl:21.14100/" + myuuid
   UserInput["calendar"             ] = "noleap"
   UserInput["timeunits"            ] = "days since 0001-01-01 00:00:00"


   # get the cmorjson version
   jsonversion = subprocess.check_output("cd cmorjson && git describe --always", shell=True)

   with open("cmorjson/jsonfiles/cmor/{}_{}_{}.json".format(UserInput["mip_name"], UserInput["tab_name"], UserInput["mod_name"])) as rdf:
       cmor_json = json.load(rdf)
   cmordicts = [ cm for cm in cmor_json["variables"] if int(cm["confidence"]) >= 90.0 ]
   
   
   # get all variable definitons
   
   for jsnf in glob.glob("cmorjson/jsonfiles/{}/*".format(mod_name.lower())):
       if "mon" in jsnf and UserInput["rel_name"] in jsnf:
           print (jsnf)
           with open (jsnf) as rdf:
                modeljson = json.load(rdf)
           compdicts = modeljson["variables"]
   

   fmiss = open("missingvars.txt", "w")

   
   for i, cmordict in enumerate(cmordicts):

       #if i > 5:
       #    break

       if cmordict["cmvar"] != "gpp":
           continue

       lwc = lwcmor(cmordict, compdicts, UserInput)
       lwc.UnitsConversion()
       lwc.CmorvarCompute()

       if len(lwc.missingvars) > 0:
          fmiss.write(','.join(lwc.missingvars)+",")

       variant_label = "r{}i{}p{}f{}".format(UserInput["realization_index"], UserInput["initialization_index"], UserInput["physics_index"], UserInput["forcing_index"])
       optgblattrs={}
       optgblattrs.update({"further_info_url":"https://furtherinfo.es-doc.org/CMIP6.{}.{}.{}.{}.{}".format(\
                          UserInput["institution_id"],UserInput["mod_name"],UserInput["experiment_id"],'none', variant_label)})      
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
    
