#!/usr/bin/env python

'''
   Parse all cmip6 CV information to be used by lmt-lwcmor
'''

class cmip6cv():
    def __init__(self):
        import json
        with open ('./CMIP6_CVs/CMIP6_activity_id.json', 'r') as fact:
             self.dict_act = json.load(fact)
        with open ('./CMIP6_CVs/CMIP6_experiment_id.json', 'r') as fexp:
             self.dict_exp = json.load(fexp)
        with open ('./CMIP6_CVs/CMIP6_frequency.json', 'r') as ffrq:
             self.dict_frq = json.load(ffrq)
        with open ('./CMIP6_CVs/CMIP6_grid_label.json', 'r') as fgrd:
             self.dict_grd = json.load(fgrd)
        with open ('./CMIP6_CVs/CMIP6_license.json', 'r') as flic:
             self.dict_lic = json.load(flic)
        with open ('./CMIP6_CVs/CMIP6_nominal_resolution.json', 'r') as fnom:
             self.dict_nom = json.load(fnom)
        with open ('./CMIP6_CVs/CMIP6_required_global_attributes.json', 'r') as fgbl:
             self.dict_gbl = json.load(fgbl)
        with open ('./CMIP6_CVs/CMIP6_source_id.json', 'r') as fsid:
             self.dict_sid = json.load(fsid)
        with open ('./CMIP6_CVs/CMIP6_source_type.json', 'r') as fsty:
             self.dict_sty = json.load(fsty)
        with open ('./CMIP6_CVs/CMIP6_sub_experiment_id.json', 'r') as fsub:
             self.dict_sub = json.load(fsub)
        with open ('./CMIP6_CVs/CMIP6_table_id.json', 'r') as ftid:
             self.dict_tid = json.load(ftid)
        with open ('./CMIP6_CVs/CMIP6_institution_id.json', 'r') as fins:
             self.dict_ins = json.load(fins)
        with open ('./CMIP6_CVs/CMIP6_realm.json', 'r') as frel:
             self.dict_rel = json.load(frel)

        # outputs
        self.ReqGblAttrs = {}
        self.OptGblAttrs = {}

    def retrievebymip(self):
        pass

    def retrievebyinstitute(self):
        pass


    def GetGlobalAttributes(self, UserInput, DreqInput):
        import re
        # default user defined values
        #-UserInput["experiment_id"        ] = '1pctCO2-bgc'
        #-UserInput["institution_id"       ] = 'RUBISCO'
        #-UserInput["grid_label"           ] = "gr"
        #-UserInput["nominal_resolution"   ] = "100 km"
        #-UserInput["further_info_url"     ] = "https://bgc-feedbacks.org"
        #-UserInput["creation_date"        ] = "2020-06-11"
        #-UserInput["source_component"     ] = ["atmos", "atmosChem", "land", "ocean", "ocnBgchem", "seaIce"]
        #-UserInput["source_id"            ] = "E3SM-1-1"
        #-UserInput["forcing_index"        ] = 1
        #-UserInput["initialization_index" ] = 1
        #-UserInput["physics_index"        ] = 1
        #-UserInput["realization_index"    ] = 1 
        #-UserInput["variant_info"         ] = "p = 1, for CTC with CNP"
        #-UserInput["branch_method"        ] = "50 yrs spinup on DOE Compy from the 969 yrs piControl spinup on NERSC Edison \
        #-                                      to account for the machine differences"
        #-UserInput["branch_time_in_child" ] = 0.0
        #-UserInput["branch_time_in_parent"] = 371935.0
        #-UserInput["parent_time_units"    ] = "days since 0001-01-01"
        #-UserInput["tracking_id"          ] = "none"

        #-# from dreq 
        #-DreqInput["frequency"] = "mon"
        #-DreqInput["realm"    ] = "atmos"
        #-DreqInput["table_id" ] = "Amon"
        #-DreqInput["variable_id"] = "none"


        members = [attr for attr in dir(self) 
                 if not callable(getattr(self, attr)) and not attr.startswith("__")]

        members.remove('dict_gbl')
        members.remove('ReqGblAttrs')
        members.remove('OptGblAttrs')

        # constants
        self.ReqGblAttrs["Conventions"         ] = "CF-1.7 CMIP-6.2"
        self.ReqGblAttrs["data_specs_version"  ] = '01.00.32'
        self.ReqGblAttrs["mip_era"             ] = "CMIP6"
        self.ReqGblAttrs["product"             ] = "model-output"


        self.ReqGblAttrs["forcing_index"       ] = UserInput["forcing_index"       ]
        self.ReqGblAttrs["initialization_index"] = UserInput["initialization_index"]
        self.ReqGblAttrs["physics_index"       ] = UserInput["physics_index"       ]
        self.ReqGblAttrs["realization_index"   ] = UserInput["realization_index"   ]
        self.ReqGblAttrs["variant_label"       ] = "r{}i{}p{}f{}".format(UserInput["realization_index"], 
                                                                         UserInput["initialization_index"], 
                                                                         UserInput["physics_index"], 
                                                                         UserInput["forcing_index"])
        self.ReqGblAttrs["variant_info"         ] = UserInput["variant_info"] + " (Information provided by this attribute "\
                    "may in some cases be flawed.  Users can find more comprehensive and up-to-date documentation "\
                    "via the further_info_url global attribute.)"

        # TBD
        self.ReqGblAttrs["branch_method"        ] = UserInput["branch_method"        ]
        self.ReqGblAttrs["branch_time_in_child" ] = UserInput["branch_time_in_child" ]
        self.ReqGblAttrs["branch_time_in_parent"] = UserInput["branch_time_in_parent"] 
        self.ReqGblAttrs["parent_time_units"    ] = UserInput["parent_time_units"    ]
 
        # online determined when it is running
        self.ReqGblAttrs["creation_date"        ] = UserInput["creation_date"]
        self.ReqGblAttrs["tracking_id"          ] = UserInput["tracking_id"]
        self.ReqGblAttrs["variable_id"          ] = DreqInput["variable_id"]

        for rtem in DreqInput['realm'].split(' '):
            if rtem.strip() not in self.dict_rel["realm"]:
               raise ValueError
        self.ReqGblAttrs["realm"] = DreqInput["realm"]   # from drq

        if DreqInput["table_id"] in self.dict_tid["table_id"]:
            self.ReqGblAttrs["table_id" ] = DreqInput["table_id"]    # from drq
        else:
            raise ValueError

        if DreqInput["frequency"] in self.dict_frq["frequency"].keys():
            self.ReqGblAttrs["frequency"] = DreqInput["frequency"]
        else:
            raise ValueError

        #derved from experiment_id
        ExpId=UserInput["experiment_id"]
        self.ReqGblAttrs["activity_id"         ] = " ".join(self.dict_exp["experiment_id"][ExpId]["activity_id"])
        self.ReqGblAttrs["source_type"         ] = " ".join(self.dict_exp["experiment_id"][ExpId]["required_model_components"])
        self.ReqGblAttrs["experiment"          ] = self.dict_exp["experiment_id"][ExpId]["experiment"]
        self.ReqGblAttrs["experiment_id"       ] = self.dict_exp["experiment_id"][ExpId]["experiment_id"]

        self.ReqGblAttrs["sub_experiment_id"   ] = " ".join(self.dict_exp["experiment_id"][ExpId]["sub_experiment_id"])
        self.ReqGblAttrs["sub_experiment"      ] = self.dict_sub["sub_experiment_id"][self.ReqGblAttrs["sub_experiment_id"]]

        self.ReqGblAttrs["parent_activity_id"  ] = " ".join(self.dict_exp["experiment_id"][ExpId]["parent_activity_id"])
        self.ReqGblAttrs["parent_experiment_id"] = " ".join(self.dict_exp["experiment_id"][ExpId]["parent_experiment_id"])

        if "parent_mip_era" in self.dict_exp["experiment_id"][ExpId].keys():
            self.ReqGblAttrs["parent_mip_era"] = self.dict_exp["experiment_id"][ExpId]["parent_mip_era"]
        else:
            #self.ReqGblAttrs["parent_mip_era"] = "CMIP6"
            self.ReqGblAttrs["parent_mip_era"] = "no parent"

        # need to consider again
        if UserInput["grid_label"] in self.dict_grd["grid_label"].keys():
            self.ReqGblAttrs["grid_label"] = UserInput["grid_label"]
            self.ReqGblAttrs["grid"      ] = self.dict_grd["grid_label"][UserInput["grid_label"]] + ":ne30np4 to 1 by 1 degree"
        else:
            raise ValueError


        if UserInput["institution_id"] in self.dict_ins["institution_id"].keys():
            self.ReqGblAttrs["institution_id"] = UserInput["institution_id"]
            self.ReqGblAttrs["institution"   ] = self.dict_ins["institution_id"][UserInput["institution_id"]]
        else:
            raise ValueError

        if UserInput["nominal_resolution"] in self.dict_nom["nominal_resolution"]:
            self.ReqGblAttrs["nominal_resolution"] = UserInput["nominal_resolution"]
        else:
            raise ValueError

        self.ReqGblAttrs["license"] = self.dict_lic["license"][0]
        self.ReqGblAttrs["license"] = self.ReqGblAttrs["license"].replace("<Your Centre Name>", self.ReqGblAttrs["institution_id"])
        #self.ReqGblAttrs["license"] = self.ReqGblAttrs["license"].replace("<some URL maintained by modeling group>", UserInput["further_info_url"])
        self.ReqGblAttrs["license"] = self.ReqGblAttrs["license"].replace("[ and at <some URL maintained by modeling group>].", ".")
        #we choose sharealike licence
        self.ReqGblAttrs["license"] = self.ReqGblAttrs["license"].replace("[NonCommercial-]", "")



        if UserInput["source_id"] in self.dict_sid["source_id"].keys():
            self.ReqGblAttrs["source_id"] = UserInput["source_id"]

            if UserInput["institution_id"] not in self.dict_sid["source_id"][UserInput["source_id"]]["institution_id"]:
                raise "Error institution is not in the source id"
        else:
            raise ValueError

        srcrec = self.dict_sid["source_id"][UserInput["source_id"]]

        mysource = "{} ({}):".format(srcrec["source_id"], srcrec["release_year"])

        srcstring = mysource
        for mc in srcrec["model_component"].keys():
            desc = srcrec["model_component"][mc]["description"]

            if desc != "none":
               mn = re.findall('([A-Z0-9-]{3,16})', desc)[0]
               if '-' in mn:
                   mn = re.findall('([A-Z0-9-]{3,16}\w+)', desc)[0]
               if len(re.findall('(v[0-9]+.[0-9]+)', desc)) > 0: 
                  mv = re.findall('(v[0-9]+.[0-9]+)', desc)[0]
                  mvfx = mv.replace('.','_')
               else:
                  mvfx = ''
                  mv = ''

               descnew = re.sub('^[,;]', '', desc.replace(mn, '').replace('(','').replace(')','').replace(mv, '').strip()).strip()
               #print (mn, mvfx, descnew)
               #print (mn)

               srcstring += " {}: {} ({}{}, {});".format(mc, mn, mn, mv, descnew)
        #print (srcstring)
            #if srcrec["model_component"][mc]["description"]:
        self.ReqGblAttrs["source"] = srcstring

if __name__ == "__main__":



   cv = cmip6cv()

   cv.GetGlobalAttributes({}, {})

   print (cv.ReqGblAttrs)


   sys.exit()

   print (cv.dict_exp["experiment_id"][ExperimentId])

   GblAttrs={}
   for gk in cv.dict_gbl['required_global_attributes']:
       GblAttrs[gk] = "None"
       for attr in members:
           cvdict = getattr(cv, attr)
           if gk in cvdict.keys():
               #print (gk, cvdict[gk])
               GblAttrs[gk] = attr


   for gk, gval in GblAttrs.items():

       print (gk, gval)

   sys.exit()
   
   
   print (dict_gbl['required_global_attributes'])
   req_gbl_attrs=dict_gbl['required_global_attributes']
   for att in req_gbl_attrs:
       if att == 'activity_id':
          print(dict_act[att][mipname])
   
       if att == 'frequency':
          for frq in dict_frq[att].keys():
             if frq in miptable:
                print (frq)
   
       if att == 'experiment':
          expid = dict_exp[att + '_id'].keys()
          for eid in expid:
   
              #print (dict_exp[att + '_id'][eid]['activity_id'])
              if mipname in dict_exp[att + '_id'][eid]['activity_id']:
                 pass
                 #-print (eid)
                 #-print (dict_exp[att + '_id'][eid]['experiment'])
   
          #print ((dict_ex[patt + '_id'].values())['activity_id'])
       if att == 'institution_id' or att == 'institution':
          print (dict_ins['institution_id'][institute])
   
       if att == 'license':
          print (dict_lic['license'][0].replace('<Your Centre Name>', institute).replace('<some URL maintained by modeling group>', 'https://www.bgc-feedbacks.org'))
          #print (dict_lic['license'][0].replace('<Your Centre Name>', institute))
   
       if att == 'nominal_resolution':
   
          print (dict_nom['nominal_resolution'])
        
          nomgrid=[]
          for nom in dict_nom['nominal_resolution']:
              if 'degree' in nom.split(' ')[1]:
                 nomgrid.append(110)
              else:
                 print (float(nom.split(' ')[0]) - gridsize)
                 nomgrid.append(float(nom.split(' ')[0]))
          print( nomgrid.index(min(nomgrid, key=lambda x:abs(x-gridsize))))
   
          print(dict_nom['nominal_resolution'][nomgrid.index(min(nomgrid, key=lambda x:abs(x-gridsize)))])
       print (att)
   
   
   #print (dict_frq)

