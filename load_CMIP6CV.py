#!/usr/bin/env python

'''
   Parse all cmip6 CV information to be used by lmt-lwcmor
'''


import json

with open ('../CMIP6_CVs/CMIP6_activity_id.json', 'r') as fact:
     dict_act = json.load(fact)

with open ('../CMIP6_CVs/CMIP6_experiment_id.json', 'r') as fexp:
     dict_exp = json.load(fexp)

with open ('../CMIP6_CVs/CMIP6_frequency.json', 'r') as ffrq:
     dict_frq = json.load(ffrq)

with open ('../CMIP6_CVs/CMIP6_grid_label.json', 'r') as fgrd:
     dict_grd = json.load(fgrd)

with open ('../CMIP6_CVs/CMIP6_license.json', 'r') as flic:
     dict_lic = json.load(flic)

with open ('../CMIP6_CVs/CMIP6_nominal_resolution.json', 'r') as fnom:
     dict_nom = json.load(fnom)

with open ('../CMIP6_CVs/CMIP6_required_global_attributes.json', 'r') as fgbl:
     dict_gbl = json.load(fgbl)

with open ('../CMIP6_CVs/CMIP6_source_id.json', 'r') as fsid:
     dict_sid = json.load(fsid)

with open ('../CMIP6_CVs/CMIP6_source_type.json', 'r') as fsty:
     dict_sty = json.load(fsty)

with open ('../CMIP6_CVs/CMIP6_sub_experiment_id.json', 'r') as fsub:
     dict_sub = json.load(fsub)

with open ('../CMIP6_CVs/CMIP6_table_id.json', 'r') as ftid:
     dict_tid = json.load(ftid)

with open ('../CMIP6_CVs/CMIP6_institution_id.json', 'r') as fins:
     dict_ins = json.load(fins)

mipname='LS3MIP'
miptable='Lmon'
institute='RUBISCO'

gridsize=55





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

