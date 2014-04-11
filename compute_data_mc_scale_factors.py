#! /bin/env python

import glob
import os
import re
import json
from uncertainties import ufloat
from pprint import pprint
from FinalStateAnalysis.Utilities.struct import struct
from FinalStateAnalysis.StatTools.quad import quad

########################################################################################################################################################
######################################                                                                            ######################################
######################################                         General Purpose Functions                          ######################################
######################################                                                                            ######################################
########################################################################################################################################################

def convert(input):
    '''converts json unicode strings into bytecode strings'''
    if isinstance(input, dict):
        return dict([(convert(key), convert(value)) for key, value in input.iteritems()])
    elif isinstance(input, list):
        return [convert(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input  


def dict_ufloat(dictionary):
    ret = {}
    for name, values in dictionary.iteritems():
        value     = values['value']
        stat_rel  = ufloat(1., values['stat'] / values['value'], name+'_stat')
        sys_rel   = reduce(
            lambda x, y: x*y,
            [ ufloat(
                1.,
                sys_val / values['value'],
                name+'_'+sysname)
                for sysname, sys_val in values.iteritems()
                if sysname.startswith('sys_')
                ],
            1.
            )
        ret[name] = value*stat_rel*sys_rel
    return ret

def get_err(value, tag):
    #print value.error_components()
    return quad(
        *[ j for i, j in value.error_components().iteritems() if tag in i.tag]
        )
    

########################################################################################################################################################
######################################                                                                            ######################################
######################################                        Loads summary tables from json                      ######################################
######################################                                                                            ######################################
########################################################################################################################################################

jobid              = os.environ['jobid']
mumu_table_file    = [i for i in glob.glob('results/%s/plots/mm/summary_table*.json' % jobid) if 'json' in i][0]
mutau_table_files  = [i for i in glob.glob('results/%s/plots/mt/*/summary_table*.json' % jobid) if 'json' in i]
output_dir         = 'results/%s/plots/' % jobid
mumu_table         = dict_ufloat(
    json.loads(
        open(mumu_table_file).read(),
        object_hook=convert)
        )

iso_regex          = re.compile(r'.+\/mt\/(?P<iso_name>\w+)/\w+')
mutau_tables       = dict(
    [ (
        iso_regex.match(i).group('iso_name'), #key of a super-dictionary --> iso name
        dict_ufloat( #Makes ufloat dictionary
                     json.loads( #Loads json
                                 open(i).read(), #Txt to load from
                                 object_hook=convert
                                )
            )
        )
     for i in mutau_table_files
     ]
    )

for iso in mutau_tables:
    mutau_tables[iso]['QCD'] *= ufloat(1, 0.1, 'qcd_extrapolation_stat')

########################################################################################################################################################
######################################                                                                            ######################################
######################################                           Output table functions                           ######################################
######################################                                                                            ######################################
########################################################################################################################################################

columns            = ('process','Obs. Events','Stat. Err','Exp. Evesnts','Stat. Err.','Sys. Error','Ratio Data/MC (%)', 'Stat. Error', 'Sys. Error','Tau Efficiency', 'Stat. Error', 'Sys. Error')
cell_spacing       = '%'+str(max([min([int(len(i)*1.5) for i in columns]), 20]))+'s' # % 
process_spacing    = '%30s'
line_format        = process_spacing+cell_spacing*(len(columns)-1)+'\n'
separator          = '-'*(20*(len(columns)-1) + 30)+'\n'
toprint            = line_format % columns
toprint           += separator

def events_format(entry, isData=True, form='%.1f', multiplier=1.):
    value = entry.nominal_value
    stat  = get_err(entry, 'stat')
    sys   = get_err(entry, 'sys')
    if isData:
        return  (form % (value*multiplier),
                 form % (stat*multiplier))
    else:
        return  form % (value*multiplier), form % (stat*multiplier), form % (sys*multiplier)
    
def get_table_entry(proc_name, dictionary, ratio, efficiency=None):
    data = dictionary['data']
    mc   = dictionary['bkg_sum']
    entries = [proc_name]
    entries.extend(events_format(data))
    entries.extend(events_format(mc, False))
    entries.extend(events_format(ratio, False, multiplier=100))
    if efficiency:
        entries.extend(events_format(efficiency, False, multiplier=100))
    else:
        entries.extend([cell_spacing%'-']*3)
    return line_format % tuple(entries)


########################################################################################################################################################
######################################                                                                            ######################################
######################################                                Data / MC ratio                             ######################################
######################################                                                                            ######################################
########################################################################################################################################################

def compute_data_MC_Ratio(dictionary, refsample = 'Z_jets'):
    #pprint(dictionary)
    return dictionary['data'] / dictionary['bkg_sum']

mumu_ratio = compute_data_MC_Ratio(mumu_table)
toprint   += get_table_entry('Zmumu',mumu_table,mumu_ratio)
toprint   += separator

mutau_data_mc_ratios = {}
for iso_name, table in mutau_tables.iteritems():
    data_MC_ratio = compute_data_MC_Ratio(table)
    mutau_data_mc_ratios[iso_name] = data_MC_ratio
    #toprint      += get_table_entry(iso_name,table,data_MC_ratio)

########################################################################################################################################################
######################################                                                                            ######################################
######################################                                Efficiency measurement                      ######################################
######################################                                                                            ######################################
########################################################################################################################################################

data_mumu   = mumu_table['data']
mc_mumu     = mumu_table['Z_jets']
mc_bkg_mumu = sum(
    [ j for i, j in mumu_table.iteritems()
       if i != 'data' and i != 'bkg_sum' and i != 'Z_jets'
       ]
    )

## muon_id     = 1.*ufloat((1,0.001),'muon_id_sys')
## muon_iso    = 0.984*ufloat((1,0.006/0.984),'muon_iso_sys')
## muon_trig   = 0.981*ufloat((1,0.006/0.981),'muon_trg_sys')

for iso_name, table in mutau_tables.iteritems():
    data_mutau   = table['data']
    mc_mutau     = table['Z_tautau']
    mc_bkg_mutau = sum(
        [ j for i, j in table.iteritems()
           if i != 'data' and i != 'bkg_sum' and i != 'Z_tautau'
           ]
        )
    ## print '(data_mutau - mc_bkg_mutau)',(data_mutau - mc_bkg_mutau)
    ## print '(data_mumu - mc_bkg_mumu)', (data_mumu - mc_bkg_mumu)
    ## print '(mc_mumu/mc_mutau)', (mc_mumu/mc_mutau)
    eff = ((data_mutau - mc_bkg_mutau) / (data_mumu - mc_bkg_mumu)) * \
        (mc_mumu/mc_mutau)
       
    ## if iso_name == 'LooseIso':
    ##     for (var, error) in eff.error_components().items():
    ##         print "%s: %f" % (var.tag, error)
    toprint      += get_table_entry(iso_name,table,mutau_data_mc_ratios[iso_name],eff)
    

    


########################################################################################################################################################
######################################                                                                            ######################################
######################################                                  writes table                              ######################################
######################################                                                                            ######################################
########################################################################################################################################################

toprint   += separator
with open(os.path.join(output_dir,'final_table.raw_txt'),'w') as outfile:
    outfile.write(toprint)
