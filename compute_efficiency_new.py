#! /bin/env python

import glob
import os
import re
import json
from uncertainties import ufloat
from pprint import pprint
from FinalStateAnalysis.StatTools.quad import quad

########################################################################################################################################################
######################################                         General Purpose Functions                          ######################################
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
        stat_rel  = ufloat((1.,values['stat'] / values['value']),name+'_stat')
        sys_rel   = 1.
        for sys_err, sys_val in values.iteritems():
            if not sys_err.startswith('sys_'):
                continue
            sys_rel *= ufloat((1., sys_val/value), sys_err)
        ret[name] = value*stat_rel*sys_rel
    return ret

def get_err(value, tag):
    #print value.error_components()
    return quad(
        *[ j for i, j in value.error_components().iteritems() if tag in i.tag]
        )


########################################################################################################################################################
######################################                        Loads summary tables from json                      ######################################
########################################################################################################################################################

jobid       = os.environ['jobid']
mumu_table  = dict_ufloat(
        json.loads(
            open('results/%s/plots/mm/full_content_dump.json' % jobid).read(),
            object_hook=convert)
    )
mutau_table = dict_ufloat(
    json.loads(
        open('results/%s/plots/mt/full_content_dump.json' % jobid).read(),
        object_hook=convert)
    )

mutau_table = dict_ufloat( json.loads( open("results/2013-Jul-05-TauIDEff-8TeV/plots/mt/LooseIso/summary_table_m_t_Mass.json").read(), object_hook=convert) )


ids =  [
    'LooseIso'    ,
    'MediumIso'   ,
    'TightIso'    , 
    'LooseMVAIso' ,
    'MediumMVAIso',
    'TightMVAIso' ,
    'LooseIso3Hits',
    'LooseMVA2Iso',
    'MediumIso3Hits',
    'MediumMVA2Iso',
    'TightIso3Hits',
    'TightMVA2Iso',
    'VLooseIso',
    ]
