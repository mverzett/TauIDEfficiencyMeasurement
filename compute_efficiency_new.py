#! /bin/env python

import glob
import os
import re
import json
from uncertainties import ufloat
from pprint import pprint
from pdb import set_trace
from fnmatch import fnmatch
from FinalStateAnalysis.StatTools.quad import quad
import FinalStateAnalysis.Utilities.prettyjson as prettyjson
import FinalStateAnalysis.Utilities.floatformatting as floatformatting
from FinalStateAnalysis.Utilities.struct import RecursiveStruct

##############
#helper
##############

class DefaultDict(object):
    '''Automatically creates an entry when it's missing'''
    def __init__(self, creator, values = {}):
        self.creator = creator
        self.vals    = values
    
    def keys(self):
        return self.vals.keys()

    def items(self):
        return self.vals.items()

    def values(self):
        return self.vals.values()

    def iteritems(self):
        return self.vals.iteritems()
    
    def __getitem__(self, key):
        if key not in self.vals:
            self.vals[key] = self.creator(key)
        return self.vals[key]
            

########################################################################################################################################################
######################################                             Systematics table                              ######################################
########################################################################################################################################################

sys_map = {
    'wz'    : ['sys_wz_xsec'],
    'zz'    : ['sys_zz_xsec'],
    'ww'    : ['sys_ww_xsec'],
    'ttbar' : ['sys_ttbar_xsec'],
}

sys_map_mt = { 'ztt' : ['sys_tauES']}
sys_map_mt.update(sys_map)

sys_map_mm = {
    ('wz','zz','ww','ttbar','ztt','zmm', 'wjets') : ['sys_eff_m', 'sys_eff_trg'],
}
sys_map_mm.update(sys_map)

sys_value_map = {
    #multiplicative
    'sys_wz_xsec'    : ufloat(1   , quad(0.04, 0.04), 'sys_wz_xsec'   ), #QCD Scale, hadronization
    'sys_zz_xsec'    : ufloat(1   , quad(0.04, 0.04), 'sys_zz_xsec'   ), #QCD Scale, hadronization
    'sys_ww_xsec'    : ufloat(1   , quad(0.04, 0.04), 'sys_ww_xsec'   ), #QCD Scale, hadronization
    'sys_ttbar_xsec' : ufloat(1   , 0.1             , 'sys_ttbar_xsec'),
    'sys_eff_m'      : ufloat(1   , 0.02            , 'sys_eff_m'     ),
    'sys_eff_trg'    : ufloat(1   , 0.01            , 'sys_eff_trg'   ),
    'sys_tauES'      : ufloat(1   , 0.03            , 'sys_tauES'     ),
    'sys_qcdSS/OS'   : ufloat(1.06, 0.1             , 'sys_qcdSS/OS'  ),
    #additive
    'sys_ues' : ufloat(0., 1., 'sys_ues'),
    'sys_jes' : ufloat(0., 1., 'sys_jes'),
    'sys_tes' : ufloat(0., 1., 'sys_tes'),
    'sys_mes' : ufloat(0., 1., 'sys_mes'),
}

sys_groups = {
    'MET'               : re.compile('sys_\wes'    ),
    'QCD Extrapolation' : re.compile('sys_qcdSS/OS'),
    'Tau ES'            : re.compile('sys_tauES'   ),
    'Mu ID/Iso Eff'     : re.compile('sys_eff_m'   ),
    'Trigger Unc.'      : re.compile('sys_eff_trg' ),
    'MC xsections'      : re.compile('sys_\w+_xsec'),
    'MC Stats'          : re.compile('stat_\w+/\w+/\w+/\w+/(?!data)'),
}
data_stat_re = re.compile('stat_.*data')

########################################################################################################################################################
######################################                         General Purpose Functions                          ######################################
########################################################################################################################################################

def apply_systematics(sample, systematics_map):
    ret = 1.
    for sys_key, sys_vals in systematics_map.iteritems():
        matches = fnmatch(sample, sys_key) if isinstance(sys_key, str) else \
                  any(fnmatch(sample, pattern) for pattern in sys_key)
        if matches:
            for sys_tag in sys_vals:
                ret *= sys_value_map[sys_tag]
    return ret

def dict_ufloat(dictionary, name, sample, systematics_map):
    value    = ufloat(dictionary['val'], dictionary['stat'],'stat_'+name)
    value   *= apply_systematics(sample, systematics_map)
    for sys_err, sys_val in dictionary.iteritems():
        if not sys_err.startswith('sys_'):
            continue
        value += sys_value_map[sys_err]*sys_val
    return value

def get_err(value, tag):
    #print value.error_components()
    return quad(
        *[ j for i, j in value.error_components().iteritems() if fnmatch(i.tag, tag)]
        )

def get_err_re(value, regex, invert=False):
    if invert:
        return quad(
            *[ j for i, j in value.error_components().iteritems() if not regex.match(i.tag)]
        )        
    else:
        return quad(
            *[ j for i, j in value.error_components().iteritems() if regex.match(i.tag)]
        )

def format_ufloat(value, _format='%.0f', show_sys=True, errsign='+/-'):
    form = [value.n]
    form.append(get_err_re(value, data_stat_re))
    if show_sys:
        form.append(get_err_re(value, data_stat_re, True))
    form = [_format % i for i in form]
    return errsign.join(form)

def convert_table(table, systematics_map, name = '', sample=''):
    if 'val' in table:
        return dict_ufloat(table, name, sample, systematics_map)
    else:
        ret = {}
        for location, subtab in table.iteritems():
            ret[location] = convert_table(subtab, systematics_map, os.path.join(name,location), location)
        return ret

def dump_unc_breakdown(value, thr=0.3):
    err        = value.s
    eff_thr    = err*thr
    components = value.error_components().items()
    components.sort(reverse=True, key=lambda x: x[1])
    return [(var.tag, val) for var, val in components if val > eff_thr]

def mkdir(path):
    if not os.path.isdir(path):
        os.mkdir(path)

########################################################################################################################################################
##################################                        MT contributions (duplicate from plotter)                      ###############################
########################################################################################################################################################
        
def compute_signal_contribution(yields):
    # uses ABCD method to extrapolate the WJets and QCD yield in SS region. Regions:
    # A - AntiIso objects LowMt
    # B - AntiIso objects HihgMt
    # C - Iso pass LowMt
    # D - Iso pass HighMt'''
    
    qcd_ss_noIso_loMt = yields.muAntiIso.ss.LoMT.data
    qcd_ss_noIso_hiMt = yields.muAntiIso.ss.HiMT.data
    qcd_ratio_lo_hi_mt= qcd_ss_noIso_loMt / qcd_ss_noIso_hiMt

    data_ss_Iso_loMt = yields.muIso.ss.LoMT.data
    data_ss_Iso_hiMt = yields.muIso.ss.HiMT.data
    
    wjet_ss_Iso_loMt = yields.muIso.ss.LoMT.wjets
    wjet_ss_Iso_hiMt = yields.muIso.ss.HiMT.wjets
    w_ratio_lo_hi_mt = wjet_ss_Iso_loMt / wjet_ss_Iso_hiMt

    qcd_yield_hiMt = (data_ss_Iso_hiMt*w_ratio_lo_hi_mt - data_ss_Iso_loMt)/(w_ratio_lo_hi_mt - qcd_ratio_lo_hi_mt)
    qcd_yield_loMt = qcd_yield_hiMt*qcd_ratio_lo_hi_mt

    yields.muIso.ss.LoMT.qcd = qcd_yield_loMt
    yields.muIso.ss.HiMT.qcd = qcd_yield_hiMt

    #add 10% uncertainty on ss/os ratio
    yields.muIso.os.LoMT.qcd = qcd_yield_loMt * sys_value_map['sys_qcdSS/OS']
    yields.muIso.os.HiMT.qcd = qcd_yield_hiMt * sys_value_map['sys_qcdSS/OS']

    #WJets in signal region from HiMT
    wjet_region   = yields.muIso.os.HiMT
    signal_region = yields.muIso.os.LoMT

    wjet_os_Iso_loMt = signal_region.wjets
    wjet_os_Iso_hiMt = wjet_region.wjets
    w_ratio_os       = wjet_os_Iso_loMt / wjet_os_Iso_hiMt

    wjet_estimate = wjet_region.data - (
        wjet_region.ztt   +
        wjet_region.zmm   +
        wjet_region.ttbar +
        wjet_region.qcd   +
        wjet_region.zz    +
        wjet_region.wz    +
        wjet_region.ww   
    )        

    signal_region.wjets_estimate = wjet_estimate * w_ratio_os
    signal_region.bkg_sum = (
        signal_region.wjets_estimate +
        signal_region.qcd   +
        signal_region.zmm   +
        signal_region.ttbar +
        signal_region.wz    +
        signal_region.ww    +
        signal_region.zz    
        )

    signal_region.mc_sum = signal_region.bkg_sum + signal_region.ztt
    signal_region.ztt_estimate = signal_region.data - signal_region.bkg_sum
    signal_region.ztt_ratio  = signal_region.ztt_estimate / signal_region.ztt


def compute_mm_contribution(yields):
    qcd_est = yields.ss.data - (
        yields.ss.wjets +
        yields.ss.ttbar +
        yields.ss.zmm   +
        yields.ss.wz    +
        yields.ss.ww    +
        yields.ss.zz    
        )
    qcd_est *= sys_value_map['sys_qcdSS/OS']

    yields.os.qcd = qcd_est
    yields.os.bkg_sum = (
        yields.os.qcd   +
        yields.os.wjets +
        yields.os.ttbar +
        yields.os.wz    +
        yields.os.ww    +
        yields.os.zz    
        )

    signal_region = yields.os
    signal_region.mc_sum = signal_region.bkg_sum + signal_region.zmm
    signal_region.zmm_estimate = signal_region.data - signal_region.bkg_sum
    signal_region.zmm_ratio  = signal_region.zmm_estimate / signal_region.zmm


########################################################################################################################################################
######################################                        Loads summary tables from json                      ######################################
########################################################################################################################################################

jobid    = os.environ['jobid']
mt_table = convert_table(
    prettyjson.loads(
        open('results/%s/plots/mt/yield_dump.json' % jobid).read()
    ),
    sys_map_mt
)
#convert to struct, easier access
mt_table = RecursiveStruct(**mt_table)

mm_table = convert_table(
    prettyjson.loads(
        open('results/%s/plots/mm/yield_dump.json' % jobid).read()
    ),
    sys_map_mm
)
#convert to struct, easier access
mm_table = RecursiveStruct(**mm_table)

ids =  [
    'VLooseIso' 	           ,
    'LooseIso'  	           ,
    'MediumIso' 	           ,
    'TightIso'  	           ,
    'LooseIso3Hits'                ,
    'LooseIso3HitsAntiEleLoose'    ,
    'LooseIso3HitsAntiEleMVAVLoose',
    'LooseIso3HitsAntiEleMVALoose' ,
    'LooseIso3HitsAntiEleMVAMedium',
    'LooseIso3HitsAntiEleMVATight' ,
    'LooseIso3HitsAntiMuon3Tight'  ,
    'LooseIso3HitsAntiMuonMVATight',
    'MediumIso3Hits' 	           ,
    'TightIso3Hits'                ,
    'VLooseIsoMVA3OldDMNoLT'       ,
    'LooseIsoMVA3OldDMNoLT'        ,
    'MediumIsoMVA3OldDMNoLT'       ,
    'TightIsoMVA3OldDMNoLT'        ,
    'VTightIsoMVA3OldDMNoLT'       ,
    'VVTightIsoMVA3OldDMNoLT'      ,
    'VLooseIsoMVA3OldDMLT'         ,
    'LooseIsoMVA3OldDMLT'          ,
    'MediumIsoMVA3OldDMLT'         ,
    'TightIsoMVA3OldDMLT'          ,
    'VTightIsoMVA3OldDMLT'         ,
    'VVTightIsoMVA3OldDMLT'        ,
    ]

for tid in ids:
    compute_signal_contribution(
        getattr(mt_table, tid)
    )

compute_mm_contribution(mm_table.h2Tau)


########################################################################################################################################################
######################################                                                                            ######################################
######################################                           Output table functions                           ######################################
######################################                                                                            ######################################
########################################################################################################################################################

rows = [('process','Obs. Events','Exp. Events', 'Obs. w/o bkg', 'Exp. Z evts', 'Ratio Data/MC (%)','Tau Efficiency SF')]

#Zmm table
rows.append((
    'Zmm', 
    mm_table.h2Tau.os.data.format('.2e'),    
    format_ufloat(mm_table.h2Tau.os.mc_sum            ),
    format_ufloat(mm_table.h2Tau.os.zmm_estimate      ),
    format_ufloat(mm_table.h2Tau.os.zmm               ),
    format_ufloat(mm_table.h2Tau.os.zmm_ratio, '%.3f' ),
    '--'
))

for tid in ids:
    tid_table = getattr(mt_table, tid)
    signal_region = tid_table.muIso.os.LoMT
    signal_region.eff_ratio = signal_region.ztt_ratio / mm_table.h2Tau.os.zmm_ratio
    rows.append((
        tid, 
        signal_region.data.format('.2e'),    
        format_ufloat(signal_region.mc_sum            ),
        format_ufloat(signal_region.ztt_estimate      ),
        format_ufloat(signal_region.ztt               ),
        format_ufloat(signal_region.ztt_ratio, '%.3f' ),
        format_ufloat(signal_region.eff_ratio, '%.3f' ), 
    ))

transposed = [[len(i) for i in column ] for column in zip(*rows)]
max_space  = [int(max(i)*1.2) for i in transposed]
tot_space = sum(max_space)
separator = '-'*tot_space
formatter = ''.join(['%'+str(i)+'s' for i in max_space])
tex_form  = ' &'.join(['%'+str(i)+'s' for i in max_space])+r' \\'
tex_sep   = r'\hline'

def tex_preamble(cols,ttype='table'): 
    return r'''
\begin{%s}
\begin{center}
\begin{tabular}{|%s|}
''' % (ttype, '|'.join(['c']*cols))

tex_end = r'''
\end{tabular}
\end{center}
\end{%s}
'''


with open('results/%s/plots/efficiencies.raw_txt' % jobid,'w') as outfile:
    outfile.write( separator+'\n')
    outfile.write( (formatter % rows[0])+'\n')
    outfile.write( separator+'\n')
    outfile.write( (formatter % rows[1])+'\n')
    outfile.write( separator+'\n')
    for row in rows[2:]:
        outfile.write( formatter % row+'\n')
    outfile.write( separator+'\n')

with open('results/%s/plots/efficiencies.tex' % jobid,'w') as outfile:
    outfile.write( tex_preamble(len(rows[0]), 'sidewaystable') )
    outfile.write( tex_sep+'\n')
    outfile.write( (tex_form % rows[0]).replace('%', r'\%')+'\n')
    outfile.write( tex_sep+'\n')
    outfile.write( (tex_form % rows[1]).replace('+/-', '$\\pm$')+'\n')
    outfile.write( tex_sep+'\n')
    for row in rows[2:]:
        outfile.write( (tex_form % row).replace('+/-', '$\\pm$')+'\n')
    outfile.write( tex_sep+'\n')
    outfile.write( tex_end % 'sidewaystable')

###########################################
##       Single eff tables               ##
###########################################

tables_dir = 'results/%s/plots/mt/' % jobid
separator  = '-'*50

signal_samples = [
    'ztt'  ,
    'wjets_estimate',
    'qcd'  , 
    'zmm'  , 
    'ttbar', 
    'wz'   , 
    'ww'   , 
    'zz'   , 
]

def compute_sys_breakdown(value):
    sys_breakdown = []
    for sys_group, regex in sys_groups.iteritems():
        sys_breakdown.append(
            (sys_group, 
             get_err_re(
                 value, 
                 regex
             )
            )
        )
    sys_breakdown.sort(reverse=True, key=lambda x: x[1])
    return sys_breakdown
    
for tid in ids:
    signal_region = getattr(mt_table, tid).muIso.os.LoMT
    sys_breakdown = []
    for sys_group, regex in sys_groups.iteritems():
        sys_breakdown.append(
            (sys_group, 
             get_err_re(
                 signal_region.eff_ratio, 
                 regex
             )
            )
        )
    sys_breakdown.sort(reverse=True, key=lambda x: x[1])
    tex_print = tex_preamble(2)
    to_print = []
    to_print.append(separator)
    tex_print += '\\hline\n'
    to_print.append(tid)
    tex_print += (r'\multicolumn{2}{|c|}{%s} \\' % tid)+'\n'
    to_print.append(format_ufloat(signal_region.eff_ratio, '%.3f' ))
    tex_print += (
        (r'\multicolumn{2}{|c|}{$%s$} \\' % format_ufloat(
            signal_region.eff_ratio, '%.3f' 
        )
     )+'\n').replace('+/-', '\\pm')
    to_print.append(separator)
    tex_print += '\\hline\n'
    to_print.append('Systematics breakdown')
    tex_print += (r'\multicolumn{2}{|c|}{Systematics breakdown} \\')+'\n'
    tex_print += '\\hline\n'
    for label, value in sys_breakdown:
        to_print.append('%20s%8s' % (label, '%.3f' % value))
        tex_print += (r'%20s & $%s$ \\' % (label, '%.3f' % value) +'\n').replace('+/-', '\\pm')
    to_print.append(separator)
    tex_print += '\\hline\n'
    to_print.append('Sample breakdown')
    tex_print += (r'\multicolumn{2}{|c|}{Sample breakdown} \\')+'\n'
    tex_print += '\\hline\n'
    for sample in signal_samples:
        to_print.append(
            '%20s%24s' % (
                sample.replace('_',' '), 
                format_ufloat(
                    getattr(signal_region, sample)
                )
            )
        )
        tex_print += (r'%20s & $%s$ \\' % (
            sample.replace('_',' '), 
            format_ufloat(
                getattr(signal_region, sample)
            )
        )).replace('+/-', '\\pm')
        tex_print += '\n'

    to_print.append(separator)
    tex_print += '\\hline\n'
    to_print.append('')
    tex_print += tex_end % 'table'

    with open(os.path.join(tables_dir, '%s.raw_txt' % tid), 'w') as tid_summary:
        tid_summary.write('\n'.join(to_print))

    with open(os.path.join(tables_dir, '%s.tex' % tid), 'w') as tid_summary:
        tid_summary.write(tex_print)
