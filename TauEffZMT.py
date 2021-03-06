'''

Make plots of Z->tau tau -> mu tau events.

Author: Mauro Verzetti, UniZh

'''

import MuTauTree
from TauEffBase import TauEffBase
from FinalStateAnalysis.PlotTools.decorators import memo
import mcCorrectors
import baseSelections as selections
import glob
import os
import ROOT
import itertools
import pprint

@memo
def getHist(dire , name):
    return '/'.join([dire,name])

sys_name_mapping = {
    'NOSYS': '',
    ''     : '',
    'mes_p': 'mes',
    'tes_p': 'tes',
    'jes_p': 'jes',
    'ues_p': 'ues',
}

@memo
def getsysattrname(name, sys):
    return name + '_' + sys_name_mapping[sys]

def getsysattr(row, name, sys):
    return getattr(row, getsysattrname(name, sys) )


################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################


class TauEffZMT(TauEffBase):
    tree = 'mt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(TauEffZMT, self).__init__(tree, outfile,MuTauTree.MuTauTree, **kwargs)
        self.pucorrector = mcCorrectors.make_puCorrector('singlemu')

        #Not Used (yet)
        def make_iso_functor(name):
            attr = 't'+name
            def _f(row):
                return bool( getattr(row, attr) )
            return _f

        self.objId = [
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

        self.systematics  = ['NOSYS',"RAW"]+[i+j for i,j in itertools.product(['mes','tes','jes','ues'],['_p'])]#,'_m' add _m if needed also scaled down ['NOSYS']#['NOSYS']#
        self.id_functions = {
            'VLooseIso' 	            : lambda row: bool(row.tDecayFinding) and bool(row.tVLooseIso),
            'LooseIso'  	            : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso), 
            'MediumIso' 	            : lambda row: bool(row.tDecayFinding) and bool(row.tMediumIso),
            'TightIso'  	            : lambda row: bool(row.tDecayFinding) and bool(row.tTightIso), 
            'LooseIso3Hits'                 : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits),	
            'LooseIso3HitsAntiEleLoose'     : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits) and bool(row.tAntiElectronLoose),
            'LooseIso3HitsAntiEleMVAVLoose' : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits) and bool(row.tAntiElectronMVA5VLoose),
            'LooseIso3HitsAntiEleMVALoose'  : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits) and bool(row.tAntiElectronMVA5Loose), 
            'LooseIso3HitsAntiEleMVAMedium' : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits) and bool(row.tAntiElectronMVA5Medium),
            'LooseIso3HitsAntiEleMVATight'  : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits) and bool(row.tAntiElectronMVA5Tight), 
            'LooseIso3HitsAntiMuon3Tight'   : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits) and bool(row.tAntiMuon3Tight),
            'LooseIso3HitsAntiMuonMVATight' : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIso3Hits) and bool(row.tAntiMuonMVATight),
            'MediumIso3Hits' 	            : lambda row: bool(row.tDecayFinding) and bool(row.tMediumIso3Hits),
            'TightIso3Hits'                 : lambda row: bool(row.tDecayFinding) and bool(row.tTightIso3Hits), 
            'VLooseIsoMVA3OldDMNoLT'        : lambda row: bool(row.tDecayFinding) and bool(row.tVLooseIsoMVA3OldDMNoLT), 	
            'LooseIsoMVA3OldDMNoLT'         : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIsoMVA3OldDMNoLT),  	
            'MediumIsoMVA3OldDMNoLT'        : lambda row: bool(row.tDecayFinding) and bool(row.tMediumIsoMVA3OldDMNoLT), 	
            'TightIsoMVA3OldDMNoLT'         : lambda row: bool(row.tDecayFinding) and bool(row.tTightIsoMVA3OldDMNoLT),  	
            'VTightIsoMVA3OldDMNoLT'        : lambda row: bool(row.tDecayFinding) and bool(row.tVTightIsoMVA3OldDMNoLT), 	
            'VVTightIsoMVA3OldDMNoLT'       : lambda row: bool(row.tDecayFinding) and bool(row.tVVTightIsoMVA3OldDMNoLT),	
            'VLooseIsoMVA3OldDMLT'          : lambda row: bool(row.tDecayFinding) and bool(row.tVLooseIsoMVA3OldDMLT), 
            'LooseIsoMVA3OldDMLT'           : lambda row: bool(row.tDecayFinding) and bool(row.tLooseIsoMVA3OldDMLT),  
            'MediumIsoMVA3OldDMLT'          : lambda row: bool(row.tDecayFinding) and bool(row.tMediumIsoMVA3OldDMLT), 
            'TightIsoMVA3OldDMLT'           : lambda row: bool(row.tDecayFinding) and bool(row.tTightIsoMVA3OldDMLT),  
            'VTightIsoMVA3OldDMLT'          : lambda row: bool(row.tDecayFinding) and bool(row.tVTightIsoMVA3OldDMLT), 
            'VVTightIsoMVA3OldDMLT'         : lambda row: bool(row.tDecayFinding) and bool(row.tVVTightIsoMVA3OldDMLT),            'sign_cut'      : self.sign_cut,
            'muon_id'       : self.muon_id,
            'is_mu_anti_iso': self.is_mu_anti_iso,
            }
        
        self.id_functions_with_sys = {
            'HiMT'     : self.HiMT    ,
            'LoMT'     : self.LoMT    , 
            'VHiMT'    : self.VHiMT   , 
            'MT70_120' : self.MT70_120,
            'MTLt40'   : self.MTLt40  ,
            }


    def build_folder_structure(self):
        flag_map = {}
        systematics = self.systematics if not os.environ['megatarget'].startswith('data_') else self.systematics[0] #for data no shifting sys applied, save time!
        for systematic in self.systematics:
            for obj_id_name in self.objId + ['QCD']:
                for sign in ['ss', 'os']:
                    for mt in ['HiMT', 'LoMT', 'MTLt40', 'VHiMT', 'MT70_120']:
                        region_name = '/'.join((systematic,obj_id_name, sign, mt))
                        region_dict = { 
                            'sign_cut' : (sign == 'os'), 
                            mt         : True,
                        }
                        if obj_id_name == 'QCD':
                            region_dict['LooseIso'] = False
                            region_dict['is_mu_anti_iso'] = True
                        else:
                            region_dict[obj_id_name] = True
                            region_dict['muon_id']   = True
                        flag_map[region_name] = region_dict
        return flag_map


    def book_histos(self, directory):
        self.book(directory, "weight", "Event weight", 100, 0, 5)
        #self.book(directory, "weight_nopu", "Event weight without PU", 100, 0, 5)
        self.book(directory, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(directory, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        #self.book(directory, "prescale", "HLT prescale", 21, -0.5, 20.5)
        self.book(directory, "mPt", "Muon 1 Pt", 100, 0, 100)
        self.book(directory, "tPt", "Muon 2 Pt", 100, 0, 100)
        self.book(directory, "mAbsEta", "Muon 1 eta", 100, 0, 5)
        self.book(directory, "mMtToPfMet_Ty1", " ", 200, 0, 400)
        self.book(directory, "tAbsEta", "Muon 2 eta", 100, 0, 5)
        self.book(directory, "m_t_Mass", "Muon 1-2 Mass", 150, 0, 150)

        # Vetoes
        self.book(directory, 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'muVetoPt5', 'Number of extra muons', 5, -0.5, 4.5)
        self.book(directory, 'tauVetoPt20Loose3HitsVtx', 'Number of extra taus', 5, -0.5, 4.5)
        self.book(directory, 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)

        if 'TauSpinned' in os.environ['megatarget']:
            #print 'booking tauSpinnerWeight'
            self.book(directory, 'tauSpinnerWeight', 'Tau spinner weight', 250, 0, 5)
            
    def event_weight(self, row):
        if row.run > 2:
            return 1.
        weight = row.tauSpinnerWeight if 'TauSpinned' in os.environ['megatarget'] else 1.

        return weight *\
            self.pucorrector(row.nTruePU) * \
            mcCorrectors.muon_pog_PFTight(row.mPt, row.mEta) * \
            mcCorrectors.muon_pog_Iso(row.mPt, row.mEta) * \
            mcCorrectors.muon_pog_IsoMu24eta2p1(row.mPt, row.mEta)

    def sign_cut(self, row):
        return not row.m_t_SS

    def is_MT_Less_than(self, row, value, currect_systematic):
        if currect_systematic == 'NOSYS' or currect_systematic == '':
            return row.mMtToPfMet_Ty1 < value
        elif currect_systematic == 'RAW':
            return row.mMtToMET < value
        elif currect_systematic == 'mes_p':
            return row.mMtToPfMet_mes < value
        elif currect_systematic == 'tes_p':
            return row.mMtToPfMet_tes < value
        elif currect_systematic == 'jes_p':
            return row.mMtToPfMet_jes < value
        elif currect_systematic == 'ues_p':
            return row.mMtToPfMet_ues < value
        else:
            raise KeyError("the current systematic, %s is not recognized" % currect_systematic)

    def pZeta(self, row, value, currect_systematic):
        if currect_systematic == 'NOSYS' or currect_systematic == '':
            return row.mMtToPfMet_Ty1 < value
        elif currect_systematic == 'RAW':
            return row.mMtToMET < value
        elif currect_systematic == 'mes_p':
            return row.mMtToPfMet_mes < value
        elif currect_systematic == 'tes_p':
            return row.mMtToPfMet_tes < value
        elif currect_systematic == 'jes_p':
            return row.mMtToPfMet_jes < value
        elif currect_systematic == 'ues_p':
            return row.mMtToPfMet_ues < value
        else:
            raise KeyError("the current systematic, %s is not recognized" % currect_systematic)

    def LoMT(self, row, sys):
        return getsysattr(row, 'mMtToPfMET', sys) < 20

    def HiMT(self, row, currect_systematic):
        return getsysattr(row, 'mMtToPfMET', sys) >= 20

    def VHiMT(self, row, currect_systematic):
        return getsysattr(row, 'mMtToPfMET', sys) >= 70

    def MT70_120(self, row, currect_systematic):
        return 70 <= getsysattr(row, 'mMtToPfMET', sys) < 120

    def MTLt40(self, row, currect_systematic):
        return getsysattr(row, 'mMtToPfMET', sys) < 40

    def muon_id(self, row):
        return selections.mu_idIso(row, 'm') 

    def is_mu_anti_iso(self, row):
        return row.mRelPFIsoDBDefault > 0.2 and row.mRelPFIsoDBDefault < 0.5

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not (row.isoMu24eta2p1Pass and \
                row.mMatchesIsoMu24eta2p1):          return False #if not row.mMatchesIsoMu24eta2p1:         return False #
        if not selections.muSelection(row, 'm'):     return False
        if not (row.mPt > 25 and row.mAbsEta < 2.1): return False
        if not bool(row.mPFIDTight):                 return False
        if not selections.tauSelection(row, 't'):    return False
        if row.tMuOverlap > 0:                       return False
        if not selections.vetos(row):                return False
        if row.m_t_Mass > 150:                    return False

        #separate Z->tautau Z->mumu
        if 'Zjets' in os.environ['megatarget']:
            if 'ZToMuMu' in os.environ['megatarget']:
                if row.isGtautau or row.isZtautau: return False
            else:
                if not (row.isGtautau or row.isZtautau):  return False
                    
        return True




