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

################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################


class TauEffZMT(TauEffBase):
    tree = 'mt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(TauEffZMT, self).__init__(tree, outfile,MuTauTree.MuTauTree, **kwargs)
        self.pucorrector = mcCorrectors.make_puCorrector('singlemu')
        self.objId = [
            'LooseIso'    ,
            'MediumIso'   ,
            'TightIso'    ,
            'LooseMVAIso' ,
            'MediumMVAIso',
            'TightMVAIso' ,
            ]
        self.systematics  = ['NOSYS',"RAW"]+[i+j for i,j in itertools.product(['mes','tes','jes','ues'],['_p'])]#,'_m' add _m if needed also scaled down
        self.id_functions = {
            'LooseIso'      : lambda row: bool(row.tLooseIso    ) ,
            'MediumIso'     : lambda row: bool(row.tMediumIso   ) ,
            'TightIso'      : lambda row: bool(row.tTightIso    ) ,
            'LooseMVAIso'   : lambda row: bool(row.tLooseMVAIso ) ,
            'MediumMVAIso'  : lambda row: bool(row.tMediumMVAIso) ,
            'TightMVAIso'   : lambda row: bool(row.tTightMVAIso ) ,
            'sign_cut'      : self.sign_cut,
            'muon_id'       : self.muon_id,
            'is_mu_anti_iso': self.is_mu_anti_iso,
            }
        
        self.id_functions_with_sys = {
            'is_MT_Low'  : self.is_MT_Low,
            'is_MT_really_high': self.is_MT_really_high,
            }

    def build_folder_structure(self):
        flag_map = {}
        for systematic in self.systematics:
            for obj_id_name in self.objId:
                for sign in ['ss', 'os']:
                    for mt in ['HiMT','LoMT']:
                        flag_map['/'.join((systematic,obj_id_name, sign, mt))] = {
                            obj_id_name : True            ,
                            'sign_cut'  : (sign == 'os')  ,
                            'is_MT_Low' : (mt   == 'LoMT'),
                            'muon_id'   : True            ,
                            }

                    #Make two region QCD dominated (anti-isolate the muon)
                    flag_map['/'.join((systematic,obj_id_name, sign, 'QCD'))] = {
                            obj_id_name      : True          ,
                            'sign_cut'       : (sign == 'os'),
                            'is_MT_Low'      : True          ,
                            'is_mu_anti_iso' : True,
                            }
                        
                #make a sideband with very high mt > 70 to have a pure sample of w+jets
                flag_map['/'.join((systematic,obj_id_name, 'os', 'VHiMT'))] = {
                    obj_id_name : True,
                    'sign_cut'    : True,
                    'is_MT_really_high' : True,
                    }

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
        self.book(directory, "mMtToMET", "Muon 1 eta", 120, 0, 120)
        self.book(directory, "tAbsEta", "Muon 2 eta", 100, 0, 5)
        self.book(directory, "m_t_Mass", "Muon 1-2 Mass", 300, 0, 300)

        # Vetoes
        self.book(directory, 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'muVetoPt5', 'Number of extra muons', 5, -0.5, 4.5)
        self.book(directory, 'tauVetoPt20', 'Number of extra taus', 5, -0.5, 4.5)
        self.book(directory, 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
            
    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m') * \
            mcCorrectors.muon_pog_Mu17Mu8_Mu17_2012(row.mPt, row.mEta)

    def sign_cut(self, row):
        return not row.m_t_SS

    def is_MT_Less_than(self, row, value):
        if self.currect_systematic == 'NOSYS' or self.currect_systematic == '':
            return row.mMtToPfMet_Ty1 < value
        elif self.currect_systematic == 'RAW':
            return row.mMtToMET < value
        elif self.currect_systematic == 'mes_p':
            return row.mMtToPfMet_mes < value
        elif self.currect_systematic == 'tes_p':
            return row.mMtToPfMet_tes < value
        elif self.currect_systematic == 'jes_p':
            return row.mMtToPfMet_jes < value
        elif self.currect_systematic == 'ues_p':
            return row.mMtToPfMet_ues < value
        else:
            raise KeyError("the current systematic, %s is not recognized" % self.currect_systematic)

    def is_MT_Low(self, row):
        return self.is_MT_Less_than(row, 20)

    def is_MT_really_high(self, row):
        return not self.is_MT_Less_than( row, 70)

    def muon_id(self, row):
        return ( row.mRelPFIsoDB < 0.1 or (row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479))

    def is_mu_anti_iso(self, row):
        return row.mRelPFIsoDB > 0.2 and row.mRelPFIsoDB < 0.5

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.mMatchesIsoMu24eta2p1:        return False
        if not selections.muSelection(row, 'm'):  return False
        if not bool(row.mPFIDTight):              return False
            #if not self.muon_id(row):                 return False    
        if not selections.tauSelection(row, 't'): return False
        if not selections.vetos(row):             return False
        if not row.tAntiElectronLoose:            return False
        if not row.tAntiMuonTight:                return False
        return True
