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
        self.objId = {
            'LooseIso'     : lambda row: bool(row.tLooseIso    ) ,
            'MediumIso'    : lambda row: bool(row.tMediumIso   ) ,
            'TightIso'     : lambda row: bool(row.tTightIso    ) ,
            'LooseMVAIso'  : lambda row: bool(row.tLooseMVAIso ) ,
            'MediumMVAIso' : lambda row: bool(row.tMediumMVAIso) ,
            'TightMVAIso'  : lambda row: bool(row.tTightMVAIso ) ,
            }


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
            mcCorrectors.get_muon_corrections(row,'m')

    def sign_cut(self, row):
        return not row.m_t_SS

    def is_MT_Low(self, row):
        return row.mMtToMET < 20

    def muon_id(self, row):
        return bool(row.mPFIDTight) and ( row.mRelPFIsoDB < 0.1 or (row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479))

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.isoMuPass:                     return False
        if not selections.muSelection(row, 'm'):  return False
        if not self.muon_id(row):                 return False    
        if not selections.tauSelection(row, 't'): return False
        if not selections.vetos(row):             return False
        if not row.tAntiElectronLoose:            return False
        if not row.tAntiMuonTight:                return False
        return True
