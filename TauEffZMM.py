'''

Make plots of Z->tau tau -> mu tau events.

Author: Mauro Verzetti, UniZh

'''

import MuMuTree
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

class TauEffZMM(TauEffBase):
    tree = 'mm/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(TauEffZMM, self).__init__(tree, outfile,MuMuTree.MuMuTree, **kwargs)
        self.pucorrector = mcCorrectors.make_puCorrector('singlemu')
        self.objId = {
            'h2Tau' : lambda row: bool(row.m2PFIDTight and ( row.m2RelPFIsoDB < 0.1 or (row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479)))
            }


    def book_histos(self, directory):
        self.book(directory, "weight", "Event weight", 100, 0, 5)
        self.book(directory, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(directory, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        #self.book(directory, "prescale", "HLT prescale", 21, -0.5, 20.5)
        self.book(directory, "m1Pt", "Muon 1 Pt", 100, 0, 100)
        self.book(directory, "m2Pt", "Muon 2 Pt", 100, 0, 100)
        self.book(directory, "m1MtToMET", "Muon 1 eta", 120, 0, 120)
        self.book(directory, "m1AbsEta", "Muon 1 eta", 100, 0, 5)
        self.book(directory, "m2AbsEta", "Muon 2 eta", 100, 0, 5)
        self.book(directory, "m1_m2_Mass", "Muon 1-2 Mass", 300, 0, 300)

        # Vetoes
        self.book(directory, 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'muVetoPt5', 'Number of extra muons', 5, -0.5, 4.5)
        self.book(directory, 'tauVetoPt20', 'Number of extra taus', 5, -0.5, 4.5)
        self.book(directory, 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
            

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU)*\
            mcCorrectors.get_muon_corrections(row,'m1','m2')

    def sign_cut(self, row):
        return not row.m1_m2_SS

    def is_MT_Low(self, row):
        return row.m1MtToMET < 20

    def muon1_id(self, row):
        return bool(row.m1PFIDTight and ( row.m1RelPFIsoDB < 0.1 or (row.m1RelPFIsoDB < 0.15 and row.m1AbsEta < 1.479)))

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.isoMuPass:                         return False
        if not selections.muSelection(row, 'm1'):     return False
        if not selections.muSelection(row, 'm2', 10): return False
        if not self.muon1_id(row):                    return False    
        if not selections.vetos(row):                 return False
        return True
