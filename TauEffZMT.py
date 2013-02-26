'''

Make plots of Z->tau tau -> mu tau events.

Author: Mauro Verzetti, UniZh

'''

import MuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
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

class TauEffZMT(MegaBase):
    tree = 'mt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(TauEffZMT, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuTauTree.MuTauTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.pucorrector = mcCorrectors.make_puCorrector('singlemu')
        self.tauIds = [
            'LooseIso',
            'MediumIso',
            'TightIso',
            'LooseMVAIso',
            'MediumMVAIso',
            'TightMVAIso',
            ]


    def begin(self):
        def book_hist(directory):
            self.book(directory, "weight", "Event weight", 100, 0, 5)
            self.book(directory, "weight_nopu", "Event weight without PU", 100, 0, 5)
            self.book(directory, "rho", "Fastjet #rho", 100, 0, 25)
            self.book(directory, "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(directory, "prescale", "HLT prescale", 21, -0.5, 20.5)
            self.book(directory, "mPt", "Muon 1 Pt", 100, 0, 100)
            self.book(directory, "tPt", "Muon 2 Pt", 100, 0, 100)
            self.book(directory, "mAbsEta", "Muon 1 eta", 100, -2.5, 2.5)
            self.book(directory, "tAbsEta", "Muon 2 eta", 100, -2.5, 2.5)
            self.book(directory, "m_t_Mass", "Muon 1-2 Mass", 240, 60, 120)

            # Vetoes
            self.book(directory, 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
            self.book(directory, 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
            self.book(directory, 'muVetoPt5', 'Number of extra muons', 5, -0.5, 4.5)
            self.book(directory, 'tauVetoPt20', 'Number of extra taus', 5, -0.5, 4.5)
            self.book(directory, 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
            
        for iso in self.tauIds:
            for region in ['lowMT','highMT']:
                book_hist('/'.join([iso,region]))

    def correction(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU)

    def fill_histos(self, directory, row):
        histos = self.histograms
        weight = self.correction(row)
        histos[getHist(directory, 'weight')].Fill(weight)
        histos[getHist(directory, 'weight_nopu')].Fill(self.correction(row))
        histos[getHist(directory, 'rho')].Fill(row.rho, weight)
        histos[getHist(directory, 'nvtx')].Fill(row.nvtx, weight)
        histos[getHist(directory, 'prescale')].Fill(row.doubleMuPrescale, weight)
        histos[getHist(directory, 'mPt')].Fill(row.mPt, weight)
        histos[getHist(directory, 'tPt')].Fill(row.tPt, weight)
        histos[getHist(directory, 'mAbsEta')].Fill(row.mAbsEta, weight)
        histos[getHist(directory, 'tAbsEta')].Fill(row.tAbsEta, weight)
        histos[getHist(directory, 'm_t_Mass')].Fill(row.m_t_Mass, weight)

        histos[getHist(directory, 'bjetVeto')].Fill(row.bjetVeto, weight)
        histos[getHist(directory, 'bjetCSVVeto')].Fill(row.bjetCSVVeto, weight)
        histos[getHist(directory, 'muVetoPt5')].Fill(row.muVetoPt5, weight)
        histos[getHist(directory, 'tauVetoPt20')].Fill(row.tauVetoPt20, weight)
        histos[getHist(directory, 'eVetoCicTightIso')].Fill(row.eVetoCicTightIso, weight)

    def muon_id(self, row):
        return bool(row.mPFIDTight) and ( row.mRelPFIsoDB < 0.1 or (row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479))

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not selections.muSelection(row, 'm'):  return False
        if not self.muon_id(row):                 return False    
        if not selections.tauSelection(row, 't'): return False
        if not selections.vetos(row):             return False
        if not row.tAntiElectronLoose:            return False
        if not row.tAntiMuonTight:                return False
        return True

    def process(self):
        for row in self.tree:
            if not self.preselection(row):
                continue
            region = 'lowMT' if row.mMtToMET < 20 else 'highMT'
            for iso in self.tauIds:
                if bool(getattr(row, selections.getVar('t', iso))):
                    self.fill_histos('/'.join([iso,region]), row)

    def finish(self):
        self.write_histos()
