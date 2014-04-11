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
        self.objId = [
            'h2Tau', 
            ]
        self.id_functions = {
            'h2Tau'    : lambda row: selections.mu_idIso(row, 'm2'),
            'sign_cut' : self.sign_cut,
            }
        
        self.id_functions_with_sys = {
            }
        self.hfunc['MET_Z_perp'] = lambda row, weight: (row.type1_pfMetEt*ROOT.TMath.Cos(row.m1_m2_ToMETDPhi_Ty1), weight)
        self.hfunc['MET_Z_para'] = lambda row, weight: (row.type1_pfMetEt*ROOT.TMath.Sin(row.m1_m2_ToMETDPhi_Ty1), weight)

    def build_folder_structure(self):
        flag_map = {}
        for obj_id_name in self.objId:
            for sign in ['ss', 'os']:
                flag_map['/'.join((obj_id_name, sign))] = {
                    'h2Tau'    : True,          
                    'sign_cut' : (sign == 'os'),
                    }
        return flag_map

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

        #MET "certification"
        self.book(directory, "type1_pfMetEt" , "type1_pfMetEt" , 200, 0, 200)
        self.book(directory, "type1_pfMetPhi", "type1_pfMetPhi", 100, -7, 7)
        self.book(directory, "MET_Z_perp", "MET_Z_perp", 200, 0, 200)
        self.book(directory, "MET_Z_para", "MET_Z_para", 200, 0, 200)

        # Vetoes
        self.book(directory, 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
        self.book(directory, 'muVetoPt5', 'Number of extra muons', 5, -0.5, 4.5)
        self.book(directory, 'tauVetoPt20Loose3HitsVtx', 'Number of extra taus', 5, -0.5, 4.5)
        self.book(directory, 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
            

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        trigger_weight  = 1.
        trigger_weight *= mcCorrectors.muon_pog_IsoMu24eta2p1(row.m1Pt, row.m1Eta) \
            if row.m1MatchesIsoMu24eta2p1 else \
            1.
        trigger_weight *= mcCorrectors.muon_pog_IsoMu24eta2p1(row.m2Pt, row.m2Eta) \
            if row.m2MatchesIsoMu24eta2p1 else \
            1.
        return self.pucorrector(row.nTruePU)*\
            mcCorrectors.muon_pog_PFTight(row.m1Pt, row.m1Eta) * \
            mcCorrectors.muon_pog_Iso(row.m1Pt, row.m1Eta) * \
            mcCorrectors.muon_pog_PFTight(row.m2Pt, row.m2Eta) * \
            mcCorrectors.muon_pog_Iso(row.m2Pt, row.m2Eta) * \
            trigger_weight

    def sign_cut(self, row):
        return not row.m1_m2_SS

    def is_MT_Low(self, row):
        return True #row.m1MtToMET < 20 #This cut is not used in this channel

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.isoMu24eta2p1Pass                : return False
        if not selections.muSelection(row, 'm1', 10): return False
        if not selections.muSelection(row, 'm2', 10): return False
        m1matches = row.m1MatchesIsoMu24eta2p1 and \
                    (row.m1Pt > 25) and (row.m1AbsEta < 2.1) 
        m2matches = row.m2MatchesIsoMu24eta2p1 and \
                    (row.m2Pt > 25) and (row.m2AbsEta < 2.1) 
        if not (m1matches or m2matches):       return False
        if not selections.mu_idIso(row, 'm1'): return False    
        if not selections.vetos(row):          return False
        return True
