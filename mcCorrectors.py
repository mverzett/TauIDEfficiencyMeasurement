import os
import glob
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight

# Determine MC-DATA corrections
is7TeV = bool('7TeV' in os.environ['jobid'])
print "Is 7TeV:", is7TeV

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions  = {
    'singlemu' : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleMu*pu.root')),
    }
mc_pu_tag                  = 'S6' if is7TeV else 'S10'

####################
#2012 Corrections
###################
muon_pog_IsoMu24eta2p1_2012 = MuonPOGCorrections.make_muon_pog_IsoMu24eta2p1_2012() #trigger
muon_pog_PFTight_2012       = MuonPOGCorrections.make_muon_pog_PFTight_2012()  #ID
muon_pog_Iso_2012           = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2012() #Iso


muon_pog_IsoMu24eta2p1 = muon_pog_IsoMu24eta2p1_2012
muon_pog_PFTight       = muon_pog_PFTight_2012      
muon_pog_Iso           = muon_pog_Iso_2012          

#####################
#  PU Corrections
#####################
def make_puCorrector(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    if dataset in pu_distributions:
        return PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributions[dataset]))
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')


