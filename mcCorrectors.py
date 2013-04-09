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
    'doublemu' : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_DoubleMu*pu.root')),
    'singlemu' : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleMu*pu.root')),
    }
mc_pu_tag                  = 'S6' if is7TeV else 'S10'
# Hack to use S6 weights for the HWW 7TeV sample we use in 8TeV
if 'HWW3l' in os.environ.get('megatarget', 'NOTSET') and not is7TeV:
    print "Using S6_600bins PU weights for HWW3l"
    mc_pu_tag = 'S6_600bins'

#Makes appropriate correction function for electrons or muons according to run period
mu_pog_2011_id             = MuonPOGCorrections.make_muon_pog_PFTight_2011()
mu_pog_2011_iso            = MuonPOGCorrections.make_muon_pog_PFRelIsoDB02_2011()
muon_pog_IsoID             = (lambda pt, eta: mu_pog_2011_id(pt,eta)*mu_pog_2011_iso(pt,eta)) if is7TeV else H2TauCorrections.correct_mu_idiso_2012

#scale factors for single muon trigger
muon_pog_IsoMu24eta2p1_2012 = MuonPOGCorrections.make_muon_pog_IsoMu24eta2p1_2012()

#Double electrons does NOT need scale factors    

def make_puCorrector(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    if dataset in pu_distributions:
        return PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributions[dataset]))
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')
## def force_pu_distribution(kind):
##     pu_corrector = PileupWeight.PileupWeight( kind, *pu_distributions)

def get_muon_corrections(row,*args):
    'makes corrections to iso and id of muons'
    ret = 1.
    for arg in args:
        eta = getattr(row, '%sEta' % arg)
        pt  = getattr(row, '%sPt'  % arg)
        #print muon_pog_IsoID.__name__
        ret *= muon_pog_IsoID( pt, eta)
    return ret



