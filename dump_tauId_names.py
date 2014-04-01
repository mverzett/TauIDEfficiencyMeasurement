#! /bin/env python

import ROOT
from DataFormats.FWLite import Events, Handle
ROOT.gROOT.SetBatch()

events = Events('/hdfs/store/user/tapas/DoubleMu/Run2012A-13Jul2012-v1/AOD/2013-04-01-8TeV-53X-PatTuple_Master/patTuple_cfg-0048B245-B9D2-E111-A8DC-0018F3D096BE.root')
handle = Handle('std::vector<pat::Tau>')
label  = "cleanPatTaus"

for evt in events:
    evt.getByLabel(label, handle)
    taus = handle.product()
    for tauid in taus[0].tauIDs():
        print tauid.first
    break
        
