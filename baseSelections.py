from FinalStateAnalysis.PlotTools.decorators import memo

@memo
def getVar(name, var):
    return name+var

def muSelection(row, name, pt_thr=20):
    if getattr( row, getVar(name,'Pt')) < pt_thr:   return False
    if getattr( row, getVar(name,'AbsEta')) > 2.1:  return False
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2: return False
    return True
    ## if not getattr( row, getVar(name,'PixHits')):   return False
    ## if getattr( row, getVar(name,'JetBtag')) > 3.3: return False

    
def tauSelection(row, name):
    if getattr( row, getVar(name,'Pt')) < 20:          return False
    if getattr( row, getVar(name,'AbsEta')) > 2.3:     return False
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2:    return False
    return True

def vetos(row):
    if row.muVetoPt5:         return False
    if row.bjetCSVVeto:       return False
    if row.tauVetoPt20Loose3HitsVtx:  return False
    if row.eVetoCicTightIso:  return False
    return True

def mu_idIso(row, name):
    return bool(getattr(row, getVar(name, 'PFIDTight'))) \
        and bool(getattr(row, getVar(name, 'RelPFIsoDBDefault')) < 0.12)




