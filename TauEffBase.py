'''

Generic base class for Z->tau tau // mu mu analysis.


'''

from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import os
import pprint
import ROOT
                        
class TauEffBase(MegaBase):
    def __init__(self, tree, outfile, wrapper, **kwargs):
        super(TauEffBase, self).__init__(tree, outfile, **kwargs)
        # Cython wrapper class must be passed
        self.tree = wrapper(tree)
        self.out = outfile
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.histo_locations = {} #just a mapping of the histograms we have to avoid changing self.histograms indexing an screw other files
        self.hfunc   = { #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in fill_histos later
            'nTruePU' : lambda row, weight: (row.nTruePU,None),
            'weight'  : lambda row, weight: (weight,None) if weight is not None else (1.,None),
            'Event_ID': lambda row, weight: (array.array("f", [row.run,row.lumi,int(row.evt)/10**5,int(row.evt)%10**5] ), None),
            }
        self.objId = {}
        self.systematics = ['']
        self.currect_systematic = ''


    def fill_histos(self, histos, folder, row, weight):
        '''fills histograms'''
        #find all keys matching
        folder_str = folder
        for attr in self.histo_locations[folder_str]:
            value = self.histograms[folder_str+'/'+attr]
            if isinstance(value, ROOT.TH2):
                if attr in self.hfunc:
                    result, weight = self.hfunc[attr](row, weight)
                    r1, r2 = result
                    if weight is None:
                        value.Fill( r1, r2 ) #saves you when filling NTuples!
                    else:
                        value.Fill( r1, r2, weight )
                else:
                    attr1, attr2 = tuple(attr.split('#'))
                    v1 = getattr(row,attr1)
                    v2 = getattr(row,attr2)
                    value.Fill( v1, v2, weight ) if weight is not None else value.Fill( v1, v2 )
            else:
                if attr in self.hfunc:
                    result, weight = self.hfunc[attr](row, weight)
                    if weight is None:
                        value.Fill( result ) #saves you when filling NTuples!
                    else:
                        value.Fill( result, weight )
                else:
                    value.Fill( getattr(row,attr), weight ) if weight is not None else value.Fill( getattr(row,attr) )
        return None

    def count_bjets(self, row):
        return row.bjetCSVVeto
    

    def begin(self):
        # Loop over regions, book histograms
        for folder in self.build_folder_structure():
            self.book_histos(folder) # defined in subclass
        for key in self.histograms:
            charpos  = key.rfind('/')
            location = key[ : charpos]
            name     = key[ charpos + 1 :]
            if location in self.histo_locations:
                self.histo_locations[location].append(name)
            else:
                self.histo_locations[location] = [name]
                #pprint.pprint(self.histo_locations) 

    def process(self):
        # For speed, map the result of the region cuts to a folder path
        # string using a dictionary
        # key = (sign, obj1, obj2, obj3)
        cut_region_map = self.build_folder_structure()

        # Reduce number of self lookups and get the derived functions here
        histos       = self.histograms
        preselection = self.preselection
        id_functions = self.id_functions
        id_functions_with_sys = self.id_functions_with_sys
        fill_histos  = self.fill_histos
        weight_func  = self.event_weight
        systematics  = self.systematics

        for row in self.tree:
            # Apply basic preselection
            if not preselection(row):
                continue

            constant_id_map = dict( [ (name, fcn(row) )  for name, fcn in id_functions.iteritems()] )
            for systematic in systematics:
                self.currect_systematic = systematic
                sys_id_map = dict( [ (name, fcn(row) )  for name, fcn in id_functions_with_sys.iteritems()] )
                sys_id_map.update(constant_id_map)
                # Get the generic event weight
                event_weight = weight_func(row)

                # Figure out which folder/region we are in, multiple regions allowed
                for folder, selection in cut_region_map.iteritems():
                    if folder.startswith(systematic): #Folder name starts with the systematic (if there are any)
                        if all( [ (sys_id_map[name] == info ) for name, info in selection.iteritems() ] ):
                            #all cuts match the one of the region, None means the cut is not needed
                            fill_histos(histos, folder, row, event_weight)

    def finish(self):
        self.write_histos()

if __name__ == "__main__":
    import pprint
    pprint.pprint(TauEffBase.build_folder_structure())
