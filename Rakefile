# Get common recipes
recipes = ENV['CMSSW_BASE'] + '/src/FinalStateAnalysis/PlotTools/rake/recipes.rake'
import recipes
require ENV['CMSSW_BASE'] + '/src/FinalStateAnalysis/PlotTools/rake/tools.rb'

$jobid = ENV['jobid']

# Figure out what run period we are in
$period = '7TeV'
if $jobid.include? '8TeV'
  $period = '8TeV'
end

################################################################################
## Sample names ################################################################
################################################################################
#
# Get sample names containing a substring
def get_sample_names(substring)
  inputs = Dir.glob("inputs/#{$jobid}/*.txt").select {|x| x.include? substring}
  inputs = inputs.map{|x| File.basename(x).sub(".txt", "")}
  return inputs
end

#puts get_sample_names('data_DoubleMu')

samples = Hash[
  "ewk" => Array['WplusJets_madgraph', 'TTplusJets_madgraph', 'Zjets_M50'],
  "wjets" => Array['WplusJets_madgraph'],
  "zjets_clone" => Array['Zjets_ZToMuMu_M50'],
  "zjets_spinned" => Array['Zjets_TauSpinned_M50'],               
  "zjets" => get_sample_names('Zjets'),
  "diboson" => get_sample_names('ZZ')+get_sample_names('WZ') + get_sample_names('WW'),  # Automagically figure out what data samples we have
  "data_m" =>  get_sample_names("data_SingleMu"),
  "data_mm" =>  get_sample_names("data_DoubleMu"),
]

# Function to get the .root files for an analyzer and samples
def get_analyzer_results(analyzer, the_samples)
  output = Array.new
  analyzer_base = analyzer.sub('.py', '')
  the_samples.each do |sample|
    output << "results/#{$jobid}/#{analyzer_base}/#{sample}.root"
  end
  return output
end


################################################################################
## Recipes to analyze TauIDEff
##  targets:
##     mt
##     mm
################################################################################

file "TauEffZMM.py" => ["TauEffBase.py"] do |t|
  sh "touch #{t.name}"
end
file "TauEffZMT.py" => ["TauEffBase.py"] do |t|
  sh "touch #{t.name}"
end



mm_results = get_analyzer_results("TauEffZMM.py", samples['ewk'] + samples['data_m'] + samples['diboson'])
mt_results = get_analyzer_results("TauEffZMT.py", samples['ewk'] + samples['data_m'] + samples['diboson'] + samples["zjets_clone"] + samples["zjets_spinned"])


task :mm => mm_results
task :mt => mt_results

