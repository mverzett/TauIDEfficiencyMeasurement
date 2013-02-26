# Get common recipes
recipes = ENV['CMSSW_BASE'] + '/src/FinalStateAnalysis/PlotTools/rake/recipes.rake'
import recipes

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
  "ewk" => Array['Zjets_M50', 'WplusJets_madgraph', 'TTplusJets_madgraph'],
  "wjets" => Array['WplusJets_madgraph'],
  "zjets" => Array['Zjets_M50'],
  # Automagically figure out what data samples we have
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




mm_results = get_analyzer_results("TauEffZMM.py", samples['ewk'] + samples['data_mm'] + samples['data_m'])
mt_results = get_analyzer_results("TauEffZMT.py", samples['ewk'] + samples['data_m'] )


task :mm => mm_results
task :mt => mt_results

################################################################################
# Recipes to make data cards (plots come for free)
#  targets:
#      mmt_shapes
#      emt_shapes
#      eet_shapes
#      cards
#      copycards -> move cards to official HTT CVS area
################################################################################


# Output directories
$mmtdir = "results/#{$jobid}/plots/mmt/"
$emtdir = "results/#{$jobid}/plots/emt/"
$eetdir = "results/#{$jobid}/plots/eet/"

directory $mmtdir
directory $eetdir
directory $emtdir

# Recipes for adding stat. error shapes.  Makes a new file task:
# input_file_stat_errors.root => input_file.root
def add_fake_errors(input_file, prefix)
  output_file = input_file.sub('.root', '_statshapes.root')
  output_sys_list = input_file.sub('.root', '_statshapes.txt')
  file output_file => [input_file] do |t|
    sh "add_stat_shapes.py #{input_file} #{output_file} --filter '#{prefix}/fakes' --prefix CMS_vhtt_#{$period}_#{prefix}_fakeshape > #{output_sys_list}"
  end
end

# The .root files with the shapes
file "#{$mmtdir}/mmt_shapes_#{$period}.root" => ['WHPlotterMMT.py', 'WHPlotterBase.py'] do |t|
  sh "python WHPlotterMMT.py"
end
add_fake_errors("#{$mmtdir}/mmt_shapes_#{$period}.root", "mmt")
task :mmt_shapes => ["#{$mmtdir}/mmt_shapes_#{$period}_statshapes.root" ]

file "#{$emtdir}/emt_shapes_#{$period}.root" => ['WHPlotterEMT.py', 'WHPlotterBase.py'] do |t|
  sh "python WHPlotterEMT.py"
end
add_fake_errors("#{$emtdir}/emt_shapes_#{$period}.root", "emt")
task :emt_shapes => ["#{$emtdir}/emt_shapes_#{$period}_statshapes.root" ]

file "#{$eetdir}/eet_shapes_#{$period}.root" => ['WHPlotterEET.py', 'WHPlotterBase.py'] do |t|
  sh "python WHPlotterEET.py"
end
add_fake_errors("#{$eetdir}/eet_shapes_#{$period}.root", "eet")
task :eet_shapes => ["#{$eetdir}/eet_shapes_#{$period}_statshapes.root" ]

$carddir = "results/#{$jobid}/cards"

# Combine all category shape files into the datacard project
file "#{$carddir}/shapes_unmorphed.root" => [
  "#{$mmtdir}/mmt_shapes_#{$period}_statshapes.root", 
  "#{$emtdir}/emt_shapes_#{$period}_statshapes.root"] do |t|
  sh "mkdir -p #{$carddir}"
  sh "hadd -f #{t.name} #{t.prerequisites.join(' ')}"
end

# Masses where actually have signal samples
pivot_masses = "120,130,140"

# Now do the horizontal morphing
file "#{$carddir}/shapes.root" => ["#{$carddir}/shapes_unmorphed.root"] do |t|
  puts "Doing horizontal morphing"
  sh "mkdir -p #{$carddir}"
  sh "cp #{t.prerequisites[0]} #{t.name}"
  sh "horizontal-morphing.py --categories='emt,mmt' --samples='WH{MASS},WH_hww{MASS}' --uncerts='' --masses='110,120,130,140' --step-size=1 -i #{t.name}"
end

stat_shape_lists = Dir.glob("results/#{$jobid}/plots/*/*_statshapes.txt")
# We make these the dependency, since they are always rpoduced w/ the .txt lists
stat_shapes = Dir.glob("results/#{$jobid}/plots/*/*_statshapes.root")

# We need to build the unc.conf and unc.vals with our stat shapes
file "#{$carddir}/unc.conf" => ["card_config/unc.conf.base.#{$period}"] + stat_shapes do |t|
  # Copy the basic template
  sh "mkdir -p #{$carddir}"
  sh "cp #{t.prerequisites[0]} #{t.name}"
  # Append all the stat shape types
  sh "echo '' >> #{t.name}"
  sh "echo '# Stat shape uncertainties' >> #{t.name}"
  stat_shape_lists.each do |list|
    sh "cat #{list} | xargs -n 1 -I {} echo '{} shape' >> #{t.name}"
  end
end

file "#{$carddir}/unc.vals" => ["card_config/unc.vals.base.#{$period}", 
  "#{$carddir}/shapes.root"] + stat_shapes do |t|
  # Copy the basic template
  sh "mkdir -p #{$carddir}"
  puts t.investigation
  sh "cp #{t.prerequisites[0]} #{t.name}"
  # Append all the stat shape types
  sh "echo '' >> #{t.name}"
  sh "echo '# Stat shape uncertainties' >> #{t.name}"
  puts stat_shape_lists.join(", ")
  stat_shape_lists.each do |list|
    sh "cat #{list} | grep mmt | xargs -n 1 -I {} echo 'mmt fakes {} 1.0' >> #{t.name}"
    sh "cat #{list} | grep emt | xargs -n 1 -I {} echo 'emt fakes {} 1.0' >> #{t.name}"
  end
  sh "python get_fake_systematic.py #{$carddir}/shapes.root mmt CMS_vhtt_mmt_fakes_#{$period} >> #{t.name}"
  sh "python get_fake_systematic.py #{$carddir}/shapes.root emt CMS_vhtt_emt_fakes_#{$period} >> #{t.name}"
  # For testing a full correlated error.
  #sh "echo 'emt fakes CMS_vhtt_llt_fakes 1.3' >> #{t.name}"
  #sh "echo 'mmt fakes CMS_vhtt_llt_fakes 1.3' >> #{t.name}"
end

file "#{$carddir}/cgs.conf" => ["card_config/cgs.conf.#{$period}"] do |t|
  sh "mkdir -p #{$carddir}"
  sh "cp #{t.prerequisites[0]} #{t.name}"
end

def make_datacard_task(mass, channel, categories)
  card = "#{$carddir}/#{mass}/vhtt_#{channel}.txt"  
  file card => ["#{$carddir}/shapes.root", "#{$carddir}/unc.conf", "#{$carddir}/unc.vals", "#{$carddir}/cgs.conf"] do |t|
    sh "mkdir -p #{$carddir}"
    chdir($carddir) do
      sh "mkdir -p #{channel}/#{mass}"
      # Link the shape file in the card directory so combine can run from anywhere
      chdir("#{channel}/#{mass}") do
        sh "rm -f shapes.root"
        sh "ln -s ../../shapes.root"
      end
      if categories == ''
        sh "create-datacard.py -i shapes.root -o #{channel}/#{mass}/vhtt_#{channel}.txt #{mass}"
      else
        sh "create-datacard.py -i shapes.root -o #{channel}/#{mass}/vhtt_#{channel}.txt #{mass} --categories '#{categories}'"
      end
    end
  end
  return card
end

task :cards => []

cardmasses = Array[110, 120, 125, 130, 135, 140]

cardmasses.each do |mass|
  task :cards => make_datacard_task(mass, '2lt', '')
  task :cards => make_datacard_task(mass, 'emt', 'emt')
  task :cards => make_datacard_task(mass, 'mmt', 'mmt')
end

################################################################################
### Copying card configuration to official place ###############################
################################################################################

#$httcombodir="/afs/hep.wisc.edu/cms/efriis/HIG-12-051/src/HiggsAnalysis/HiggsToTauTau/setup/vhtt"
$httcombodir="/afs/cern.ch/work/f/friis/HttLimits/src/HiggsAnalysis/HiggsToTauTau/setup/vhtt"

file "#{$httcombodir}/cgs-sm-#{$period}-00.conf" => ["#{$carddir}/cgs.conf"] do |t|
    puts t.investigation
  sh "cp -v #{t.prerequisites[0]} #{t.name}"
end

file "#{$httcombodir}/unc-sm-#{$period}-00.conf" => ["#{$carddir}/unc.conf"] do |t|
  sh "cp -v #{t.prerequisites[0]} #{t.name}"
end

file "#{$httcombodir}/unc-sm-#{$period}-00.vals" => ["#{$carddir}/unc.vals"] do |t|
  sh "cp -v #{t.prerequisites[0]} #{t.name}"
end

file "#{$httcombodir}/vhtt_llt.inputs-sm-#{$period}.root" => ["#{$carddir}/shapes.root"] do |t|
  sh "cp -v #{t.prerequisites[0]} #{t.name}"
end

task :copycards => [
  "#{$httcombodir}/cgs-sm-#{$period}-00.conf",
  "#{$httcombodir}/unc-sm-#{$period}-00.conf",
  "#{$httcombodir}/unc-sm-#{$period}-00.vals",
  "#{$httcombodir}/vhtt_llt.inputs-sm-#{$period}.root"]  do |t|
    puts t.investigation
  end
