# b-taggingEfficienciesMeasurment
Measurment of the b-tagging efficiencies with 2TagCount method in tt dilepton events at 13 TeV using NanoAOD
FrameWork Bamboo: A high-level HEP analysis library for ROOT::RDataFrame
for more information, take a look here: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html
- The current code work with ; Fall17, Autmun 2018 and Summer19UL17 nanoAOD data. 

# Bamboo Installation(1st time):
Important notes:
- Only possible on CC7
- Follow the instructions here: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/install.html
- Use the LCG_96bpython3 environment
- You can ignore everything related to "SAMADhi"
- However I summed all what you need to start working with bamboo in the commands just below: 
```
mkdir bamboodev
cd bamboodev
# make a virtualenv
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
# clone and install bamboo
git clone -o upstream https://gitlab.cern.ch/cp3-cms/bamboo.git
pip install ./bamboo
# clone and install plotIt
git clone -o upstream https://github.com/cp3-llbb/plotIt.git
mkdir build-plotit
cd build-plotit
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
make -j2 install
cd -
```
# Environment setup (always *):
- In your ~/.bashrc add:
```
function cms_env() {
    module purge
    module load grid/grid_environment_sl6
    /cvmfs/cms.cern.ch/cmsset_default.sh
    module load crab/crab3
    module load slurm/slurm_utils
    module load cms/cmssw
}
alias bamboo_env="source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh"
alias bambooenv="source $HOME/bamboodev/bamboovenv/bin/activate"
```
- then every time you want to setup your bamboo enviroment:
```
cms_env
voms-proxy-init --voms cms
bamboo_env
bambooenv
```
# Update Bamboo :
```
cd bamboodev/bamboo
git checkout master
git pull upstream master
pip install --upgrade .
```
# How to run:
- I recommend to test locally first with --maxFiles=1, after you can submit to slurm with --distributed=driver

```
-s : --systematics add to your plots PSweight (FSR , ISR) , PDFs and six QCD scale variations, ele_id, ele_reco, pu, BtagWeight, DY, top reweighting, all what you have as sytematics in your .yml configuration
-v : --verbose     give you more print out for debugging 
-m : --module      your analysis script
-y : --yeilds      to get the latex table with yield at different stage of selection

bambooRun --distributed=driver -v -s -m TTbarAnalysis.py:TTbarDileptonMeasurment config/choose_One_.yml -o ~/path_to_your_Output_dir/
```
- In case you want to run plotIt again (after changing few options such fill color, legend position, unable systematics, etc...)
```
plotIt -i /path_to_your_dir/ -o /path_to_your_dir/plots_{add_era: 2016, 2017, 2018} -y -e era /In_your_Output_dir/plots.yml
```
- Or --onlypost
```
bambooRun --onlypost -v -s -m TTbarAnalysis.py:TTbarDileptonMeasurment config/choose_One_.yml -o ~/path_to_your_Output_dir/
```
- N.B:
- In case you have a small surprise from slurm : whatever error but the jobs are submitted, wait for them to finish and run with --distributed=finalize (afterwards, to merge the outputs)
- to re-submit only the failed jobs :
```
sbatch --array=_your_failed_jobs_seperated_by_comma --partition=cp3 --qos=cp3 --export=ALL --license=cms_storage:3 --exclude=mb-sky002 ~/path_to_your/batch/slurmSubmission.sh
```

# For NanoAOD the content of the branches 
- [Fall17 V2 94X data](https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/data94Xv2_doc.html)
- [Fall17 V2 94X MC ](https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc94Xv2_doc.html)

- [Summer19UL17 106X data](https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/data106Xul17_doc.html)
- [Summer19UL17 106X MC](https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc106Xul17_doc.html#MET)
