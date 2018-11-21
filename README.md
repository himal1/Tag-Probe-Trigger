# Tag-Probe-Trigger
#TnP analysis of 2017 Dataset
#commands
#place for production root file for TnP
/uscms_data/d3/hacharya/TnP_Himal_Stefan/CMSSW_9_4_0/src/
#place for doing fitter tag and probe 
/uscms_data/d3/hacharya/TnP_Himal_Stefan/CMSSW_9_4_0/src/FitterTree/NEW_FItterProgram/FinalWork/VeryFinalCalculation
#The MC DAS link is
https://cmsweb.cern.ch/das/request?instance=prod/global&input=file+dataset%3D%2FJpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8%2FRunIIFall17DRPremix-RECOSIMstep_9\
4X_mc2017_realistic_v10-v1%2FAODSIM


#MC data set name is
file dataset=/JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-RECOSIMstep_94X_mc2017_realistic_v10-v1/AODSIM


#skimtree commands

 ./skimTree JPsiMuMu_MC_JPsi_TnP_Trigger_Tree_Combined.root JPsiMuMu_MC_JPsi_TnP_Trigger_Tree_CombinedSkimmed.root -r "all" -k "TMOST tkTrackerLay tkPixelLa\
y dzPV dB  Mu7p5_Track2_PathMacHA Mu7p5_L2Mu2_PathMacHA  Mu7p5_L2Mu2_Leg2MacHA Mu7p5_Track2_Leg2MacHA Mu7p5_L2Mu2_Leg3MacHA Mu7p5_Track2_Leg3MacHA Dimuon0_J\
psi_L2MacHA Dimuon0_Jpsi3p5_Muon2_L2MacHA  tag_Mu7p5_Track2_PathMacHIA tag_Mu7p5_L2Mu2_PathMacHIA  tag_Mu7p5_L2Mu2_Leg2MacHIA tag_Mu7p5_Track2_Leg2MacHIA ta\
g_Mu7p5_L2Mu2_Leg3MacHIA tag_Mu7p5_Track2_Leg3MacHIA tag_Dimuon0_Jpsi_L2MacHIA tag_Dimuon0_Jpsi3p5_Muon2_L2MacHIA" -c "tag_ Mu7p5_Track2_Leg2MacHIA==1 && Mu\
7p5_Track2_Leg2MacHA==1 && pair_probeMultiplicity == 1 && TMOST==1 && tkTrackerLay>5 && tkPixelLay>0 && abs(dzPV)<20 && abs(dB)< 0.3"

#TnP fitter commands 
#L1L2 efficiency the command line to run is
cmsRun TriggerEff_94X.py data Dimuon0_Jpsi25_L2MacHA

#For L3

cmsRun TriggerEff_94X.py data Mu_L3
#(the tag and probe cut should be adjusted in fitter)
