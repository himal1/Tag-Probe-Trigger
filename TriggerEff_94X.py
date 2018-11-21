import FWCore.ParameterSet.Config as cms
import os
import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario

### USAGE:
###    cmsRun TriggerEff.py <scenario> [ <id> [ <binning1> ... <binningN> ] ]
###    ex> cmsRun TriggerEff.py mc IsoMu20_from_Tight2012
### scenarios:
###   - data_25ns (default)  
###   - mc_weight

PassingProbe = ""
ProbeCondition = ""
if "_from_" in args[1]:
    PassingProbe = args[1].split("_from_")[0]
    ProbeCondition = args[1].split("_from_")[1]
else:
    PassingProbe = args[1]
    ProbeCondition = "None"

print "PassingProbe: " + PassingProbe
print "On top of " + ProbeCondition

def Add_AdditionalBranches( _PassingProbe, _ProbeCondition):
    # -- common probe conditions --  delta R #
    
    Template.Variables.pair_deltaR = cms.vstring("pair_deltaR", "0", "999", "")
    Template.Variables.tkTrackerLay = cms.vstring("tkTrackerLay", "-1" , "999" , "")
    Template.Variables.tkPixelLay =  cms.vstring("tkPixelLay","-1","999" , "")
    Template.Variables.dzPV = cms.vstring("dzPV","-1000","1000", "")
    Template.Variables.dB = cms.vstring("dB","-1000","1000", "")


    # -- collect all variables used for efficiency defintion (passing probe, probe) -- #
    List_VarName = []
    if "_and_" in PassingProbe:
        for varName in PassingProbe.split("_and_"): List_VarName.append( varName )
    else:
        List_VarName.append( PassingProbe )

    if "_and_" in ProbeCondition:
        for varName in ProbeCondition.split("_and_"): List_VarName.append( varName )
    elif "None" != ProbeCondition:
        List_VarName.append( ProbeCondition )

    # -- add them in Categories if Categories didn't have them -- #
    # -- so, all variables in efficiency definition should be 0 or 1, not "Variables"! -- #
    for varName in List_VarName:
        if getattr(Template.Categories, varName, "NotExist") == "NotExist": # -- if "Template" didn't have the variable -- #
            setattr(Template.Categories, varName, cms.vstring(varName, "dummy[pass=1,fail=0]") )


def Determine_Shape( _PassingProbe, _ProbeCondition ):
    _shape = cms.vstring( "CBPlusExpo" ) # -- default -- #
    return _shape

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-muon Mass", "2.8", "3.35", "GeV/c^{2}"),
        pt = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        phi = cms.vstring("muon #phi", "-3.14", "3.14", ""),
        tag_pt = cms.vstring("tag muon p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("muon #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        tag_abseta = cms.vstring("tag muon |#eta|", "0", "2.5", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        # combRelIsoPF04dBeta = cms.vstring("pf relative isolation", "0", "999", ""),
        # tkIso = cms.vstring("Tracker Isolation", "0", "999", ""),
        # relTkIso = cms.vstring("Relative Tracker Isolation", "0", "999", ""),
    ),

    Categories = cms.PSet(
        # Tight2012 = cms.vstring("Tight2012", "dummy[pass=1,fail=0]"),
        # tag_IsoMu24 = cms.vstring("tag_IsoMu24", "dummy[pass=1,fail=0]"),
        tag_Mu7p5_L2Mu2_PathMacHIA = cms.vstring("tag_Mu7p5_L2Mu2_PathMacHIA", "dummy[pass=1,fail=0]"),
        Mu7p5_Track2_Jpsi_TK = cms.vstring("Mu7p5_Track2_Jpsi_TK", "dummy[pass=1,fail=0]"),
        tag_Mu7p5_Track2_Jpsi_MU =cms.vstring("tag_Mu7p5_Track2_Jpsi_MU", "dummy[pass=1,fail=0]"),
        TMOST = cms.vstring("TMOST", "dummy[pass=1,fail=0]"),
        Dimuon0_Jpsi25_L2MacHA = cms.vstring("Dimuon0_Jpsi25_L2MacHA","dummy[pass=1,fail=0]")
        #tkTrackerLay = cms.vstring("tkTrackerLay","-1","999","")
    ),

    Expressions = cms.PSet(),
    Cuts = cms.PSet(),

    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        CBPlusExpo = cms.vstring(
            "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.06], alpha[3., 0.5, 5.], n[1, 0.1, 100.])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),

        # -- PDF Sets for the case of failing on error calculation by MINOS -- #

        vpvPlusExpo5 = cms.vstring(
            "Voigtian::signalPass1(mass, meanPass1[91.1,86.3,96.1], width[2.495], sigmaPass1[2.52,1.4,5.1])",
            "Voigtian::signalPass2(mass, meanPass2[91.2,78.4,104.5], width,        sigmaPass2[5.12,1.3,8.2])",
            "SUM::signalPass(vFracPass[0.7,0,1]*signalPass1, signalPass2)",
            "Voigtian::signalFail1(mass, meanFail1[91.3,86.1,96.3], width[2.495], sigmaFail1[2.4,1,5.3])",
            "Voigtian::signalFail2(mass, meanFail2[91.2,82.5,100.5], width,        sigmaFail2[5.2,1,10.1])",
            "SUM::signalFail(vFracFail[0.7,0,1]*signalFail1, signalFail2)",
            "Exponential::backgroundPass(mass, lp[-0.11,-1.3,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.12,-1.2,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        
        vpvPlusExpoMin70 = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,3,10])",
            "SUM::signal(vFrac[0.8,0.5,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0.7,1]",
            "signalFractionInPassing[0.9]"
        ),

        ###############################################
        # -- for systematic uncertainty estimation -- #
        ###############################################
        # -- signal shape variation -- #
        cpvPlusCMS = cms.vstring(
            "CBShape::signal1(mass, m01[90,89,92], sigma1[2,1,3], alpha1[1,0,5], n1[1,0,10])",
            "Voigtian::signal2(mass, mean2[90,89,92], width2[2.495], sigma2[2,1,3])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "RooCMSShape::backgroundPass(mass, alphaPass[70.,60.,90.], betaPass[0.02, 0.01,0.1], gammaPass[0.001, 0.,0.1], peakPass[90.0])",
            "RooCMSShape::backgroundFail(mass, alphaFail[70.,60.,90.], betaFail[0.02, 0.01,0.1], gammaFail[0.001, 0.,0.1], peakPass)",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),

        cpvPlusExpo2 = cms.vstring(
            "CBShape::signal1(mass, m01[90,89,92], sigma1[2,1,3], alpha1[1,0,5], n1[1,0,10])",
            "Voigtian::signal2(mass, mean2[90,89,92], width2[2.495], sigma2[2,1,3])",
            "SUM::signalPass(vFrac[0.8,0,1]*signal1, signal2)",
            "CBShape::signal3(mass, m03[90,89,92], sigma3[2,1,3], alpha3[1,0,5], n3[1,0,10])",
            "Voigtian::signal4(mass, mean4[90,89,92], width4[2.495], sigma4[2,1,3])",
            "SUM::signalFail(vFrac[0.8,0,1]*signal3, signal4)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        # -- background shape variation -- #
        vpvPlusQuad = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Chebychev::backgroundPass(mass, {cPass1[0,-1,1], cPass2[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {cFail1[0,-1,1], cFail2[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),


    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(False),

    Efficiencies = cms.PSet(), # will be filled later
)

Add_AdditionalBranches( PassingProbe, ProbeCondition)

print "\nTemplate.Variables: \n", Template.Variables
print "\nTemplate.Categories: \n", Template.Categories

# print Template.Categories

PtMin = 3
EtaMax = 2.4
PT_ETA_BINS = cms.PSet(
    pt = cms.vdouble( 0, 9999 ), # -- Will be set later -- #
    abseta = cms.vdouble( 0.0, 0.9, 1.2, 2.1, 2.4 ),
)

# -- refer to the binning used in ID+Iso case -- #
# -- https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonTagAndProbeTreesRun2#Scale_factors_review -- #
PT_ETA_BINS.pt = cms.vdouble( 3., 3.25, 3.5, 4., 4.5, 5., 7., 9., 11., 15., 20., 25., 30., 35., 40., 50., 200.)
PT_BINS = cms.PSet(
    pt = cms.vdouble( 0, 9999 ), #Will be set later
    abseta = cms.vdouble(0.0, EtaMax),
)

PT_BINS.pt = cms.vdouble( 2,4,6,8,10,14, 18, 22, 24, 26, 30, 40, 50, 60, 120, 200)

ETA_BINS = cms.PSet(
	eta = cms.vdouble(-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4),
    pt = cms.vdouble( PtMin, 9999 ),
)
phi = 3.141592
degree15 = phi / 12;
PHI_BINS = cms.PSet(
    phi = cms.vdouble( (-1)*degree15*12, (-1)*degree15*11, (-1)*degree15*9, (-1)*degree15*7, (-1)*degree15*5, (-1)*degree15*3, (-1)*degree15*1, degree15*1, degree15*3, degree15*5, degree15*7, degree15*9, degree15*11, degree15*12),
                    
    pt = cms.vdouble( PtMin, 9999 ),
    abseta = cms.vdouble(0.0, EtaMax),
)

VTX_BINS  = cms.PSet(
	tag_nVertices = cms.vdouble(
        0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 
        20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5, 36.5, 38.5, 
        40.5, 42.5, 44.5, 46.5, 48.5, 50.5),
    pt = cms.vdouble(  PtMin, 9999 ),
    abseta = cms.vdouble(  0.0, EtaMax),
)

SINGLE_BIN  = cms.PSet(
    pt = cms.vdouble( PtMin, 9999 ),
    abseta = cms.vdouble( 0.0, EtaMax),
)


process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring("JPsiMuMu_MC_JPsi_TnP_Trigger_Tree_FinalCombined1.root",
                                 "JPsiMuMu_MC_JPsi_TnP_Trigger_Tree_FinalCombined2.root",
                                 "JPsiMuMu_MC_JPsi_TnP_Trigger_Tree_FinalCombined3.root",
                                 #"Run2017C_17Nov2017_TnP_Trigger_Tree_CombinedFinal.root"
                                 ),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_MuonTrigger_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

#Add the variables for PU reweighting
if "_weight" in scenario:
    process.TnP_MuonID.WeightVariable = cms.string("weight")
    process.TnP_MuonID.Variables.weight = cms.vstring("weight","-10","10","")

IDS = [args[1]] #here the id is taken from the arguments provided to cmsRun 
# ALLBINS = [ ("pt",PT_BINS), ("eta",ETA_BINS), ("phi",PHI_BINS), ("vtx",VTX_BINS), ("pteta",PT_ETA_BINS), ("single",SINGLE_BIN) ]
# ALLBINS = [ ("pt",PT_BINS), ("pteta",PT_ETA_BINS) ]
# ALLBINS = [ ("pteta",PT_ETA_BINS) ]
ALLBINS = [ ("pt",PT_ETA_BINS)]

if len(args) > 1 and args[1] not in IDS: IDS += [ args[1] ]
for ID in IDS:
    print "now doing ",ID
    if len(args) > 1 and ID != args[1]: continue
    for X,B in ALLBINS:
        if len(args) > 2 and X not in args[2:]: continue
        #Add the information about ID and binning into the outputFileName
        module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_MuonTrigger_%s_%s_%s.root" % (scenario, ID, X)))
        
        shape = Determine_Shape( PassingProbe, ProbeCondition)

        #DEN: Binning
        DEN = B.clone(); num = ID; numstate = "pass"

        # -- tag condition for all binning -- //
        DEN.tag_Mu7p5_Track2_Jpsi_MU = cms.vstring("pass")
        #DEN.tag_Mu7p5_L2Mu2_PathMacHIA=cms.vstring("pass")
        DEN.tag_pt = cms.vdouble(7,9999)
        DEN.tag_abseta = cms.vdouble(0, 2.4)
        
        # condition between tag & probe muons for all binning  -- #
        DEN.pair_deltaR = cms.vdouble(0.5, 999)#delta R cut larger than 0.5 acc to CMS AN -2015/297
        DEN.TMOST = cms.vstring("pass")
        DEN.tkTrackerLay = cms.vdouble(5, 999)
        DEN.tkPixelLay =  cms.vdouble(0,999)
        DEN.dzPV = cms.vdouble(-20., 20.)
        DEN.dB = cms.vdouble(-0.3,0.3)
        DEN.Mu7p5_Track2_Jpsi_TK = cms.vstring("pass")
        #DEN.Dimuon0_Jpsi25_L2MacHA = cms.vstring("pass")

        
        
        if "_from_" in ID:
            parts = ID.split("_from_")
            num = parts[0]
            # add Additional ID conditions to the binning ... 
            # ex> cmsRun TriggerEff.py mc IsoMu20_and_Tight2012_from_SIP4_and_PFIso25 => SIP4 and PFIso25 info. is added to the binning definition
            for D in parts[1].split("_and_"):
                if D == "SIP4":      DEN.SIP = cms.vdouble(0,4.0)
                elif D == "PFIso25": DEN.pfCombRelIso04EACorr  = cms.vdouble(0,0.25)
                elif D == "PFIso40": DEN.pfCombRelIso04EACorr  = cms.vdouble(0,0.40)
                #elif D == "TMOST": DEN.TMOST = cms.vstring("pass")
                #elif D == "tkTrackerLay": DEN.tkTrackerLay = cms.vdouble(5,9999)
                #elif D == "tkPixelLay" : DEN.tkPixelLay =  cms.vdouble(0,9999)
                #elif D == "dzPV" :   DEN.dzPV = cms.vdouble(-20.,20.)
                #elif D == "dB" :  DEN.dB = cms.double(-0.3,0.3)

                # elif D == "dBeta": DEN.combRelIsoPF04dBeta = cms.vdouble(0, 0.12)
                # elif D == "dBeta_015": DEN.combRelIsoPF04dBeta = cms.vdouble(0, 0.15)
                # elif D == "dBeta_025": DEN.combRelIsoPF04dBeta = cms.vdouble(0, 0.25)
            	# elif D == "RelTrkIso_010": DEN.relTkIso = cms.vdouble(0, 0.10)

                # Set D as the variable of DEN ... DEN.D = cms.vstring("pass")
                else: setattr(DEN, D, cms.vstring("pass"))

        print "#" * 100
        print "Binning variable: ", X
        print "Binning: ", DEN
        print "PDF: ", shape

        # numString: EfficiencyCategoryState variable. 
        # ex> cmsRun TriggerEff.py mc IsoMu20_and_Tight2012_from_SIP4_and_PFIso25 => numString = cms.vstring("IsoMu20", "pass", "Tight2012", "pass")
        numString = cms.vstring()
        for N in num.split("_and_"):
            numString += [N, "pass"]
                
        print "Passing probe condition: ", numString
        print "#" * 100
        print "\n"
        
        #Set Efficiency
        setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
            EfficiencyCategoryAndState = numString,
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = DEN,
            BinToPDFmap = cms.vstring(shape)
        ))

        #Add mcTure Efficinecy when mc fitting
        if scenario.find("mc") != -1:
            # setattr(module.Efficiencies, ID+"_"+X+"_mcTrue", cms.PSet(
            #     EfficiencyCategoryAndState = numString,
            #     UnbinnedVariables = cms.vstring("mass"),
            #     BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
            # ))
            #When mc is PU-weighted, "weight" variable is added to UnbinnedVariables
            if "_weight" in scenario:
                getattr(module.Efficiencies, ID+"_"+X          ).UnbinnedVariables.append("weight")
                # getattr(module.Efficiencies, ID+"_"+X+"_mcTrue").UnbinnedVariables.append("weight")

        #Add module to the process
        setattr(process, "TnP_MuonID_"+ID+"_"+X, module)        
        #Add a path of module to the process
        setattr(process, "run_"+ID+"_"+X, cms.Path(module))
