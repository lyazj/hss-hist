#!/usr/bin/env python3

import numpy as np

xs = {  # https://github.com/StephenChao/boostedHWW/blob/01dd15d5d6192a648c24a47a2710706a9d6a4033/preprocessing/ntuple_to_tree/Scripts/TransMergedMC.py#L47-L68
    "QCD_HT500to700": 30310000,
    "QCD_HT700to1000": 6444000,
    "QCD_HT1000to1500": 1092000,
    "QCD_HT1500to2000": 99760,
    "QCD_HT2000toInf": 20350,
    "ST_s-channel_4f_hadronicDecays": 11240,
    "ST_t-channel_antitop": 80950,
    "ST_t-channel_top": 136020,
    "ST_tW_antitop": 35850,
    "ST_tW_top": 35850,
    "TTToHadronic": 380094,
    "TTToSemiLeptonic": 364350.8,
    "WJetsToQQ_HT-400to600": 276500,
    "WJetsToQQ_HT-600to800": 59250,
    "WJetsToQQ_HT-800toInf": 28750,
    "WW_TuneCP5": 76250,
    "WZ_TuneCP5": 27550,
    "ZZ_TuneCP5_13TeV-pythia8": 12230,
    "ZJetsToQQ_HT-400to600": 114200,
    "ZJetsToQQ_HT-600to800": 25340,
    "ZJetsToQQ_HT-800toInf": 13100,
}

xs['QCD'] = sum([xs[key] for key in xs if 'QCD' in key])
xs['WJets'] = sum([xs[key] for key in xs if 'WJets' in key])
xs['Wcs'] = xs['WJets'] * 0.6741 * 0.46  # https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S043.8

weights = np.array([xs['QCD'], xs['Wcs']])
weights = weights / np.sum(weights)
print(*['%.6f' % weight for weight in weights])
