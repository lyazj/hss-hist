#!/usr/bin/env python3

import ROOT
import uproot
import glob
import numpy as np
import awkward as ak
import mplhep as hep
import matplotlib.pyplot as plt

plt.style.use(hep.style.CMS)

rootfiles = []
for rootdir in glob.glob('samples/*/*'):
    print(rootdir)
    rootfiles.extend([file + ':Events' for file in glob.glob(rootdir + '/*.root')])  # Can be truncated here.
print(*rootfiles, sep='\n')

expressions = [
    'FatJet_inclParTMDV2_probHbc',
    'FatJet_inclParTMDV2_probHbs + FatJet_inclParTMDV2_probHcc + FatJet_inclParTMDV2_probHcs + FatJet_inclParTMDV2_probHee + FatJet_inclParTMDV2_probHgg + FatJet_inclParTMDV2_probHmm + FatJet_inclParTMDV2_probHqq + FatJet_inclParTMDV2_probHss + FatJet_inclParTMDV2_probHtauhtaue + FatJet_inclParTMDV2_probHtauhtauh + FatJet_inclParTMDV2_probHtauhtaum',
    'FatJet_inclParTMDV2_probQCDb + FatJet_inclParTMDV2_probQCDbb + FatJet_inclParTMDV2_probQCDc + FatJet_inclParTMDV2_probQCDcc + FatJet_inclParTMDV2_probQCDothers',
]
bins = np.linspace(0, 1, 51)

print('Loading events...')
ak8_probHbc, ak8_probHothers, ak8_probQCD = uproot.concatenate(rootfiles, expressions, how=tuple)
print('%s events loaded.' % len(ak8_probHbc))
ak8_probHbc, ak8_probHothers, ak8_probQCD = map(
    lambda x: ak.flatten(x).to_numpy(),
    (ak8_probHbc, ak8_probHothers, ak8_probQCD)
)
print('%s AK8 jets.' % len(ak8_probHbc))
ak8_HbcVSQCD = 1 / (1 + ak8_probQCD / ak8_probHbc)

ak8_probHbc, ak8_probHothers, ak8_probQCD, ak8_HbcVSQCD = map(
    lambda x: np.histogram(x, bins=bins),
    (ak8_probHbc, ak8_probHothers, ak8_probQCD, ak8_HbcVSQCD)
)

plt.figure(figsize=(12, 9), dpi=150)
try:
    hep.cms.label(data=False, paper=False, supplementary=False, year=2018, lumi=59.7)
except Exception:
    hep.cms.label(data=False, label='Preliminary', year=2018, lumi=59.7)
hep.histplot(ak8_probHbc, histtype='step', label='probHbc')
hep.histplot(ak8_probHothers, linestyle='--', histtype='step', label='probHothers')
hep.histplot(ak8_probQCD, linestyle='--', histtype='step', label='probQCD')
hep.histplot(ak8_HbcVSQCD, histtype='step', label='HbcVSQCD')
plt.yscale('log')
plt.xlabel('ParT Score')
plt.ylabel('Unweighted Events')
plt.legend(loc='upper right')
plt.grid()
plt.tight_layout()
plt.savefig('plot.pdf')
