#!/usr/bin/env python3

import uproot
import glob
import numpy as np
import awkward as ak
import mplhep as hep
import matplotlib.pyplot as plt

rootfiles = [file + ':Events' for file in sorted(glob.glob('samples/*/*/*.root'))][:10]
expressions = [
    'AK15Puppi_mass',
    'AK15Puppi_mass * (1 - AK15Puppi_rawFactor) * AK15Puppi_inclParTMDV2_resonanceMassCorr',
    'AK15Puppi_mass * (1 - AK15Puppi_rawFactor) * AK15Puppi_inclParTMDV2_visiableMassCorr',
    'FatJet_mass',
    'FatJet_particleNet_mass',
]
bins = np.linspace(20, 170, 51)

ak15_mass, ak15_resmass, ak15_vismass, ak8_mass, ak8_resmass = map(
    lambda x: np.histogram(ak.flatten(x).to_numpy(), bins=bins, density=True),
    uproot.concatenate(rootfiles, expressions, how=tuple)
)
hep.histplot(ak15_mass, histtype='step', label='AK15 Mass')
hep.histplot(ak15_resmass, histtype='step', label='AK15 Res. Mass')
hep.histplot(ak15_vismass, histtype='step', label='AK15 Vis. Mass')
hep.histplot(ak8_mass, linestyle='--', histtype='step', label='AK8 Mass')
hep.histplot(ak8_resmass, linestyle='--', histtype='step', label='AK8 Res. Mass')
plt.xlabel('Mass [GeV]')
plt.ylabel('Density')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('plot.pdf')
