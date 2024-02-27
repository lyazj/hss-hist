#!/usr/bin/env python3

import ROOT
import uproot
import glob
import numpy as np
import awkward as ak
import mplhep as hep
import matplotlib.pyplot as plt

rootfiles = [file + ':Events' for file in sorted(glob.glob('samples/*/*/*.root'))][:10]
expressions = [
    'GenPart_pt',
    'GenPart_eta',
    'GenPart_phi',
    'GenPart_mass',
    'GenPart_genPartIdxMother',
    'GenPart_pdgId',
    'AK15Puppi_pt',
    'AK15Puppi_eta',
    'AK15Puppi_phi',
    'AK15Puppi_mass',
    'AK15Puppi_mass * (1 - AK15Puppi_rawFactor) * AK15Puppi_inclParTMDV2_resonanceMassCorr',
    'AK15Puppi_mass * (1 - AK15Puppi_rawFactor) * AK15Puppi_inclParTMDV2_visiableMassCorr',
    'FatJet_pt',
    'FatJet_eta',
    'FatJet_phi',
    'FatJet_mass',
    'FatJet_particleNet_mass',
]
bins = np.linspace(20, 170, 51)

# Match nearest Z within radius dr_max.
def match_z(pt, eta, phi, mass, pid, jet_pt, jet_eta, jet_phi, jet_mass, dr_max):
    p4 = ROOT.TLorentzVector(); p4.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass)
    iz, dr = None, None
    for i in range(len(pid)):
        if pid[i] != 23: continue  # PID_Z
        p4i = ROOT.TLorentzVector(); p4i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i])
        dri = p4i.DeltaR(p4)
        if dri > dr_max: continue
        if dr is None or dri < dr:
            iz, dr = i, dri
    return iz

# Match exact 2 quarks as daughters of mother im within radius dr_max.
def match_qq(pt, eta, phi, mass, mom, pid, im, dr_max):
    p4 = ROOT.TLorentzVector(); p4.SetPtEtaPhiM(pt[im], eta[im], phi[im], mass[im])
    iq = []
    for i in range(len(pid)):
        if pid[i] == 0 or abs(pid[i]) > 6: continue  # PID_Q
        if mom[i] != im: continue
        p4i = ROOT.TLorentzVector(); p4i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i])
        dri = p4i.DeltaR(p4)
        if dri > dr_max: continue
        iq.append(i)
    if len(iq) != 2: return None
    return iq

# Match nearest Z, then with exact 2 daughter quarks, both within radius dr_max.
def match_zqq(pt, eta, phi, mass, mom, pid, jet_pt, jet_eta, jet_phi, jet_mass, dr_max):
    iz = match_z(pt, eta, phi, mass, pid, jet_pt, jet_eta, jet_phi, jet_mass, dr_max)
    if iz is None:
        #print('DEBUG: no matching Z')
        return None
    iq = match_qq(pt, eta, phi, mass, mom, pid, iz, dr_max)
    if iq is None:
        #print('DEBUG: no matching qq')
        return None
    return [iz, *iq]

ak15_mass_list = []
ak15_resmass_list = []
ak15_vismass_list = []
ak8_mass_list = []
ak8_resmass_list = []
ievent = 0
iak8 = 0
iak15 = 0
print('Loading events...')
for gp_pt, gp_eta, gp_phi, gp_mass, gp_mom, gp_pid, \
        ak15_pt, ak15_eta, ak15_phi, ak15_mass, ak15_resmass, ak15_vismass, \
        ak8_pt, ak8_eta, ak8_phi, ak8_mass, ak8_resmass in \
        zip(*uproot.concatenate(rootfiles, expressions, how=tuple)):
    if ievent == 0:
        print('Performing truth matching...')
    elif ievent % 1000 == 0:
        print('%d events processed, %d/%d AK15 jets, %d/%d AK8 jets.'
              % (ievent, len(ak15_mass_list), iak15, len(ak8_mass_list), iak8))
    ievent += 1
    for jet_pt, jet_eta, jet_phi, jet_mass, jet_resmass, jet_vismass in \
            zip(ak15_pt, ak15_eta, ak15_phi, ak15_mass, ak15_resmass, ak15_vismass):
        iak15 += 1
        if not match_zqq(gp_pt, gp_eta, gp_phi, gp_mass, gp_mom, gp_pid,
                         jet_pt, jet_eta, jet_phi, jet_mass, 1.5): continue
        ak15_mass_list.append(jet_mass)
        ak15_resmass_list.append(jet_resmass)
        ak15_vismass_list.append(jet_vismass)
    for jet_pt, jet_eta, jet_phi, jet_mass, jet_resmass in \
            zip(ak8_pt, ak8_eta, ak8_phi, ak8_mass, ak8_resmass):
        iak8 += 1
        if not match_zqq(gp_pt, gp_eta, gp_phi, gp_mass, gp_mom, gp_pid,
                         jet_pt, jet_eta, jet_phi, jet_mass, 0.8): continue
        ak8_mass_list.append(jet_mass)
        ak8_resmass_list.append(jet_resmass)
print('Truth matching done!')

ak15_mass, ak15_resmass, ak15_vismass, ak8_mass, ak8_resmass = map(
    lambda x: np.histogram(x, bins=bins, density=True),
    (ak15_mass_list, ak15_resmass_list, ak15_vismass_list, ak8_mass_list, ak8_resmass_list)
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
