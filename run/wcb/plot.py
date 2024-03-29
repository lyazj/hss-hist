#!/usr/bin/env python3

import os
import ROOT
import uproot
import glob
import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt

plt.style.use(hep.style.CMS)
expressions = [
    'GenPart_pt',
    'GenPart_eta',
    'GenPart_phi',
    'GenPart_mass',
    'GenPart_genPartIdxMother',
    'GenPart_pdgId',
    'FatJet_pt',
    'FatJet_eta',
    'FatJet_phi',
    'FatJet_mass',
    'FatJet_inclParTMDV2_probHbc',
    'FatJet_inclParTMDV2_probHbs + FatJet_inclParTMDV2_probHcc + FatJet_inclParTMDV2_probHcs + FatJet_inclParTMDV2_probHee + FatJet_inclParTMDV2_probHgg + FatJet_inclParTMDV2_probHmm + FatJet_inclParTMDV2_probHqq + FatJet_inclParTMDV2_probHss + FatJet_inclParTMDV2_probHtauhtaue + FatJet_inclParTMDV2_probHtauhtauh + FatJet_inclParTMDV2_probHtauhtaum',
    'FatJet_inclParTMDV2_probQCDb + FatJet_inclParTMDV2_probQCDbb + FatJet_inclParTMDV2_probQCDc + FatJet_inclParTMDV2_probQCDcc + FatJet_inclParTMDV2_probQCDothers',
]
weights = {
    'samples/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/HIG-RunIISummer20UL18MiniAODv2-01120-wcb-merged': 2.1215 * 59700,
    'samples/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/HIG-RunIISummer20UL18MiniAODv2-00715-wcb-merged': 0.23338 * 59700,
    'samples/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/HIG-RunIISummer20UL18MiniAODv2-00714-wcb-merged': 0.04983 * 59700,
    'samples/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/HIG-RunIISummer20UL18MiniAODv2-00716-wcb-merged': 0.024404 * 59700,
}
bins = np.linspace(0, 1, 51)
hists = [np.histogram([], bins=bins)[0].astype('float') for i in range(4)]
weighted_hists = [np.histogram([], bins=bins)[0].astype('float') for i in range(4)]
weighted_hists_400toInf = [np.histogram([], bins=bins)[0].astype('float') for i in range(4)]

# Match nearest W within radius dr_max.
def match_w(pt, eta, phi, mass, pid, jet_pt, jet_eta, jet_phi, jet_mass, dr_max):
    p4 = ROOT.TLorentzVector(); p4.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass)
    iw, dr = None, None
    for i in range(len(pid)):
        if abs(pid[i]) != 24: continue  # PID_W
        p4i = ROOT.TLorentzVector(); p4i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i])
        dri = p4i.DeltaR(p4)
        if dri > dr_max: continue
        if dr is None or dri < dr:
            iw, dr = i, dri
    return iw

# Match cb quark pair as daughters of mother im within radius dr_max.
def match_cb(pt, eta, phi, mass, mom, pid, im, dr_max):
    p4 = ROOT.TLorentzVector(); p4.SetPtEtaPhiM(pt[im], eta[im], phi[im], mass[im])
    iq = []
    for i in range(len(pid)):
        if abs(pid[i]) not in [4, 5]: continue  # PID_C PID_B
        if mom[i] != im: continue
        p4i = ROOT.TLorentzVector(); p4i.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i])
        dri = p4i.DeltaR(p4)
        if dri > dr_max: continue
        iq.append(i)
    if len(iq) != 2: return None
    return iq

# Match nearest W, then with cb daughter quark pair, both within radius dr_max.
def match_wcb(pt, eta, phi, mass, mom, pid, jet_pt, jet_eta, jet_phi, jet_mass, dr_max):
    iw = match_w(pt, eta, phi, mass, pid, jet_pt, jet_eta, jet_phi, jet_mass, dr_max)
    if iw is None:
        #print('DEBUG: no matching W')
        return None
    iq = match_cb(pt, eta, phi, mass, mom, pid, iw, dr_max)
    if iq is None:
        #print('DEBUG: no matching cb')
        return None
    return [iw, *iq]

for rootdir in glob.glob('samples/*/*'):

    print(rootdir)
    rootfiles = [file + ':Events' for file in glob.glob(rootdir + '/*.root')]
    print(*rootfiles, sep='\n')

    print('Loading events...')
    events = uproot.concatenate(rootfiles, expressions, how=tuple)
    print('%s events loaded.' % len(events[0]))

    ak8_probHbc_list = []
    ak8_probHothers_list = []
    ak8_probQCD_list = []
    ievent = 0
    iak8 = 0

    print('Performing truth matching...')
    for gp_pt, gp_eta, gp_phi, gp_mass, gp_mom, gp_pid, \
        ak8_pt, ak8_eta, ak8_phi, ak8_mass, \
        ak8_probHbc, ak8_probHothers, ak8_probQCD in zip(*events):
        if ievent % 1000 == 0:
            print('%d events processed, %d/%d AK8 jets.' % (ievent, len(ak8_probHbc_list), iak8))
        ievent += 1
        for jet_pt, jet_eta, jet_phi, jet_mass, jet_probHbc, jet_probHothers, jet_probQCD \
            in zip(ak8_pt, ak8_eta, ak8_phi, ak8_mass, ak8_probHbc, ak8_probHothers, ak8_probQCD):
            iak8 += 1
            if not match_wcb(gp_pt, gp_eta, gp_phi, gp_mass, gp_mom, gp_pid,
                             jet_pt, jet_eta, jet_phi, jet_mass, 0.8): continue
            ak8_probHbc_list.append(jet_probHbc)
            ak8_probHothers_list.append(jet_probHothers)
            ak8_probQCD_list.append(jet_probQCD)
    print('Truth matching done!')

    ak8_probHbc, ak8_probHothers, ak8_probQCD = map(
        np.array,
        (ak8_probHbc_list, ak8_probHothers_list, ak8_probQCD_list)
    )
    print('%d events processed, %d/%d AK8 jets.' % (ievent, len(ak8_probHbc_list), iak8))
    ak8_HbcVSQCD = 1 / (1 + ak8_probQCD / ak8_probHbc)
    ak8_probHbc, ak8_probHothers, ak8_probQCD, ak8_HbcVSQCD = map(
        lambda x: np.histogram(x, bins=bins)[0].astype('float'),
        (ak8_probHbc, ak8_probHothers, ak8_probQCD, ak8_HbcVSQCD)
    )

    plt.figure(figsize=(12, 9), dpi=150)
    try:
        hep.cms.label(data=False, paper=False, supplementary=False, year=2018, lumi=59.7)
    except Exception:
        hep.cms.label(data=False, label='Preliminary', year=2018, lumi=59.7)
    hep.histplot((ak8_probHbc, bins), histtype='step', label='probHbc')
    hep.histplot((ak8_probHothers, bins), linestyle='--', histtype='step', label='probHothers')
    hep.histplot((ak8_probQCD, bins), linestyle='--', histtype='step', label='probQCD')
    hep.histplot((ak8_HbcVSQCD, bins), histtype='step', label='HbcVSQCD')
    #plt.yscale('log')
    plt.xlabel('ParT Score')
    plt.ylabel('Unweighted Events')
    plt.legend(loc='upper right')
    plt.grid()
    plt.tight_layout()
    plt.savefig('plot-%s.pdf' % os.path.basename(os.path.dirname(rootdir)))
    plt.close()

    for hist, h in zip(hists, (ak8_probHbc, ak8_probHothers, ak8_probQCD, ak8_HbcVSQCD)):
        hist += h
    for hist, h in zip(weighted_hists, (ak8_probHbc, ak8_probHothers, ak8_probQCD, ak8_HbcVSQCD)):
        hist += h * (weights[rootdir] / h.sum())
    if rootdir.find('200to400') >= 0: continue
    for hist, h in zip(weighted_hists_400toInf, (ak8_probHbc, ak8_probHothers, ak8_probQCD, ak8_HbcVSQCD)):
        hist += h * (weights[rootdir] / h.sum())

plt.figure(figsize=(12, 9), dpi=150)
try:
    hep.cms.label(data=False, paper=False, supplementary=False, year=2018, lumi=59.7)
except Exception:
    hep.cms.label(data=False, label='Preliminary', year=2018, lumi=59.7)
hep.histplot((hists[0], bins), histtype='step', label='probHbc')
hep.histplot((hists[1], bins), linestyle='--', histtype='step', label='probHothers')
hep.histplot((hists[2], bins), linestyle='--', histtype='step', label='probQCD')
hep.histplot((hists[3], bins), histtype='step', label='HbcVSQCD')
#plt.yscale('log')
plt.xlabel('ParT Score')
plt.ylabel('Unweighted Events')
plt.legend(loc='upper right')
plt.grid()
plt.tight_layout()
plt.savefig('plot.pdf')
plt.close()

plt.figure(figsize=(12, 9), dpi=150)
try:
    hep.cms.label(data=False, paper=False, supplementary=False, year=2018, lumi=59.7)
except Exception:
    hep.cms.label(data=False, label='Preliminary', year=2018, lumi=59.7)
hep.histplot((weighted_hists[0], bins), histtype='step', label='probHbc')
hep.histplot((weighted_hists[1], bins), linestyle='--', histtype='step', label='probHothers')
hep.histplot((weighted_hists[2], bins), linestyle='--', histtype='step', label='probQCD')
hep.histplot((weighted_hists[3], bins), histtype='step', label='HbcVSQCD')
#plt.yscale('log')
plt.xlabel('ParT Score')
plt.ylabel('Events')
plt.legend(loc='upper right')
plt.grid()
plt.tight_layout()
plt.savefig('plot-weighted.pdf')
plt.close()

plt.figure(figsize=(12, 9), dpi=150)
try:
    hep.cms.label(data=False, paper=False, supplementary=False, year=2018, lumi=59.7)
except Exception:
    hep.cms.label(data=False, label='Preliminary', year=2018, lumi=59.7)
hep.histplot((weighted_hists_400toInf[0], bins), histtype='step', label='probHbc')
hep.histplot((weighted_hists_400toInf[1], bins), linestyle='--', histtype='step', label='probHothers')
hep.histplot((weighted_hists_400toInf[2], bins), linestyle='--', histtype='step', label='probQCD')
hep.histplot((weighted_hists_400toInf[3], bins), histtype='step', label='HbcVSQCD')
#plt.yscale('log')
plt.xlabel('ParT Score')
plt.ylabel('Events')
plt.legend(loc='upper right')
plt.grid()
plt.tight_layout()
plt.savefig('plot-weighted-400toInf.pdf')
plt.close()
