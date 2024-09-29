#!/usr/bin/env python3

import re
import uproot
import glob
import builtins
import numpy as np
import awkward as ak
import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Transcript printed content to a same-name log file.
logfile = open(re.sub(r'\.py$', '.log', __file__), 'w')
def print(*args, **kwargs): builtins.print(*args, **kwargs); builtins.print(*args, **{**kwargs, 'file': logfile})

plt.style.use(hep.style.CMS)
gs = gridspec.GridSpec(3, 1, height_ratios=[4, 1, 1])

NEVENT_MAX = None
#NEVENT_MAX = 1000000

event_expressions = list(map(lambda x: (x[0], re.sub(r'\s+', ' ', x[1])), [
    ('weight',   '''weight'''),
    ('puWeight', '''puWeight'''),
    ('nb',       '''nb_t_deep_ex''')  # Tight. Exclusive or inclusive?
]))
jet_expressions   = list(map(lambda x: (x[0], re.sub(r'\s+', ' ', x[1])), [
    #('%s_pT',          '''PTj_V2_%s'''),
    #('%s_eta',         '''Etaj_V2_%s'''),
    #('%s_phi',         '''Phij_V2_%s'''),
    ('%s_sdmass',      '''Mj_V2_%s'''),
    ('%s_probHbc',     '''%s_Hbc'''),
    ('%s_probHcs',     '''%s_Hcs'''),
    #('%s_probHothers', '''%s_Hbb +
    #                   %s_Hbs +
    #                   %s_Hcc +
    #                   %s_Hqq +
    #                   %s_Hss'''),
    ('%s_probQCD',     '''%s_QCDb +
                       %s_QCDbb +
                       %s_QCDc +
                       %s_QCDcc +
                       %s_QCDothers'''),
]))
#jet_labels = ['a', 'b', 'c']
jet_labels = ['a']
labels = {
    'Wcb':   r'W + Jets and Top ($W \to cb$)',
    'QCD':   r'QCD',
    'WJets': r'W + Jets ($W \to \mathrm{others}$)',
    'Top':   r'Top',
    #'TT':    r'TTbar',
    #'ST':    r'SingleTop',
    'Rest':  r'Others',
}
mcfiles = glob.glob('samples/2016APV/mc/*.root')
datafiles = glob.glob('samples/2016APV/data/SlimmedTree_*.root')
mcfile_pattern = re.compile(r'^.*Tree_(.*)\.root$')

# Compute expressions to be evaluated on input ROOT files.
expressions = [expression[1] for expression in event_expressions]
for jet_label in jet_labels:
    expressions.extend([expression[1].replace('%s', jet_label) for expression in jet_expressions])
print('Expressions:', *expressions, sep='\n  - ')

# Evaluate expressions on specified input ROOT files. Read at most n events.
def concatenate(files, expressions, n=None):
    events = []
    for file in files:
        if n is not None and n == 0: break
        with uproot.open(file) as file:
            events.append(file.arrays(expressions=expressions, entry_start=0, entry_stop=n))
        if n is not None: n -= len(events)
    events = ak.concatenate(events)
    aliased_events = ak.Array({ alias: events[origin] for alias, origin in event_expressions if origin in events.fields })
    for jet_label in jet_labels:
        for alias, origin in jet_expressions:
            aliased_events[alias.replace('%s', jet_label)] = events[origin.replace('%s', jet_label)]
        aliased_events['%s_HbcVSQCS' % jet_label] = 1.0 / (1.0 + (
            0.997032 * aliased_events['%s_probQCD' % jet_label] +
            0.002968 * aliased_events['%s_probHcs' % jet_label]
        ) / aliased_events['%s_probHbc' % jet_label])
    return aliased_events

# Categorized events.
events = { 'Wcb': [] }

for mcfile in mcfiles:
    # Extract category from pathname.
    match = mcfile_pattern.search(mcfile)
    if not match: continue
    category = match.group(1)
    if category == 'Signal': category = 'Wcb'

    # Split Wcb and non-Wcb events.
    print('Loading %s events...' % category)
    current_events = concatenate([mcfile + ':PKUTree'], expressions, NEVENT_MAX)
    print('%d events loaded from %s.' % (len(current_events), len(current_events)))
    events[category] = current_events

# Load data.
print('Loading data...')
data_expressions = expressions.copy()
data_expressions.remove('puWeight')
events['data'] = concatenate([datafile + ':PKUTree' for datafile in datafiles], data_expressions, NEVENT_MAX)
print('Blinding data...')
events['data'] = ak.concatenate([
    events['data'][events['data']['a_sdmass'] <  50],
    events['data'][events['data']['a_sdmass'] > 110],
])

categories = list(events.keys())

def figure(*args, **kwargs):
    fig = plt.figure(*args, **kwargs)
    fig.add_subplot(gs[0])
    try:
        hep.cms.label(data=False, paper=False, supplementary=False, year='2016APV', lumi=19.52)
    except Exception:
        hep.cms.label(data=False, label='Preliminary', year='2016APV', lumi=19.52)
    return fig

def histplot(hists, cates):
    data_hist = None
    hs, cs = [], []
    for hist, cate in zip(hists, cates):
        if cate == 'data': data_hist = hist; continue
        hs.append(hist); cs.append(cate)
    hists, cates = hs, cs
    counts     = [hist[0] for (hist, cate) in zip(hists, cates) if cate != 'Wcb']
    bins       = [hist[1] for (hist, cate) in zip(hists, cates) if cate != 'Wcb']
    wcb_hists  = [hist    for (hist, cate) in zip(hists, cates) if cate == 'Wcb']
    wcb_cates  = ['Wcb']
    cates      = [cate for cate in cates if cate != 'Wcb']
    count_sums = [np.sum(count) for count in counts]
    items = sorted(zip(count_sums, cates, counts, bins))
    cates      = [item[1]            for item in items]
    hists      = [(item[2], item[3]) for item in items]
    hep.histplot(hists,                stack=True,  histtype='fill',     label=[labels[cate] for cate in cates    ], edgecolor='black', linewidth=0.5)
    plt.errorbar(np.mean([bins[0][:-1], bins[0][1:]], axis=0), np.sum(counts, axis=0), yerr=np.sqrt(np.sum(counts, axis=0)), linestyle='', color='black')
    hep.histplot(wcb_hists, yerr=True, stack=False, histtype='step',     label=[labels[cate] for cate in wcb_cates], color='black')
    if not data_hist: return
    hep.histplot(data_hist, yerr=True, stack=False, histtype='errorbar', label='data', color='black')

def savefig(path, *args, **kwargs):
    print('Saving to %s...' % path)
    plt.savefig(path, *args, **kwargs)

def get_signif(s, b):
    return np.sqrt(2 * ((s + b) * np.log(1 + s / (b + (s == 0))) - s))

def signif(hists, cates):
    hs, cs = [], []
    for hist, cate in zip(hists, cates):
        if cate == 'data': continue
        hs.append(hist); cs.append(cate)
    hists, cates = hs, cs
    wcb_hists  = [hist    for (hist, cate) in zip(hists, cates) if cate == 'Wcb']
    hists      = [hist    for (hist, cate) in zip(hists, cates) if cate != 'Wcb']
    bins = wcb_hists[0][1]
    for hist in wcb_hists + hists:
        if np.any(hist[1] != bins): raise NotImplementedError('rebinning not implemented')
        if hist[0].shape[0] + 1 != hist[1].shape[0]: raise ValueError('incompatible bin count size')
    s_cumsum = np.sum([np.cumsum(np.concatenate([[0], count])) for (count, _) in wcb_hists], axis=0)
    b_cumsum = np.sum([np.cumsum(np.concatenate([[0], count])) for (count, _) in hists    ], axis=0)
    signif = np.zeros((len(bins), len(bins)))
    for imin in range(len(bins) - 1):
        for imax in range(imin + 1, len(bins)):
            s = s_cumsum[imax] - s_cumsum[imin]
            b = b_cumsum[imax] - b_cumsum[imin]
            signif[imin, imax] = get_signif(s, b)
    signif_l = signif.max(axis=1)
    signif_u = signif.max(axis=0)
    l = signif_l.argmax()
    u = signif_u.argmax()
    signif_max = signif_l[l]
    plt.plot(bins, signif_l, label='lower')
    plt.plot(bins, signif_u, label='upper')
    plt.plot([bins[l]] * 2, [0, signif_max * 1.2], 'k--')
    plt.plot([bins[u]] * 2, [0, signif_max * 1.2], 'k--')
    plt.plot([], [], 'k--', label='optimal (%.3f)' % signif_max)
    plt.legend()
    return signif_max

def data_over_mc(hists, cates):
    data_hist = None
    hs, cs = [], []
    for hist, cate in zip(hists, cates):
        if cate == 'data': data_hist = hist; continue
        hs.append(hist); cs.append(cate)
    hists, cates = hs, cs
    bins = data_hist[1]
    d = data_hist[0]
    m = np.sum([hist[0] for hist in hists], axis=0)
    dom = d / (m + (d == 0))
    err_dom = dom * np.hypot(1 / np.sqrt(d + (d == 0)), 1 / np.sqrt(m + (m == 0)))
    x = np.mean([bins[:-1], bins[1:]], axis=0)
    plt.ylim(-0.2, 1.8)
    plt.errorbar(x, dom, yerr=err_dom, fmt='ko')

print('B veto:')
cut_events = { }
for category in categories:
    cut_events[category] = events[category][events[category]['nb'] == 0]
    print('  - %s:\t%d' % (category, len(cut_events[category])))
events = cut_events

fig = figure(figsize=(12, 13.5), dpi=150)
sdmass_bins = np.linspace(20, 220, 21)
sdmass_hists = [np.histogram(events[category]['a_sdmass'].to_numpy(), sdmass_bins, weights=events[category]['weight'].to_numpy()) for category in categories]
histplot(sdmass_hists, categories)
plt.ylabel('Events'); plt.yscale('log'); plt.legend(); plt.grid()
plt.gca().set_xticklabels([]); fig.add_subplot(gs[1])
s = signif(sdmass_hists, categories)
plt.ylabel('Significance'); plt.grid()
plt.gca().set_xticklabels([]); fig.add_subplot(gs[2])
data_over_mc(sdmass_hists, categories)
plt.xlabel('Soft Dropped Mass [GeV]'); plt.ylabel('Data/MC'); plt.grid()
plt.tight_layout(); savefig('all.pdf')
plt.close()

thresholds = [
    [0.990, 0.995, 0.998, 0.999],
    [0.800, 0.900, 0.950, 0.980],
    [0.000, 0.500, 0.800, 0.900],
]

for iSR in range(len(thresholds)):

    s_best = -1.0
    for threshold in thresholds[iSR]:

        print('Cut HbcVSQCS >= %.3f:' % threshold)
        cut_events = { }
        for category in categories:
            cut_events[category] = events[category][events[category]['a_HbcVSQCS'] >= threshold]
            print('  - %s:\t%d' % (category, len(cut_events[category])))

        fig = figure(figsize=(12, 13.5), dpi=150)
        sdmass_bins = np.linspace(20, 220, 21)
        sdmass_hists = [np.histogram(cut_events[category]['a_sdmass'].to_numpy(), sdmass_bins, weights=cut_events[category]['weight'].to_numpy()) for category in categories]
        histplot(sdmass_hists, categories)
        plt.ylabel('Events'); plt.yscale('log'); plt.legend(); plt.grid()
        plt.gca().set_xticklabels([]); fig.add_subplot(gs[1])
        s = signif(sdmass_hists, categories)
        plt.ylabel('Significance'); plt.grid()
        plt.gca().set_xticklabels([]); fig.add_subplot(gs[2])
        data_over_mc(sdmass_hists, categories)
        plt.xlabel('Soft Dropped Mass [GeV]'); plt.ylabel('Data/MC'); plt.grid()
        plt.tight_layout(); savefig('sr%d-%.3f.pdf' % (iSR + 1, threshold))
        plt.close()

        fig = plt.figure(figsize=(12, 9), dpi=150)
        try: hep.cms.label(data=False, paper=False, supplementary=False, year='2016APV', lumi=19.52)
        except Exception: hep.cms.label(data=False, label='Preliminary', year='2016APV', lumi=19.52)
        weight_bins = np.logspace(-3, 3, 61)
        weight_hists = [np.histogram(cut_events[category]['weight'].to_numpy(), weight_bins) for category in categories]
        histplot(weight_hists, categories)
        plt.xlabel('Weight'); plt.ylabel('Unweighted events'); plt.xscale('log'); plt.yscale('log'); plt.legend(); plt.grid()
        plt.tight_layout(); savefig('sr%d-%.3f_weight.pdf' % (iSR + 1, threshold))
        plt.close()

        fig = plt.figure(figsize=(12, 9), dpi=150)
        try: hep.cms.label(data=False, paper=False, supplementary=False, year='2016APV', lumi=19.52)
        except Exception: hep.cms.label(data=False, label='Preliminary', year='2016APV', lumi=19.52)
        weight_bins = np.logspace(-3, 3, 61)
        weight_hists = [np.histogram(cut_events[category]['puWeight'].to_numpy(), weight_bins) for category in categories if category != 'data']
        histplot(weight_hists, categories)
        plt.xlabel('PU Weight'); plt.ylabel('Unweighted events'); plt.xscale('log'); plt.yscale('log'); plt.legend(); plt.grid()
        plt.tight_layout(); savefig('sr%d-%.3f_puweight.pdf' % (iSR + 1, threshold))
        plt.close()

        if s > s_best: s_best = s; threshold_best = threshold

    threshold = threshold_best
    print('Cut HbcVSQCS < %.3f:' % threshold)
    cut_events = { }
    for category in categories:
        cut_events[category] = events[category][events[category]['a_HbcVSQCS'] < threshold]
        print('  - %s:\t%d' % (category, len(cut_events[category])))
    events = cut_events
