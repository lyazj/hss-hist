#!/usr/bin/env python3

import re
import builtins
import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle

logfile = open(re.sub(r'\.py$', '.log', __file__), 'w')
def print(*args, **kwargs): builtins.print(*args, **kwargs); builtins.print(*args, **{**kwargs, 'file': logfile})

plt.style.use(hep.style.CMS)
gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])

def figure(*args, **kwargs):
    fig = plt.figure(*args, **kwargs)
    fig.add_subplot(gs[0])
    try:
        hep.cms.label(data=False, paper=False, supplementary=False, year=2018, lumi=59.7)
    except Exception:
        hep.cms.label(data=False, label='Preliminary', year=2018, lumi=59.7)
    return fig

def histplot(hists, cates):
    counts     = [hist[0] for (hist, cate) in zip(hists, cates) if cate != 'Wcb']
    bins       = [hist[1] for (hist, cate) in zip(hists, cates) if cate != 'Wcb']
    wcb_hists  = [hist    for (hist, cate) in zip(hists, cates) if cate == 'Wcb']
    wcb_cates  = ['Wcb']
    cates      = [cate for cate in cates if cate != 'Wcb']
    count_sums = [np.sum(count) for count in counts]
    items = sorted(zip(count_sums, cates, counts, bins))
    cates      = [item[1]            for item in items]
    hists      = [(item[2], item[3]) for item in items]
    hep.histplot(hists,     yerr=True, stack=True,  histtype='fill', label=[labels[cate] for cate in cates    ], edgecolor='black', linewidth=0.5)
    hep.histplot(wcb_hists, yerr=True, stack=False, histtype='step', label=[labels[cate] for cate in wcb_cates], color='black')

def savefig(path, *args, **kwargs):
    print('Saving to %s...' % path)
    plt.savefig(path, *args, **kwargs)

plot = pickle.load(open('plot-QCS.pkl', 'rb'))
labels = plot['labels']

def get_signif(s, b):
    return np.sqrt(2 * ((s + b) * np.log(1 + s / (b + (s == 0))) - s))

def signif(hists, cates):
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

print(*plot.keys())
for name in plot:
    if name[:6] != 'sdmass': continue
    fig = figure(figsize=(12, 11.25), dpi=150)
    sdmass_hists = plot[name]
    histplot(sdmass_hists, labels.keys())
    plt.ylabel('Events'); plt.yscale('log'); plt.legend(); plt.grid()
    plt.gca().set_xticklabels([]); fig.add_subplot(gs[1])
    signif(sdmass_hists, labels.keys())
    plt.xlabel('Soft Dropped Mass [GeV]'); plt.ylabel('Significance'); plt.grid()
    plt.tight_layout(); savefig('plot-QCS-%s-sig.pdf' % name)
    plt.close()
