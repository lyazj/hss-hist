#!/usr/bin/env python3.6

import re
import uproot
import glob
import builtins
import numpy as np
import awkward as ak
import mplhep as hep
import matplotlib.pyplot as plt

# Transcript printed content to a same-name log file.
logfile = open(re.sub(r'\.py$', '.log', __file__), 'w')
def print(*args, **kwargs): builtins.print(*args, **kwargs); builtins.print(*args, **{**kwargs, 'file': logfile})

plt.style.use(hep.style.CMS)

NEVENT_MAX = None
#NEVENT_MAX = 1000000

event_expressions = list(map(lambda x: (x[0], re.sub(r'\s+', ' ', x[1])), [
    ('isWcb',  '''isWcb'''),
    ('weight', '''weight'''),
    ('nb',     '''nb_t_deep_ex''')  # Tight. Exclusive or inclusive?
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
rootfiles = glob.glob('samples/2018/mc/*.root')
rootfile_pattern = re.compile('^.*Tree_(.*)\.root$')

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
    aliased_events = ak.Array({ alias: events[origin] for alias, origin in event_expressions })
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

for rootfile in rootfiles:
    # Extract category from pathname.
    match = rootfile_pattern.search(rootfile)
    if not match: continue
    category = match.group(1)
    if category == 'Wcb': raise RuntimeError('unexpected category: Wcb')
    if category not in labels: continue

    # Split Wcb and non-Wcb events.
    print('Loading %s events...' % category)
    current_events = concatenate([rootfile + ':PKUTree'], expressions, NEVENT_MAX)
    wcb_events    = current_events[current_events['isWcb'] == True ]
    nonwcb_events = current_events[current_events['isWcb'] == False]
    print('%d events (%d Wcb) loaded from %s.' % (len(current_events), len(wcb_events), rootfile))
    events['Wcb'].append(wcb_events)
    events[category] = nonwcb_events

# Merge Wcb events from different categories.
events['Wcb'] = ak.concatenate(events['Wcb'])
categories = list(events.keys())

def figure(*args, **kwargs):
    fig = plt.figure(*args, **kwargs)
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
    hep.histplot(hists,                stack=True,  histtype='fill', label=[labels[cate] for cate in cates    ], edgecolor='black', linewidth=0.5)
    hep.histplot(wcb_hists, yerr=True, stack=False, histtype='step', label=[labels[cate] for cate in wcb_cates], color='black')

def savefig(path, *args, **kwargs):
    print('Saving to %s...' % path)
    plt.savefig(path, *args, **kwargs)

figure(figsize=(12, 9), dpi=150)
nb_bins = np.linspace(0, 8.0, 9)
nb_hists = [np.histogram(events[category]['nb'], nb_bins, weights=events[category]['weight']) for category in labels.keys()]
histplot(nb_hists, labels.keys())
plt.xlabel('$N_b$'); plt.ylabel('Events'); plt.yscale('log')
plt.legend(); plt.grid(); plt.tight_layout(); savefig('view-nb.pdf')
plt.close()
