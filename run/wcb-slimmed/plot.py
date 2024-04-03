#!/usr/bin/env python3.6

import re
import uproot
import glob
import builtins
import numpy as np
import awkward as ak
import mplhep as hep
import matplotlib.pyplot as plt
import pickle

#NEVENT_MAX = 1000000
NEVENT_MAX = None

logfile = open(re.sub(r'\.py$', '.log', __file__), 'w')
def print(*args, **kwargs): builtins.print(*args, **kwargs); builtins.print(*args, **{**kwargs, 'file': logfile})

plt.style.use(hep.style.CMS)
event_expressions = list(map(lambda x: (x[0], re.sub(r'\s+', ' ', x[1])), [
    ('isWcb',       '''isWcb'''       ),  # Must be the first item.
    ('weight',      '''weight'''      ),
]))
jet_expressions   = list(map(lambda x: (x[0], re.sub(r'\s+', ' ', x[1])), [
    ('pT',          '''PTj_V2_%s'''   ),
    ('eta',         '''Etaj_V2_%s'''  ),
    ('phi',         '''Phij_V2_%s'''  ),
    ('sdmass',      '''Mj_V2_%s'''    ),
    ('probHbc',     '''%s_Hbc'''      ),
    #('probHcs',     '''%s_Hcs'''      ),
    #('probHothers', '''%s_Hbb +
    #                   %s_Hbs +
    #                   %s_Hcc +
    #                   %s_Hee +
    #                   %s_Hgg +
    #                   %s_Hmm +
    #                   %s_Hqq +
    #                   %s_Hss +
    #                   %s_Htauhtaue +
    #                   %s_Htauhtauh +
    #                   %s_Htauhtaum'''),
    ('probQCD',     '''%s_QCDb +
                       %s_QCDbb +
                       %s_QCDc +
                       %s_QCDcc +
                       %s_QCDothers'''),
]))
jet_labels = ['a', 'b', 'c']
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

expressions = [expression[1] for expression in event_expressions]
for jet_label in jet_labels:
    expressions.extend([expression[1].replace('%s', jet_label) for expression in jet_expressions])
print('Expressions:', *expressions, sep='\n  - ')

              #       dict      dict     Array
events = { }  # events[category][feature][ievent]

              #       dict      list  dict     Array
jets   = { }  # jets  [category][ijet][feature][ievent]

def concatenate(files, expressions, n=None):
    values = [[] for expression in expressions]
    if not values: return values
    for file in files:
        if n is not None and n == 0: break
        with uproot.open(file) as file:
            for i, value in enumerate(file.arrays(expressions=expressions, entry_start=0, entry_stop=n, how=tuple)):
                values[i].append(value)
        if n is not None: n -= len(values[0][-1])
    for i, value in enumerate(values):
        values[i] = ak.concatenate(value)
    return values

events['Wcb'] = {expression[0]: [] for expression in event_expressions}
jets['Wcb'] = [{expression[0]: [] for expression in jet_expressions} for i in range(len(jet_labels))]
for rootfile in rootfiles:
    match = rootfile_pattern.search(rootfile)
    if not match: continue
    category = match.group(1)
    if category not in labels: continue

    print('Loading %s events...' % category)
    values = concatenate([rootfile + ':PKUTree'], expressions, NEVENT_MAX)
    n = len(values[0])
    wcb_values = [value[values[0] == True ] for value in values]
    values     = [value[values[0] == False] for value in values]
    print('%d events (%d Wcb) loaded from %s.' % (n, len(wcb_values[0]), rootfile))

    events[category] = {expression[0]: values[i] for (i, expression) in enumerate(event_expressions)}
    for (i, expression) in enumerate(event_expressions): events['Wcb'][expression[0]].append(wcb_values[i])
    values = values[len(events[category]) :]
    wcb_values = wcb_values[len(events[category]) :]

    jets[category] = [ ]
    for i in range(len(jet_labels)):
        jets[category].append({expression[0]: values[i] for (i, expression) in enumerate(jet_expressions)})
        for (j, expression) in enumerate(jet_expressions): jets['Wcb'][i][expression[0]].append(wcb_values[j])
        values = values[len(jets[category][-1]) :]
        wcb_values = wcb_values[len(jets[category][-1]) :]
labels = {category: labels[category] for category in events}
for feature in events['Wcb']: events['Wcb'][feature] = ak.concatenate(events['Wcb'][feature])
for i in range(len(jet_labels)):
    for feature in jets['Wcb'][i]: jets['Wcb'][i][feature] = ak.concatenate(jets['Wcb'][i][feature])

for category in labels.keys():
    for i in range(len(jet_labels)):
        jets[category][i]['HbcVSQCD'] = 1.0 / (1.0 + jets[category][i]['probQCD'] / jets[category][i]['probHbc'])

def figure(*args, **kwargs):
    plt.figure(*args, **kwargs)
    try:
        hep.cms.label(data=False, paper=False, supplementary=False, year=2018, lumi=59.7)
    except Exception:
        hep.cms.label(data=False, label='Preliminary', year=2018, lumi=59.7)

def histplot(hists, cates):
    counts     = [hist[0] for (hist, cate) in zip(hists, cates) if cate != 'Wcb']
    bins       = [hist[1] for (hist, cate) in zip(hists, cates) if cate != 'Wcb']
    wcb_hists  = [hist    for (hist, cate) in zip(hists, cates) if cate == 'Wcb']
    wcb_cates  = ['Wcb']
    cates      = [cate for cate in cates if cate != 'Wcb']
    count_sums = [ak.sum(count) for count in counts]
    items = sorted(zip(count_sums, cates, counts, bins))
    cates      = [item[1]            for item in items]
    hists      = [(item[2], item[3]) for item in items]
    hep.histplot(hists,     stack=True,  histtype='fill', label=[labels[cate] for cate in cates    ], edgecolor='black', linewidth=0.5)
    hep.histplot(wcb_hists, stack=False, histtype='step', label=[labels[cate] for cate in wcb_cates], color='black')

plot = {'labels': labels}

figure(figsize=(12, 9), dpi=150)
sdmass_bins = np.linspace(20, 220, 51)
sdmass_hists = [np.histogram(jets[category][0]['sdmass'], sdmass_bins, weights=events[category]['weight']) for category in labels.keys()]
histplot(sdmass_hists, labels.keys())
plt.xlabel('Soft Dropped Mass [GeV]'); plt.ylabel('Events'); plt.yscale('log')
plt.legend(); plt.grid(); plt.tight_layout(); plt.savefig('plot-sdmass.pdf')
plt.close()
plot['sdmass'] = sdmass_hists

figure(figsize=(12, 9), dpi=150)
HbcVSQCD_bins = np.linspace(0.0, 1.0, 51)
HbcVSQCD_hists = [np.histogram(jets[category][0]['HbcVSQCD'], HbcVSQCD_bins, weights=events[category]['weight']) for category in labels.keys()]
histplot(HbcVSQCD_hists, labels.keys())
plt.xlabel('HbcVSQCD'); plt.ylabel('Events'); plt.yscale('log')
plt.legend(); plt.grid(); plt.tight_layout(); plt.savefig('plot-HbcVSQCD.pdf')
plt.close()
plot['HbcVSQCD'] = HbcVSQCD_hists

figure(figsize=(12, 9), dpi=150)
HbcVSQCD_bins = np.linspace(0.95, 1.0, 51)
HbcVSQCD_hists = [np.histogram(jets[category][0]['HbcVSQCD'], HbcVSQCD_bins, weights=events[category]['weight']) for category in labels.keys()]
histplot(HbcVSQCD_hists, labels.keys())
plt.xlabel('HbcVSQCD'); plt.ylabel('Events'); plt.yscale('log')
plt.legend(); plt.grid(); plt.tight_layout(); plt.savefig('plot-HbcVSQCD-0.95.pdf')
plt.close()
plot['HbcVSQCD-0.95'] = HbcVSQCD_hists

def cut(event, jet, cut):
    for feature, value in event.items(): event[feature] = value[cut]
    for i in range(len(jet_labels)):
        for feature, value in jet[i].items(): jet[i][feature] = value[cut]

for threshold in sorted([0.95, 0.98, 0.99, 0.995, 0.998, 0.999]):
    print('Cut HbcVSQCD >= %.3f:' % threshold)
    for category in labels.keys():
        cut(events[category], jets[category], jets[category][0]['HbcVSQCD'] >= threshold)
        print('  - %s:\t%d' % (category, len(events[category]['weight'])))
    figure(figsize=(12, 9), dpi=150)
    sdmass_bins = np.linspace(20, 220, 51)
    sdmass_hists = [np.histogram(jets[category][0]['sdmass'], sdmass_bins, weights=events[category]['weight']) for category in labels.keys()]
    histplot(sdmass_hists, labels.keys())
    plt.xlabel('Soft Dropped Mass [GeV]'); plt.ylabel('Events'); plt.yscale('log')
    plt.legend(); plt.grid(); plt.tight_layout(); plt.savefig('plot-sdmass-%.3f.pdf' % threshold)
    plt.close()
    plot['sdmass-%.3f' % threshold] = sdmass_hists

pickle.dump(plot, open('plot.pkl', 'wb'))
