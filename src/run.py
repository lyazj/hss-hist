import mplhep as hep
import numpy as np
import uproot
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use(hep.style.CMS)

def get_significance(s, b):

    return np.sqrt(2 * ((s + b) * np.log(1 + s / b) - s))

def run(config):

    # An expression produces values to fill a histogram.
    if not config['hists']: return
    expressions = [hist['expr'] for hist in config['hists']]
    lines_all_categories = [[] for hist in config['hists']]

    # Draw a group of histograms for each category.
    nsg, nbg, wsg, wbg = 0, 0, 0.0, 0.0
    for category in config['categories']:
        lines = [list(np.histogram([], bins=hist['nbin'], range=(hist['lb'], hist['ub']), weights=[]))
                 for hist in config['hists']]  # 0: counts, 1: bins
        windows = [hist['window'] for hist in config['hists']]
        name = category['name']
        xs = 0.0
        nevent = 0
        nvalid_sum = 0
        weight_sum = 0.0

        # Sum up all samples.
        for sample in category['samples']:
            print('Processing: %s' % sample['name'])
            xs += sample['xs']
            nevent += sample['nevent']
            weight = sample['weight']
            files = sample.get('merged-file')
            if files:
                files = [files]
            else:
                files = sample['files']
            if not files: continue
            files = [file + ':Events' for file in files]

            # Evaluate expressions.
            values = uproot.concatenate(files, expressions, library='np', how=tuple, allow_missing=True)

            # Fill histograms. Use zero weights as masks.
            weight = np.ones_like(values[0]) * weight
            for hist, line, value, window in zip(config['hists'], lines, values, windows):
                line[0] += np.histogram(value, bins=line[1], weights=weight)[0]
                nvalid = np.sum(weight.astype('bool').astype('int'))
                print('%s: %d/%d events in %s' % (sample['name'], nvalid, weight.shape[0], hist['name']))
                weight *= np.logical_and(window[0] <= value, value <= window[1])
            nvalid = np.sum(weight.astype('bool').astype('int'))
            print('%s: %d/%d events in the end' % (sample['name'], nvalid, weight.shape[0]))
            nvalid_sum += nvalid
            weight_sum += np.sum(weight)

        # Store histograms.
        for line_all_categories, line in zip(lines_all_categories, lines):
            line_all_categories.append((name, line))
        print('Summary for %s: %d/%d events scaled to %f pb\n' % (name, nvalid_sum, nevent, xs))
        if name in config['signal-categories']:
            nsg += nvalid_sum
            wsg += weight_sum
        else:
            nbg += nvalid_sum
            wbg += weight_sum
    sig = get_significance(wsg, wbg)
    print('Summary: nsg=%d nbg=%d wsg=%f wbg=%f sig=%f' % (nsg, nbg, wsg, wbg, sig))

    # Draw and export histograms.
    plt.figure(figsize=config['figsize'], dpi=config['dpi'])
    for hist, line_all_categories in zip(config['hists'], lines_all_categories):
        if not line_all_categories: continue
        bins = line_all_categories[0][1][1]
        names = [h[0] for h in line_all_categories]
        counts = [h[1][0] for h in line_all_categories]
        sg_names, bg_names, sg_counts, bg_counts = [], [], [], []
        for name, count in zip(names, counts):
            if name in config['signal-categories']:
                sg_names.append(name); sg_counts.append(count)
            else:
                bg_names.append(name); bg_counts.append(count)
        plt.clf()
        try:
            hep.cms.label(data=not config['mc'], paper=not hist['preliminary'], supplementary=hist['supplementary'],
                          year=config['year'], lumi=config['luminosity'])
        except Exception:
            label = 'Preliminary' if hist['preliminary'] else 'Supplementary' if hist['supplementary'] else ''
            hep.cms.label(data=not config['mc'], label=label, year=config['year'], lumi=config['luminosity'])
        if hist['stack'] and hist['no-stack-signal']:
            hep.histplot(bg_counts, bins, label=bg_names, stack=True, histtype='fill', edgecolor='black', linewidth=0.5)
            hep.histplot(sg_counts, bins, label=sg_names)
        else:
            bs_counts = bg_counts + sg_counts
            bs_names = bg_names + sg_names
            if hist['stack']:
                hep.histplot(bs_counts, bins, label=bs_names, stack=True, histtype='fill', edgecolor='black', linewidth=0.5)
            else:
                hep.histplot(bs_counts, bins, label=bs_names)
        plt.xlabel(hist['xlabel'])
        plt.ylabel(hist['ylabel'])
        plt.xscale(hist['xscale'])
        plt.yscale(hist['yscale'])
        if hist['grid']: plt.grid(axis=hist['grid'])
        plt.legend(**hist['legend-options'])
        plt.tight_layout()
        for extension in hist['format']:
            plt.savefig('%s-%s-%s-%s-%s.%s' % (
                hist['name'], hist['xscale'], hist['yscale'],
                'stack' if hist['stack'] else 'step',
                'stack' if hist['stack'] and not hist['no-stack-signal'] else 'step',
                extension
            ))
    plt.close()

if __name__ == '__main__':

    # A concise example.
    import os
    from config import Config, basedir
    config = Config(os.path.join(basedir, 'run', '2018', '1L', 'mc', 'zss_part', 'config.yaml'),
                    {'sample-dir': os.path.join(basedir, 'example'), 'do-not-merge': True})
    run(config)
