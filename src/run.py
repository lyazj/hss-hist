import mplhep as hep
import numpy as np
import uproot
import matplotlib.pyplot as plt

plt.style.use(hep.style.CMS)

def run(config):

    if not config['hists']: return
    plt.figure(figsize=config['figsize'], dpi=config['dpi'])

    # An expression produces values to fill a histogram.
    expressions = [hist['expr'] for hist in config['hists']]
    hists_all_categories = [[] for hist in config['hists']]

    # Draw a group of histograms for each category.
    for category in config['categories']:
        hists = [list(np.histogram([], bins=hist['nbin'], range=(hist['lb'], hist['ub']), weights=[]))
                 for hist in config['hists']]  # 0: counts, 1: bins
        thresholds = [hist['threshold'] for hist in config['hists']]
        name = category['name']
        xs = 0.0
        nevent = 0
        nevent_now = 0

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
            nevent_now += len(values[0])

            # Fill histograms. Use zero weights as masks.
            weight = np.ones_like(values[0]) * weight
            for hist, value, threshold in zip(hists, values, thresholds):
                hist[0] += np.histogram(value, bins=hist[1], weights=weight)[0]
                if threshold is not None: weight *= value >= threshold

        # Store histograms.
        for hist_all_categories, hist in zip(hists_all_categories, hists):
            hist_all_categories.append((name, hist))
        print('Summary for %s: %d/%d events scaled to %f pb' % (name, nevent_now, nevent, xs))

    # Draw and export histograms.
    for hist, hist_all_categories in zip(config['hists'], hists_all_categories):
        if not hist_all_categories: continue
        bins = hist_all_categories[0][1][1]
        names = [h[0] for h in hist_all_categories]
        counts = [h[1][0] for h in hist_all_categories]
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
        if hist['stack-background']:
            hep.histplot(bg_counts, bins, label=bg_names, stack=True, histtype='fill', edgecolor='black', linewidth=0.5)
        else:
            hep.histplot(bg_counts, bins, label=bg_names)
        hep.histplot(sg_counts, bins, label=sg_names)
        plt.xlabel(hist['xlabel'])
        plt.ylabel(hist['ylabel'])
        plt.xscale(hist['xscale'])
        plt.yscale(hist['yscale'])
        if hist['grid']: plt.grid(axis=hist['grid'])
        plt.legend(**hist['legend-options'])
        plt.tight_layout()
        for extension in hist['format']:
            plt.savefig('%s-%s-%s-%s.%s' % (
                hist['name'], hist['xscale'], hist['yscale'],
                'stack' if hist['stack-background'] else 'step', extension
            ))

    plt.close()

if __name__ == '__main__':

    # A concise example.
    import os
    from config import Config, basedir
    config = Config(os.path.join(basedir, 'run', '2018', '1L', 'mc', 'zss_part', 'config.yaml'),
                    {'sample-dir': os.path.join(basedir, 'example'), 'do-not-merge': True})
    run(config)
