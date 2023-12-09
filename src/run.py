import numpy as np
import uproot
import matplotlib.pyplot as plt

SIGNAL_COLORS = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white']

def run(config):

    # An expression produces values to fill a histogram.
    if not config['hists']: return
    expressions = [hist['expr'] for hist in config['hists']]
    hists_all_categories = [[] for hist in config['hists']]

    # Draw a group of histograms for each category.
    for category in config['categories']:
        hists = [list(np.histogram([], bins=hist['nbin'], range=(hist['lb'], hist['ub']), weights=[]))
                 for hist in config['hists']]  # 0: counts, 1: bins
        thresholds = [hist['threshold'] for hist in config['hists'][:2]] + [0.0]
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
            files = [file + ':Events' for file in sample['files']]
            if not files: continue

            # Evaluate expressions.
            values = uproot.concatenate(files, expressions, library='np', how=tuple, allow_missing=True)
            nevent_now += len(values[0])

            # Fill histograms. Use zero weights as masks.
            weight = np.ones_like(values[0]) * weight
            for figid, (hist, value, threshold) in enumerate(zip(hists, values, thresholds), 1):
                hist[0] += np.histogram(value, bins=hist[1], weights=weight)[0]
                weight *= value >= threshold

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
        if hist['stack-background']:
            plt.hist([bins[:-1]] * len(bg_counts), bins, weights=bg_counts, label=bg_names, histtype='barstacked', alpha=0.5)
            plt.hist([bins[:-1]] * len(sg_counts), bins, weights=sg_counts, label=sg_names, histtype='step', color=SIGNAL_COLORS[:len(sg_counts)])
        else:
            plt.hist([bins[:-1]] * len(bg_counts), bins, weights=bg_counts, label=bg_names, histtype='step')
            plt.hist([bins[:-1]] * len(sg_counts), bins, weights=sg_counts, label=sg_names, histtype='step')
        plt.xlabel(hist['name'])
        plt.ylabel('number')
        plt.xscale(hist['xscale'])
        plt.yscale(hist['yscale'])
        plt.grid(axis='y')
        plt.legend()
        plt.tight_layout()
        plt.savefig('%s-%s-%s-%s.%s' % (hist['name'], hist['xscale'], hist['yscale'], 'stacked' if hist['stack-background'] else 'step', hist['format']))

if __name__ == '__main__':

    # A concise example.
    import os
    from config import Config, basedir
    config = Config(os.path.join(basedir, 'run', '2018', '1L', 'mc', 'zss_part', 'config.yaml'),
                    {'sample-dir': os.path.join(basedir, 'example')})
    run(config)
