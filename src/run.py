import numpy as np
import uproot
import matplotlib.pyplot as plt

def run(config):

    # An expression produces values to fill a histogram.
    expressions = [hist['expr'] for hist in config['hists']]

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

        # Draw histograms.
        for figid, hist in enumerate(hists, 1):
            plt.figure(figid)
            plt.hist(hist[1][:-1], hist[1], weights=hist[0], label=name, histtype='step')
        print('Summary for %s: %d/%d events scaled to %f pb' % (name, nevent_now, nevent, xs))

    # Export histograms.
    for figid, hist in enumerate(config['hists'], 1):
        plt.figure(figid)
        plt.xlabel(hist['name'])
        plt.ylabel('number')
        plt.xscale(hist['xscale'])
        plt.yscale(hist['yscale'])
        plt.grid(axis='y')
        plt.legend()
        plt.tight_layout()
        plt.savefig(hist['name'] + '.pdf')

if __name__ == '__main__':

    # A concise example.
    import os
    from config import Config, basedir
    config = Config(os.path.join(basedir, 'run', '2018', '1L', 'mc', 'zss_part', 'config.yaml'),
                    {'sample-dir': os.path.join(basedir, 'example')})
    run(config)
