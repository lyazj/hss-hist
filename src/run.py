import mplhep as hep
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.style.use(hep.style.CMS)

def get_significance(s, b):

    return np.sqrt(2 * ((s + b) * np.log(1 + s / (b + (s == 0))) - s))

def sort_names_and_counts(names, counts):
    items = sorted(((np.sum(count), name, count) for (name, count) in zip(names, counts)))
    return [item[1] for item in items], [item[2] for item in items]

def run(config):

    # An expression produces values to fill a histogram.
    if not config['hists']: return
    expressions = [hist['expr'] for hist in config['hists']]
    lines_all_categories = [[] for hist in config['hists']]

    # Draw a group of histograms for each category.
    nsg, nbg, wsg, wbg, wsg_before, wbg_before = 0, 0, 0.0, 0.0, 0.0, 0.0
    for category in config['categories']:
        lines = [list(np.histogram([], bins=hist['nbin'], range=(hist['lb'], hist['ub']), weights=[]))
                 for hist in config['hists']]  # 0: counts, 1: bins
        windows = [hist['window'] for hist in config['hists']]
        name = category['name']
        xs = 0.0
        nevent = 0
        nvalid_sum = 0
        weight_sum = 0.0
        weight_sum_before = 0.0

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
            weight_sum_before += len(values[0]) * weight
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
            wsg_before += weight_sum_before
        else:
            nbg += nvalid_sum
            wbg += weight_sum
            wbg_before += weight_sum_before
    sig = get_significance(wsg, wbg)
    print('Summary: nsg=%d nbg=%d wsg=%f(%f) wbg=%f(%f) sig=%f' % (nsg, nbg, wsg, wsg / (wsg_before + (wsg == 0)), wbg, wbg / (wbg_before + (wbg == 0)), sig))

    # Draw and export histograms.
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
        sg_names, sg_counts = sort_names_and_counts(sg_names, sg_counts)
        bg_names, bg_counts = sort_names_and_counts(bg_names, bg_counts)
        fig = plt.figure(figsize=hist['figsize'], dpi=hist['dpi'])
        gs = gridspec.GridSpec(hist['nsubplot-y'], hist['nsubplot-x'],
                               wspace=hist['subplot-space-x'], hspace=hist['subplot-space-y'],
                               width_ratios=hist['subplot-ratios-x'], height_ratios=hist['subplot-ratios-y'])
        gsi = 0
        fig.add_subplot(gs[gsi])
        gsi += 1
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

        # Draw significance subplots.
        if hist['subplot-significance-lower'] or hist['subplot-significance-upper']:
            plt.xlabel('')
            plt.gca().set_xticklabels([])
            fig.add_subplot(gs[gsi])
            gsi += 1
            if hist['subplot-significance-lower']:
                csg = np.sum(sg_counts, axis=0)
                cbg = np.sum(bg_counts, axis=0)
                csg = np.concatenate([np.cumsum(csg[::-1])[::-1], [0.0]])
                cbg = np.concatenate([np.cumsum(cbg[::-1])[::-1], [0.0]])
                sig = get_significance(csg, cbg)
                plt.plot(bins, sig, '*-k')
            if hist['subplot-significance-upper']:
                csg = np.sum(sg_counts, axis=0)
                cbg = np.sum(bg_counts, axis=0)
                csg = np.concatenate([[0.0], np.cumsum(csg)])
                cbg = np.concatenate([[0.0], np.cumsum(cbg)])
                sig = get_significance(csg, cbg)
                plt.plot(bins, sig, '^-k')
            plt.xlabel(hist['xlabel'])
            plt.ylabel('significance')
            plt.xscale(hist['xscale'])
            if hist['grid']: plt.grid(axis=hist['grid'])
            plt.tight_layout()

        # Export the plots.
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
