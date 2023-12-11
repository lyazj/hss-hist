import re
import os
import yaml
import merge
import math

basedir = os.path.join(os.path.dirname(__file__), '..')

class Config(dict):

    # Load configuration from yaml file.
    def __init__(self, filename, override={}):

        dict.__init__(self, yaml.load(open(filename), yaml.Loader))
        self.update(override)
        self.preproc()

    # Pre-process configuration.
    def preproc(self):

        # Compute weight and list files for each sample on demand.
        candidate_files = None
        luminosity = self['luminosity']
        for category in self['categories']:
            for sample in category['samples']:
                if 'weight' not in sample:
                    sample['weight'] = sample['xs'] * 1e3 * luminosity / sample['nevent']
                if 'files' not in sample:
                    if candidate_files is None:
                        # Group <self['sample-dir']>/<name>_<id>_tree.root by name.
                        candidate_files = {}
                        files = map(lambda s: s.rsplit('_', 2), os.listdir(self['sample-dir']))
                        files = sorted((f[0], int(f[1])) for f in files if len(f) == 3 and f[2] == 'tree.root')
                        for fname, fid in files:
                            if fname not in candidate_files: candidate_files[fname] = []
                            candidate_files[fname].append(fid)
                    sample['files'] = [os.path.join(self['sample-dir'], '_'.join([sample['name'], str(fid), 'tree.root']))
                                       for fid in candidate_files.get(sample['name'], [])]
                if not self.get('do-not-merge') and not sample.get('do-not-merge') and 'merged-file' not in sample:
                    sample['merged-file'] = os.path.join(self['merged-sample-dir'], '_'.join([sample['name'], 'bigtree.root']))
                    merge.update_merged_root_file(sample['files'], sample['merged-file'])

        # Complete histogram attributes.
        for hist in self['hists']:
            hist['preliminary'] = hist.get('preliminary', True)
            hist['supplementary'] = hist.get('supplementary', False)
            hist['window'] = eval(hist.get('window', '[-inf, +inf]'), {'inf': math.inf})
            hist['xlabel'] = hist.get('xlabel', hist['name'])
            hist['ylabel'] = hist.get('ylabel', 'number')
            hist['xscale'] = hist.get('xscale', 'linear')
            hist['yscale'] = hist.get('yscale', 'linear')
            hist['grid'] = hist.get('grid')
            hist['legend-options'] = hist.get('legend-options', {})
            hist['stack'] = hist.get('stack', False)
            hist['no-stack-signal'] = hist.get('no-stack-signal', False)
            hist['format'] = hist.get('format', 'pdf').split(',')
            hist['subplot-significance'] = hist.get('subplot-significance', True)
            hist['nsubplot-x'] = hist.get('nsubplot-x', 1)
            hist['nsubplot-y'] = hist.get('nsubplot-y', 1 + hist['subplot-significance'])
            hist['subplot-ratios-x'] = hist.get('subplot-ratios-x', [1] * hist['nsubplot-x'])
            hist['subplot-ratios-y'] = hist.get('subplot-ratios-y', [4] + [1] * (hist['nsubplot-y'] - 1))

        # List active branches on demand.
        if 'active-branches' not in self:
            re_identifier = re.compile(r'[A-Za-z_][A-Za-z0-9_]*')
            active_branches = set()
            for hist in self['hists']:
                for branch in re_identifier.findall(hist['expr']):
                    active_branches.add(branch)
            self['active-branches'] = sorted(active_branches)

        # Global configuration.
        self['figsize'] = self.get('figsize', (12, 9))
        self['dpi'] = self.get('dpi', 150)

if __name__ == '__main__':

    # A concise example.
    config = Config(os.path.join(basedir, 'run', '2018', '1L', 'mc', 'zss_part', 'config.yaml'),
                    {'sample-dir': os.path.join(basedir, 'example'), 'do-not-merge': True})
    print(config)
