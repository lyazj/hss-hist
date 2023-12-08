import re
import os
import yaml
import glob

basedir = os.path.join(os.path.dirname(__file__), '..')

class Config(dict):

    # Load configuration from yaml file.
    def __init__(self, filename, override={}):

        dict.__init__(self, yaml.load(open(filename), yaml.Loader))
        self.update(override)
        self.preproc()

    # Pre-process configuration.
    def preproc(self):

        # Add branch prefixes on demand.
        branch_prefix = self['branch-prefix']
        for branches in 'signal-branches', 'background-branches':
            branches = self[branches]
            for ibranch, branch in enumerate(branches):
                if branch.find(branch_prefix) != 0:
                    branches[ibranch] = branch_prefix + branch

        # Compute weight and list files for each sample on demand.
        luminosity = self['luminosity']
        for category in self['categories']:
            for sample in category['samples']:
                if 'weight' not in sample:
                    sample['weight'] = sample['xs'] * 1e3 * luminosity / sample['nevent']
                if 'files' not in sample:
                    sample['files'] = glob.glob(os.path.join(self['sample-dir'], '%s_*_tree.root' % sample['name']))

        # Generate tagger determinate expression on demand.
        if 'expr' not in self['hists'][0]:
            tagger_expr = '1 / (1 + (%s) / (%s))' % \
                    (' + '.join(self['background-branches']), ' + '.join(self['signal-branches']))
            self['hists'][0]['expr'] = tagger_expr

        # Complete histogram attributes.
        for hist in self['hists']:
            hist['xscale'] = hist.get('xscale', None) or 'linear'
            hist['yscale'] = hist.get('yscale', None) or 'linear'

        # List active branches on demand.
        if 'active-branches' not in self:
            active_branches = set()
            for hist in self['hists']:
                for branch in re.findall(r'[A-Za-z_][A-Za-z0-9_]*', hist['expr']):
                    active_branches.add(branch)
            self['active-branches'] = sorted(active_branches)

if __name__ == '__main__':

    # A concise example.
    config = Config(os.path.join(basedir, 'run', '2018', '1L', 'mc', 'zss_part', 'config.yaml'),
                    {'sample-dir': os.path.join(basedir, 'example')})
    print(config)
