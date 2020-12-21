from os.path import join, exists
from os import makedirs


def check_make(curr_dir, sub):
    outdir = join(curr_dir, sub)
    if not exists(curr_dir):
        makedirs(curr_dir)
    return outdir


def logger():
    # TODO
    pass


class metaQUAST(luigi.Task):
    homedir = luigi.Parameter() # TODO luigi.cfg
    refdir = luigi.Parameter() # TODO luigi.cfg
    outdir = luigi.Parameter(default = '')
    dist_list = luigi.Parameter() # TODO luigi.cfg? 
    org = luigi.Parameter()
    seq = luigi.Parameter()
    outfile = luigi.Parameter(default = '')
    outtype = luigi.Parameter() # TODO luigi.cfg

    def __init__(self, *args, **kwargs):
        super(metaQUAST, self).__init__(*args, **kwargs)
        self.outdir = check_make(self.outdir, f'{self.org}_{self.seq}_mQ')
        self.outfile = join(outdir, '_'.join([org, seq, 'mQ_tbl.csv']))

    def output(self):
        return luigi.LocalTarget(self.outfile)

    def run(self): 
        mQ_cmd = ['python', '/home/lam4003/bin/quast-master/metaquast.py']
        for i in dist_list:
            mQ_cmd.append(f'{self.org}_{self.seq}_{i}/scaffolds.fasta')
        mQ_cmd.extend(['-o', self.outdir, '-r', self.refdir])
        subprocess.run(mQ_cmd)
        abbreviate_table(self.outdir, ','.join(dist_list), self.outfile, self.outtype) # TODO


class AllAnalyses(luigi.WrapperTask):
    date = luigi.DateParameter(default=datetime.date.today()) # TODO change?

    def requires(self):
        yield metaQUAST(self.date)
        # yield DeconvStats(self.date) # TODO fill out