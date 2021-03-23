import luigi
from os import listdir, makedirs
from os.path import join as djoin
from os.path import exists
import subprocess
import shutil 

from evaluate_support import (
    evaluate_clouds,
    tidy_quast_report
)


# PYTHONPATH='.' luigi --module evaluate_pipeline de_Novo_Assembly --local-scheduler
# ./run_pipeline.sh evaluate mock5_10x.cfg 10


def check_make(curr_dir, sub):
    outdir = djoin(curr_dir, sub)
    if not exists(outdir):
        makedirs(outdir)
    return outdir


class metaQUAST(luigi.Task):
    reference = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(djoin(gp.out_dir, 'report.tsv')) # TODO Write a processor function for this

    def run(self):
        mq_dir = djoin(gp.work_dir, 'metaQUAST')
        mq_cmd = ['python', '/home/lam4003/bin/quast-master/metaquast.py', '-l', gp.deconv_labels, \
            '-o', mq_dir, '-r', self.reference]
        for fp in gp.file_prefixes:
            mq_cmd.append(djoin(gp.work_dir, fp + '.scaffolds.fasta'))
        subprocess.run(mq_cmd)
        for f in listdir(mq_dir):
            if f.endswith('.html'):
                shutil.copy(djoin(mq_dir, f), gp.out_dir)
        shutil.copytree(djoin(mq_dir, 'icarus_viewers'), djoin(gp.out_dir, 'icarus_viewers'))
        shutil.copy(djoin(mq_dir, 'combined_reference', 'report.tsv'), gp.out_dir) 


class Tidy_metaQUAST_Report(luigi.Task):

    def requires(self):
        return metaQUAST()

    def output(self):
        return luigi.LocalTarget(djoin(gp.out_dir, gp.base_prefix + '.metaQUAST.tex'))

    def run(self):
        tidy_quast_report(gp.base_prefix, gp.out_dir)


class Aggregate_Statistics(luigi.Task):

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, 'Avg_Summary_Stats.tbl')) 

    def run(self):
        deconv_labels = gp.deconv_labels.split(',')
        deconv_labels.remove('Illumina')
        deconv_no_illumina = ','.join(deconv_labels)
        file_prefixes = gp.file_prefixes
        file_prefixes.remove('illumina')
        prefixes_no_illumina = ','.join(file_prefixes)

        evaluate_clouds(deconv_no_illumina, prefixes_no_illumina, gp.work_dir)
        for f in listdir(gp.work_dir):
            if 'Purity' in f or 'Entropy' in f or 'Size' in f or 'Status' in f:
                shutil.copy(djoin(gp.work_dir, f), gp.out_dir)
        shutil.copy(djoin(gp.work_dir, 'Avg_Summary_Stats.tbl'), gp.out_dir)


class de_Novo_Assembly(luigi.WrapperTask):
    """Top-level function calling summary statistic aggregators and graphics generators."""
    prefix = luigi.Parameter()
    master_dir = luigi.Parameter()
    deconv_methods = luigi.Parameter()

    def requires(self):
        gp.set_params(self.prefix, self.master_dir, self.deconv_methods)
        yield Tidy_metaQUAST_Report()
        yield Aggregate_Statistics()


class Global_Parameters:
    """Store parameters, prefixes, and target paths. Adapted from https://stackoverflow.com/questions/51152485/luigi-global-variables."""

    def __init__(self):
        self.base_prefix = None
        self.work_dir = None
        self.out_dir = None
        self.deconv_labels = None
        self.file_prefixes = None

    def set_params(self, prefix, master_dir, deconv_methods):
        self.base_prefix = prefix
        self.work_dir = check_make(master_dir, prefix + '_analyses')
        self.out_dir = check_make(self.work_dir, 'output')
        self.deconv_labels = deconv_methods
        deconv_lst = deconv_methods.split(',')
        
        self.file_prefixes = []
        self.file_prefixes.append('illumina') if 'Illumina' in deconv_lst else None
        self.file_prefixes.append('no_deconv') if 'No_Deconv' in deconv_lst else None
        self.file_prefixes.append(prefix + '_spc_bsort') if 'Species' in deconv_lst else None
        self.file_prefixes.append(prefix + '_frg_bsort') if 'Fragments' in deconv_lst else None
        self.file_prefixes.append(prefix + '_ema_bsort') if 'EMA' in deconv_lst else None
        self.file_prefixes.append(prefix + '_mnv_bsort') if 'Minerva' in deconv_lst else None
        self.file_prefixes.append('1000') if '5000' in deconv_lst else None
        self.file_prefixes.append('2000') if '5000' in deconv_lst else None
        self.file_prefixes.append('4000') if '5000' in deconv_lst else None
        self.file_prefixes.append('5000') if '5000' in deconv_lst else None
        self.file_prefixes.append('10000') if '10000' in deconv_lst else None
        self.file_prefixes.append('15000') if '15000' in deconv_lst else None
        self.file_prefixes.append('20000') if '20000' in deconv_lst else None


gp = Global_Parameters()

if __name__ == '__main__':
    luigi.run(main_cls_task = de_Novo_Assembly)
