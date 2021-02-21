import luigi
from os import makedirs
from os.path import join as djoin
from os.path import exists
import subprocess
import shutil

from assembly_support import (
    barcode_sorter,
    complete_reads
)

from evaluate_support import (
    fastq_to_table,
    generate_summaries
)


# spack load gcc@6.3.0
# PYTHONPATH='.' luigi --module ariadne_pipeline de_Novo_Assembly --local-scheduler


def check_make(curr_dir, sub):
    outdir = djoin(curr_dir, sub)
    if not exists(outdir):
        makedirs(outdir)
    return outdir


class Single_Sort_FastQ(luigi.Task):
    in_dir = luigi.Parameter()
    prefix = luigi.Parameter()
    direction = luigi.Parameter()
    out_dir = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(djoin(self.out_dir, '.'.join([self.prefix + '_bsort', self.direction, 'fastq'])))

    def run(self):
        barcode_sorter(djoin(self.in_dir, '.'.join([self.prefix, self.direction, 'fastq'])), self.out_dir, False)


class Sort_FastQs_by_Barcode(luigi.Task):
    """Sort raw barcoded reads by barcode."""

    def requires(self):
        return [Single_Sort_FastQ(in_dir=gp.read_dir, prefix=gp.prefix, direction='R1', out_dir=gp.read_dir),
                Single_Sort_FastQ(in_dir=gp.read_dir, prefix=gp.prefix, direction='R2', out_dir=gp.read_dir)]

    def output(self):
        return self.input()


class cloudSPAdes_init(luigi.Task):
    """Run cloudSPAdes on the sorted read clouds."""

    def requires(self):
        return Sort_FastQs_by_Barcode()

    def output(self):
        return luigi.LocalTarget(djoin(gp.cs_init_dir, 'scaffolds.fasta'))
        #[luigi.LocalTarget(djoin(gp.cs_init_dir, 'K55', gp.search_dist + '.R1.fastq')), luigi.LocalTarget(djoin(gp.cs_init_dir, 'K55', gp.search_dist + '.R2.fastq'))]

    def run(self):
        retcode = subprocess.run(['/home/lam4003/bin/spades/assembler/spades.py', '--meta', '--only-assembler', '--gemcode1-1',
                        self.input()[0].path, '--gemcode1-2', self.input()[1].path, '--search-distance', gp.search_dist,
                        '--size-cutoff', '6', '-t', gp.num_threads, '-m', gp.memory, '-o', gp.cs_init_dir])
        spades_mem = int(gp.memory) + 100
        while retcode.returncode:
            retcode = subprocess.run(['/home/lam4003/bin/spades/assembler/spades.py', '--restart-from', 'k55', '-m', str(spades_mem), '-o', gp.cs_init_dir])
            spades_mem += 100


class Single_Complete_FastQ(luigi.Task):
    """Complete the Ariadne-deconvolved read dataset i.e.: add missing reads from the original sorted FastQ."""
    direction = luigi.Parameter()

    def requires(self):
        return cloudSPAdes_init()

    def output(self):
        return luigi.LocalTarget(
            djoin(gp.ariadne_dir, '.'.join([gp.ariadne_prefix, self.direction, 'fastq'])))

    def run(self):
        complete_reads(djoin(gp.read_dir, '.'.join([gp.prefix + '_bsort', self.direction, 'fastq'])),
                       djoin(gp.cs_init_dir, 'K55', '.'.join([gp.search_dist, self.direction, 'fastq'])),
                       djoin(gp.ariadne_dir, '.'.join([gp.search_dist + '_full', self.direction, 'fastq'])))
        yield Single_Sort_FastQ(in_dir = gp.ariadne_dir, prefix = gp.search_dist + '_full', direction = self.direction,
                          out_dir = gp.ariadne_dir)


class Complete_Ariadne_FastQs(luigi.Task):
    """Individually complete the Ariadne-deconvolved FastQs, then sort in order of the first FastQ."""

    def requires(self):
        finished_fastqs = []
        for i in ['R1', 'R2']:
            finished_fastqs.append(Single_Complete_FastQ(direction=i))
        return finished_fastqs

    def output(self):
        return [luigi.LocalTarget(
            djoin(gp.ariadne_dir, gp.ariadne_prefix + '.R1.fastq')), luigi.LocalTarget(
            djoin(gp.ariadne_dir, gp.ariadne_prefix + '.R2.fastq'))]


class cloudSPAdes_final(luigi.Task):
    """Run cloudSPAdes on the completed Ariadne-deconvolved read clouds."""

    def requires(self):
        return cloudSPAdes_init()

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, gp.search_dist + '.scaffolds.fasta'))

    def run(self):
        fastq_prefix = djoin(gp.cs_init_dir, 'K55', gp.search_dist)
        retcode = subprocess.run(['/home/lam4003/bin/spades/assembler/spades.py', '--meta', '--only-assembler', '--gemcode1-1',
                        fastq_prefix + '.R1.fastq', '--gemcode1-2', fastq_prefix + '.R2.fastq', '--search-distance', '0',
                        '--size-cutoff', '6', '-t', gp.num_threads, '-m', gp.memory, '-o', gp.cs_final_dir])
                        # '-k', '55', '--assembly-graph', djoin(gp.cs_init_dir, 'assembly_graph.fastg')])
        spades_mem = int(gp.memory) + 100
        while retcode.returncode:
            retcode = subprocess.run(['/home/lam4003/bin/spades/assembler/spades.py', '--restart-from', 'k55', '-m', str(spades_mem), '-o', gp.cs_final_dir])
            spades_mem += 100
        shutil.copy(djoin(gp.cs_final_dir, 'scaffolds.fasta'), self.output().path)


class Single_FastQ_to_Table(luigi.Task):
    sorted_fastq = luigi.Parameter()
    direction = luigi.Parameter()

    def requires(self):
        return Complete_Ariadne_FastQs()

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, '.'.join([gp.search_dist, self.direction, 'csv'])))

    def run(self): 
        fastq_to_table(self.sorted_fastq, djoin(gp.enhd_cld_dir, '.'.join(['bsort', self.direction, 'csv'])), gp.analyses_dir)


class Summarize_FastQ_Statistics(luigi.Task):
    """Generate size, purity, and entropy summary statistics from enhanced reads."""

    def requires(self):
        # Match predicted read clouds to actual (reference sequence) read clouds
        fastq_tbls = []
        for i, j in enumerate(['R1','R2']):
            fastq_tbls.append(Single_FastQ_to_Table(sorted_fastq = self.input()[i].path, direction = j))
        return fastq_tbls

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, gp.search_dist + '.statistics.csv'))

    def run(self):
        # Generate read cloud quality statistics tables 
        generate_summaries(self.input()[0].path, self.input()[1].path, None, gp.analyses_dir)


class de_Novo_Assembly(luigi.WrapperTask):
    """Top-level function calling cloudSPAdes and enhancement summaries."""
    read_dir = luigi.Parameter()
    prefix = luigi.Parameter()
    search_dist = luigi.Parameter()
    master_dir = luigi.Parameter()
    num_threads = luigi.Parameter()
    memory = luigi.Parameter()

    def requires(self):
        gp.set_params(self.read_dir, self.prefix, self.search_dist, self.master_dir, self.num_threads, self.memory)
        yield cloudSPAdes_final()
        # yield Summarize_FastQ_Statistics()


class Global_Parameters:
    """Store parameters, prefixes, and target paths. Adapted from https://stackoverflow.com/questions/51152485/luigi-global-variables."""

    def __init__(self):
        self.read_dir = None
        self.work_dir = None
        self.prefix = None
        self.search_dist = None
        self.num_threads = None
        self.memory = None
        self.cs_init_dir = None
        self.ariadne_dir = None
        self.cs_final_dir = None
        self.analyses_dir = None

    def set_params(self, read_dir, prefix, search_dist, master_dir, num_threads, memory):
        self.read_dir = read_dir
        self.work_dir = check_make(master_dir, prefix + '_' + search_dist)
        self.prefix = prefix
        self.search_dist = search_dist
        self.ariadne_prefix =  search_dist + '_full_bsort'
        self.num_threads = num_threads
        self.memory = memory
        self.cs_init_dir = check_make(self.work_dir, 'cloudSPAdes_init')
        self.ariadne_dir = check_make(self.work_dir, 'ariadne_intermediate')
        self.cs_final_dir = check_make(self.work_dir, 'cloudSPAdes_final')
        self.analyses_dir = check_make(self.work_dir, 'analyses')


gp = Global_Parameters()

if __name__ == '__main__':
    luigi.run(main_cls_task=de_Novo_Assembly)
