import luigi
from os.path import basename, exists, join, makedirs
import subprocess

from assembly_support import (
    barcode_sorter,
    bam_to_annotate,
    subdiv_annotations,
    concat_annotations,
    fastq_enhance
)

from evaluate_support import (
    fastq_to_table,
    generate_summaries
)

# PYTHONPATH='.' luigi --module simulation Generate_Strain_Datasets --workers 10


def check_make(curr_dir, sub):
    outdir = join(curr_dir, sub)
    if not exists(outdir):
        makedirs(outdir)
    return outdir


class Single_Sort_FastQ(luigi.Task):
    direction = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(join(gp.read_dir, '.'.join(gp.sort_prefix, self.direction, 'fastq')))

    def run(self):
        barcode_sorter(join(gp.read_dir, '.'.join(gp.base_prefix, self.direction, 'fastq')), gp.read_dir)


class Sort_FastQs_by_Barcode(luigi.Task):
    """Sort raw barcoded reads by barcode."""

    def run(self):
        return [Single_Sort_FastQ(direction = 'R1'), Single_Sort_FastQ(direction = 'R2')]            

    def output(self): 
        return self.input()


class Map_Original_Clouds(luigi.Task):
    """Map barcoded reads to the metagenome and add barcode information."""
    reference = luigi.Parameter()

    def requires(self):
        return Sort_by_Barcode()

    def output(self):
        return luigi.LocalTarget(join(gp.work_dir, gp.sort_prefix + '.final.csv'))

    def run(self):
        bam = join(gp.work_dir, gp.sort_prefix + '.bam')
        subprocess.run(['bwa', 'mem', '-t', gp.num_threads, self.reference, self.input()[0].path, self.input()[1].path, '|', 'samtools', 'view', '-hb', '-f', 0, '-F', 256, '-', '|', 'samtools', 'sort', '-o', bam], shell = True)
        # bwa mem -M ${REFERENCE} ${READIR}.R1.fastq ${READIR}.R2.fastq | samtools view -hb -f 0 -F 256 - | samtools sort -o ${READIR}\_psort.bam -
        subprocess.run(['samtools', 'index', bam])
        bam_to_annotate(bam, gp.ids_to_names)
        add_barcodes(self.input()[0].path, join(gp.work_dir, gp.sort_prefix + '.csv'))


class Subdivide_Original_Clouds(luigi.Task):

    def requires(self): 
        return Map_Original_Clouds()

    def output():
        chunk_lst = []
        for i in range(100):
            chunk_lst.append(luigi.LocalTarget(join(gp.orig_map_dir, '.'.join(gp.sort_prefix, i,'.csv'))))
        return chunk_lst

    def run(self):
        subdiv_annotations(self.input().path, gp.orig_map_dir)


class Single_Enhanced_Chunk(luigi.Task):
    chunk_num = luigi.IntParameter()

    def requires(self):
        return Subdivide_Original_Clouds() # TODO does this method of calling a dependency work?

    def output(self):
        chunk_file_name = '.'.join(gp.sort_prefix, self.chunk_num, 'csv')
        return [luigi.LocalTarget(join(gp.enhd_cld_dir, 'R1', chunk_file_name)), luigi.LocalTarget(join(gp.enhd_cld_dir, 'R2', chunk_file_name))]

    def run(self):
        concat_annotations(self.input()[self.chunk_num].path, gp.enhd_cld_dir)


class Generate_Enhanced_Chunks(luigi.WrapperTask):
    """Convert reference sequence (and position) information to enhanced read clouds."""

    def requires(self):
        for i in range(100):
            Single_Enhanced_Chunk(i)


class Single_Enhanced_FastQ(luigi.Task):
    direction = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(join(gp.read_dir, '.'.join(gp.final_prefix, self.direction, 'fastq')))

    def run(self): 
        fastq_enhance(join(gp.read_dir, '.'.join(gp.sort_prefix, self.direction, 'fastq')), gp.read_dir, join(gp.enhd_cld_dir, self.direction), None) 
        # TODO Does this work to call external function with options? 


class Enhance_Original_FastQs(luigi.Task):
    """Create read files with enhanced read clouds."""

    def requires(self):
        return Generate_Enhanced_Chunks()

    def output(self):
        return [luigi.LocalTarget(join(gp.read_dir, gp.final_prefix + '.R1.fastq')), luigi.LocalTarget(join(gp.read_dir, gp.final_prefix + '.R2.fastq'))]

    def run(self): 
        for i in ['R1','R2']:
            yield Single_Enhanced_FastQ(direction = i)


class cloudSPAdes(luigi.Task):
    """Run cloudSPAdes on the bwa-enhanced read clouds."""

    def requires(self):
        return Enhance_Original_FastQs()

    def output(self):
        return luigi.LocalTarget(join(gp.analyses_dir, 'bwa.scaffolds.fasta'))

    def run(self):
        subprocess.run(['/home/lam4003/bin/spades/assembler/spades.py', '--meta', '--only-assembler', '--gemcode1-1', self.input()[0].path, '--gemcode1-2', self.input()[1].path, '--search-distance ', '0', '--size-cutoff', '6', '-t', gp.num_threads, '-m', gp.memory, '-o', gp.cldspades_dir])
        shutil.copy(join(gp.cldspades_dir, 'scaffolds.fasta'), self.output().path)


class Single_FastQ_to_Table(luigi.Task):
    sorted_fastq = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(join(gp.analyses_dir, basename(sorted_fastq).replace('fastq', 'csv')))

    def run(self): 
        fastq_to_table(self.sorted_fastq, join(gp.work_dir, gp.sort_prefix + '.final.csv'), gp.analyses_dir)


class Summarize_FastQ_Statistics(luigi.Task):
    """Generate size, purity, and entropy summary statistics from enhanced reads."""

    def requires(self):
        return Enhance_Original_FastQs()

    def output(self):
        return luigi.LocalTarget(join(gp.analyses_dir, gp.final_prefix + '.statistics.csv'))

    def run(self):
        # Match predicted read clouds to actual (reference sequence) read clouds
        fastq_tbls = []
        for i in self.input():
            fastq_tbls.append(Single_FastQ_to_Table(sorted_fastq = i.path))
        # Generate read cloud quality statistics tables 
        generate_summaries(fastq_tbls[0].path, fastq_tbls[1].path, gp.ids_to_names, gp.analyses_dir)


class bwa_de_Novo_Assembly(luigi.WrapperTask):
    """Top-level function calling cloudSPAdes and enhancement summaries."""
    read_dir = luigi.Parameter()
    prefix = luigi.Parameter()
    master_dir = luigi.Parameter()
    num_threads = luigi.IntParameter()
    memory = luigi.IntParameter()
    ids_to_names = luigi.Parameter()

    def requires(self):
        gp.set_params(self.read_dir, self.prefix, self.master_dir, self.num_threads, self.memory, self.ids_to_names)
        yield cloudSPAdes()
        yield Summarize_FastQ_Statistics()


class Global_Parameters:
    """Store parameters, prefixes, and target paths. Adapted from https://stackoverflow.com/questions/51152485/luigi-global-variables."""

    def __init__(self):
        self.read_dir = None
        self.work_dir = None
        self.base_prefix = None
        self.sort_prefix = None
        self.final_prefix = None
        self.num_threads = None
        self.memory = None
        self.orig_map_dir = None
        self.enhd_cld_dir = None
        self.cldspades_dir = None
        self.analyses_dir = None
        self.ids_to_names = None

    def set_params(self, read_dir, prefix, master_dir, num_threads, memory, ids_to_names):
        self.read_dir = read_dir
        self.work_dir = check_make(join(master_dir, self.prefix + '_bwa'))
        self.base_prefix = prefix
        self.sort_prefix = prefix + '_bsort'
        self.final_prefix = prefix + '_final'
        self.num_threads = num_threads
        self.memory = memory
        self.orig_map_dir = check_make(join(self.work_dir, 'original_mapping'))
        self.enhd_cld_dir = check_make(join(self.work_dir, 'enhanced_clouds'))
        check_make(self.enhd_cld_dir, 'R1')
        check_make(self.enhd_cld_dir, 'R2')
        self.cldspades_dir = check_make(join(self.work_dir, 'cloudSPAdes'))
        self.analyses_dir = check_make(join(self.master_dir, self.prefix + '_analyses'))
        self.ids_to_names = ids_to_names


gp = Global_Parameters()

if __name__ == '__main__':
    luigi.run(main_cls_task = bwa_de_Novo_Assembly)
