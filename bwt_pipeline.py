import luigi
from os import listdir, makedirs
from os.path import join as djoin
from os.path import abspath, basename, exists
import subprocess
import shutil 

from assembly_support import (
    barcode_sorter,
    bam_to_annotate,
    add_barcodes,
    subdiv_annotations,
    concat_annotations,
    fastq_enhance
)

from evaluate_support import (
    fastq_to_table,
    generate_summaries
)

# PYTHONPATH='.' luigi --module bwt_pipeline bwt_de_Novo_Assembly --local-scheduler


def check_make(curr_dir, sub):
    outdir = djoin(curr_dir, sub)
    if not exists(outdir):
        makedirs(outdir)
    return outdir


class Single_Sort_FastQ(luigi.Task):
    direction = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(djoin(gp.read_dir, '.'.join([gp.base_prefix + '_' + gp.sort_prefix, self.direction, 'fastq'])))

    def run(self):
        print(djoin(gp.read_dir, '.'.join([gp.base_prefix, self.direction, 'fastq'])))
        barcode_sorter(djoin(gp.read_dir, '.'.join([gp.base_prefix, self.direction, 'fastq'])), gp.read_dir, False)


class Sort_FastQs_by_Barcode(luigi.Task):
    """Sort raw barcoded reads by barcode."""

    def requires(self):
        return [Single_Sort_FastQ(direction = 'R1'), Single_Sort_FastQ(direction = 'R2')]            

    def output(self): 
        return self.input()


class Map_Original_Clouds(luigi.Task):
    """Map barcoded reads to the metagenome and add barcode information."""
    reference = luigi.Parameter()
    ids_to_names = luigi.Parameter()

    def requires(self):
        return Sort_FastQs_by_Barcode()

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, gp.sort_prefix + '.final.csv'))

    def run(self):
        raw_sam = djoin(gp.work_dir, gp.sort_prefix + '.sam')
        raw_bam = djoin(gp.work_dir, gp.sort_prefix + '_raw.bam')
        srt_bam = djoin(gp.work_dir, gp.sort_prefix + '.bam')

        map_step = subprocess.run(['bowtie2', '--sensitive-local', '-p', gp.num_threads, '-x', self.reference, '-1', self.input()[0].path, '-2', self.input()[1].path, '-S', raw_sam])
        bam_step = subprocess.run(['samtools', 'view', '-hb', '-f', '0', '-F', '256', raw_sam, '-o', raw_bam]) # Primary alignments only
        srt_step = subprocess.run(['samtools', 'sort', raw_bam, '-o', srt_bam])
        idx_step = subprocess.run(['samtools', 'index', srt_bam])
        bam_to_annotate(srt_bam, self.ids_to_names, gp.work_dir, '--est-fragments')
        add_barcodes(self.input()[0].path, djoin(gp.work_dir, gp.sort_prefix + '.csv'), gp.work_dir)


class Subdivide_Original_Clouds(luigi.Task):

    def requires(self): 
        return Map_Original_Clouds()

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, 'subdivide_original_clouds.out'))

    def run(self):
        subdiv_annotations(self.input().path, gp.num_chunks, gp.orig_map_dir)
        with self.output().open('w') as of:
            of.write(f'{gp.num_chunks} chunks were created')


class Single_Enhanced_Chunk(luigi.Task):
    chunk_num = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super(Single_Enhanced_Chunk, self).__init__(*args, **kwargs)
        self.chunk_file_name = '.'.join([gp.sort_prefix, str(self.chunk_num), 'csv'])

    def output(self):
        return luigi.LocalTarget(djoin(gp.enhd_cld_dir, 'R1', self.chunk_file_name)), luigi.LocalTarget(djoin(gp.enhd_cld_dir, 'R2', self.chunk_file_name))

    def run(self):
        concat_annotations(djoin(gp.orig_map_dir, self.chunk_file_name), gp.enhd_cld_dir, True)


class Generate_Enhanced_Chunks(luigi.Task):
    """Convert reference sequence (and position) information to enhanced read clouds."""

    def requires(self):
        return Subdivide_Original_Clouds()

    def output(self):
        return [luigi.LocalTarget(djoin(gp.enhd_cld_dir, '.'.join([gp.sort_prefix, 'R1', 'csv']))), luigi.LocalTarget(djoin(gp.enhd_cld_dir, '.'.join([gp.sort_prefix, 'R2', 'csv'])))]

    def run(self):
        for i in range(gp.num_chunks):
            yield Single_Enhanced_Chunk(chunk_num = i)
        for i, j in enumerate(['R1','R2']):
            with self.output()[i].open('w') as outf:
                for chunk_file in listdir(djoin(gp.enhd_cld_dir, j)):
                    with open(djoin(gp.enhd_cld_dir, j, chunk_file)) as f: 
                        for l in f:
                            outf.write(l)


class Single_Enhanced_FastQ(luigi.Task):
    enh_csv = luigi.Parameter()
    direction = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(djoin(gp.read_dir, '.'.join([gp.final_prefix, self.direction, 'fastq'])))

    def run(self): 
        fastq_enhance(djoin(gp.read_dir, '.'.join([gp.base_prefix + '_' + gp.sort_prefix, self.direction, 'fastq'])), self.enh_csv, gp.enhd_cld_dir, gp.final_prefix) 
        barcode_sorter(djoin(gp.enhd_cld_dir, '.'.join([gp.base_prefix + '_bwt', self.direction, 'fastq'])), gp.read_dir, False)


class Enhance_Original_FastQs(luigi.Task):
    """Create read files with enhanced read clouds."""

    def requires(self):
        return Generate_Enhanced_Chunks()

    def output(self):
        return [luigi.LocalTarget(djoin(gp.read_dir, gp.final_prefix + '.R1.fastq')), luigi.LocalTarget(djoin(gp.read_dir, gp.final_prefix + '.R2.fastq'))]

    def run(self): 
        for i, j in enumerate(['R1','R2']):
            yield Single_Enhanced_FastQ(enh_csv = self.input()[i].path, direction = j)


class cloudSPAdes(luigi.Task):
    """Run cloudSPAdes on the bwt-enhanced read clouds."""
    memory = luigi.Parameter()

    def requires(self):
        return Enhance_Original_FastQs()

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, 'bwt.scaffolds.fasta'))

    def run(self):
        subprocess.run(['/home/lam4003/bin/spades/assembler/spades.py', '--meta', '--only-assembler', '--gemcode1-1', self.input()[0].path, '--gemcode1-2', self.input()[1].path, '--search-distance', '0', '--size-cutoff', '6', '-t', gp.num_threads, '-m', self.memory, '-o', gp.cldspades_dir])
        shutil.copy(djoin(gp.cldspades_dir, 'scaffolds.fasta'), self.output().path)


class Single_FastQ_to_Table(luigi.Task):
    sorted_fastq = luigi.Parameter()
    direction = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, '.'.join(['bwt', self.direction, 'csv'])))

    def run(self): 
        fastq_to_table(self.sorted_fastq, djoin(gp.enhd_cld_dir, '.'.join([gp.sort_prefix, self.direction, 'csv'])), gp.analyses_dir)
        shutil.copy(basename(self.sorted_fastq).replace('fastq', 'csv'), self.output().path)


class Summarize_FastQ_Statistics(luigi.Task):
    """Generate size, purity, and entropy summary statistics from enhanced reads."""

    def requires(self):
        return Enhance_Original_FastQs()

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, 'bwt.statistics.csv'))

    def run(self):
        # Match predicted read clouds to actual (reference sequence) read clouds
        fastq_tbls = []
        for i, j in enumerate(['R1','R2']):
            yield Single_FastQ_to_Table(sorted_fastq = self.input()[i].path, direction = j)
            fastq_tbls.append(djoin(gp.analyses_dir, '.'.join(['bwt', j, 'csv'])))
        # Generate read cloud quality statistics tables 
        generate_summaries(fastq_tbls[0], fastq_tbls[1], None, gp.analyses_dir)


class bwt_de_Novo_Assembly(luigi.WrapperTask):
    """Top-level function calling cloudSPAdes and enhancement summaries."""
    read_dir = luigi.Parameter()
    prefix = luigi.Parameter()
    master_dir = luigi.Parameter()
    num_threads = luigi.Parameter()
    num_chunks = luigi.IntParameter()

    def requires(self):
        gp.set_params(self.read_dir, self.prefix, self.master_dir, self.num_threads, self.num_chunks)
        yield cloudSPAdes()
        yield Summarize_FastQ_Statistics()


class Global_Parameters:
    """Store parameters, prefixes, and target paths. Adapted from https://stackoverflow.com/questions/51152485/luigi-global-variables."""

    def __init__(self):
        self.read_dir = None
        self.work_dir = None
        self.base_prefix = None
        self.sort_prefix = 'bsort'
        self.final_prefix = None
        self.num_threads = None
        self.orig_map_dir = None
        self.enhd_cld_dir = None
        self.cldspades_dir = None
        self.analyses_dir = None
        self.num_chunks = None

    def set_params(self, read_dir, prefix, master_dir, num_threads, num_chunks):
        self.read_dir = read_dir
        self.work_dir = check_make(master_dir, prefix + '_bwt')
        self.base_prefix = prefix
        self.final_prefix = prefix + '_bwt_bsort'
        self.num_threads = num_threads
        self.orig_map_dir = check_make(self.work_dir, 'original_mapping')
        self.enhd_cld_dir = check_make(self.work_dir, 'enhanced_clouds')
        check_make(self.enhd_cld_dir, 'R1')
        check_make(self.enhd_cld_dir, 'R2')
        self.cldspades_dir = check_make(self.work_dir, 'cloudSPAdes')
        self.analyses_dir = check_make(master_dir, prefix + '_analyses')
        self.num_chunks = num_chunks


gp = Global_Parameters()

if __name__ == '__main__':
    luigi.run(main_cls_task = bwt_de_Novo_Assembly)
