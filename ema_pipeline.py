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
    enhance_mode = luigi.BoolParameter(parsing = luigi.BoolParameter.EXPLICIT_PARSING)

    def output(self):
        return luigi.LocalTarget(djoin(gp.read_dir, '.'.join([self.prefix + '_bsort', self.direction, 'fastq'])))

    def run(self):
        if self.enhance_mode:
            yield Single_Enhanced_FastQ(enh_csv = djoin(gp.enhd_cld_dir, '.'.join(['bsort', self.direction, 'csv'])), direction = self.direction)
        barcode_sorter(djoin(self.in_dir, '.'.join([self.prefix, self.direction, 'fastq'])), gp.read_dir, False) 


class Interleave_FastQs(luigi.Task):
    """Interleave paired-end FastQ files."""
    # Function based on https://github.com/ekg/interleave-fastq/blob/62e09673b18b0d5b0204b5342046be6131937d60/interleave-fastq

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, '.'.join([gp.sort_prefix, 'fastq'])))

    def run(self):
        read_prefix = djoin(gp.read_dir, gp.sort_prefix)
        with open(djoin(gp.read_dir, gp.sort_prefix) + '.R1.fastq') as fw_fq, \
             open(djoin(gp.read_dir, gp.sort_prefix) + '.R2.fastq') as rv_fq, \
             self.output().open('w') as out_fq: 
            while True:
                line = fw_fq.readline()
                if line.strip() == "":
                    break
                out_fq.write(line)
                for i in range(3):
                    out_fq.write(fw_fq.readline())
                for i in range(4):
                    out_fq.write(rv_fq.readline())


class Make_EMA_Counts(luigi.Task):
    """Generate preliminary barcode counts to bin reads."""

    def __init__(self, *args, **kwargs):
        super(Make_EMA_Counts, self).__init__(*args, **kwargs)
        self.out_prefix = djoin(gp.work_dir, gp.sort_prefix)

    def requires(self):
        return Interleave_FastQs()

    def output(self):
        return luigi.LocalTarget(self.out_prefix + '.ema-ncnt')

    def run(self):
        subprocess.run(['/home/lam4003/bin/ema/ema', 'count', '-w', gp.whitelist, '-o', self.out_prefix], \
            input = djoin(gp.work_dir, '.'.join([gp.sort_prefix, 'fastq'])).encode('utf-8'))


class Generate_EMA_Bins(luigi.Task):
    """Generate bins of reads for EMA-based alignment."""

    def requires(self):
        return Make_EMA_Counts()

    def output(self):
        return luigi.LocalTarget(djoin(gp.bin_dir, 'ema-bin-499'))

    def run(self):
        subprocess.run(['/home/lam4003/bin/ema/ema', 'preproc', '-w', gp.whitelist, '-n', gp.num_chunks, '-o', gp.bin_dir, \
            djoin(gp.work_dir, gp.sort_prefix + '.ema-ncnt')], \
            input = djoin(gp.work_dir, '.'.join([gp.sort_prefix, 'fastq'])).encode('utf-8'))


class Single_Bin_BAM(luigi.Task):
    """Generate best alignments to the reference genomes based on barcodes using EMA."""
    bin_num = luigi.Parameter()
    reference = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super(Single_Bin_BAM, self).__init__(*args, **kwargs)
        bin_zeroes = str('0' * (3 - len(self.bin_num)))
        self.out_prefix = djoin(gp.aln_dir, 'ema-bin-' + bin_zeroes + self.bin_num)

    def requires(self):
        return Generate_EMA_Bins()

    def output(self):
        return luigi.LocalTarget(self.out_prefix + '.bam') 

    def run(self):
        subprocess.run(['/home/lam4003/bin/ema/ema', 'align', '-d', '-r', self.reference, '-s', \
            djoin(djoin(gp.bin_dir, 'ema-bin-' + bin_zeroes + self.bin_num)), '-o', self.out_prefix + '.sam'])
        subprocess.run(['samtools', 'view', '-hb', '-f', '0', '-F', '256', self.out_prefix + '.sam', '-o', self.output().path])


class Generate_Bin_BAMs(luigi.Task):

    def requires(self):
        bin_bam_paths = []
        for i in range(int(gp.num_chunks)): 
            bin_bam_paths.append(Single_Bin_BAM(bin_num = str(i)))
        return bin_bam_paths

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, 'ema_bams.txt'))

    def run(self):
        with self.output().open('w') as outf:
            for i in self.input():
                outf.write(i.path + '\n')


class Merge_Bin_BAMs(luigi.Task):
    """Merge BAMs generated from each bin of reads aligned by EMA."""

    def requires(self):
        return Generate_Bin_BAMs()

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, 'ema.bam'))

    def run(self):
        merge_cmd = ['samtools', 'merge', self.output().path]
        with self.input().open('r') as bam_lst:
            for line in bam_lst:
                merge_cmd.append(line.strip())
        subprocess.run(merge_cmd)


class Map_EMA_Clouds(luigi.Task):
    """Map barcoded reads to the metagenome and add barcode information."""
    ids_to_names = luigi.Parameter()

    def requires(self):
        return Merge_Bin_BAMs() 

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, 'ema_psort.csv'))

    def run(self):
        srt_bam = djoin(gp.work_dir, 'ema_psort.bam')
        srt_step = subprocess.run(['samtools', 'sort', self.input().path, '-o', srt_bam])
        idx_step = subprocess.run(['samtools', 'index', srt_bam])
        bam_to_annotate(srt_bam, self.ids_to_names, gp.work_dir, False)


# class Map_Original_Clouds(luigi.Task):
#     """Map barcoded reads to the metagenome and add barcode information."""

#     def requires(self):
#         return Generate_Bin_BAMs()

#     def output(self):
#         return luigi.LocalTarget(djoin(gp.work_dir, 'ema-nobc_psort.final.csv'))

#     def run(self):
#         raw_sam = djoin(gp.aln_dir, 'ema-nobc.sam')
#         raw_bam = djoin(gp.work_dir, 'ema-nobc.bam')
#         srt_bam = djoin(gp.work_dir, 'ema-nobc_psort.bam')

#         map_step = subprocess.run(['bowtie2', '--sensitive-local', '-p', gp.num_threads, '-x', gp.reference, '--interleaved', djoin(gp.bin_dir, 'ema-nobc'), '-S', raw_sam])
#         bam_step = subprocess.run(['samtools', 'view', '-hb', '-f', '0', '-F', '256', raw_sam, '-o', raw_bam]) # Primary alignments only
#         srt_step = subprocess.run(['samtools', 'sort', raw_bam, '-o', srt_bam])
#         idx_step = subprocess.run(['samtools', 'index', srt_bam])
#         bam_to_annotate(srt_bam, gp.ids_to_names, gp.work_dir, False)
#         add_barcodes(djoin(gp.read_dir, gp.sort_prefix) + '.R1.fastq', djoin(gp.work_dir, 'ema-nobc_psort.csv'), gp.work_dir)


# def py_concat(f1, f2, fout):
#     data = d2 = ''
#     with open(f1) as r1: 
#         data = f1.read() 
#     with open(f2) as r2: 
#         d2 = f2.read() 
#     data += ('\n' + d2)
#     with open (fout) as outw: 
#         outw.write(data)


class Subdivide_Original_Clouds(luigi.Task):
    """Split EMA-based annotations into multiple chunks for faster cloud generation/FastQ processing."""

    def requires(self): 
        return Map_EMA_Clouds() # , Map_Original_Clouds()]

    def output(self):
        return luigi.LocalTarget(djoin(gp.work_dir, 'subdivide_original_clouds.out'))

    def run(self):
        # py_concat(djoin(gp.work_dir, 'ema_psort.csv'), djoin(gp.work_dir, 'ema-nobc_psort.final.csv'), djoin(gp.work_dir, 'ema_full.csv'))
        subdiv_annotations(self.input().path, gp.num_chunks, gp.orig_map_dir)
        with self.output().open('w') as of:
            of.write(f'{gp.num_chunks} chunks were created')


class Single_Enhanced_Chunk(luigi.Task):
    chunk_num = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super(Single_Enhanced_Chunk, self).__init__(*args, **kwargs)
        self.chunk_file_name = '.'.join(['bsort', str(self.chunk_num), 'csv'])

    def requires(self):
        return Subdivide_Original_Clouds()

    def output(self):
        return luigi.LocalTarget(djoin(gp.enhd_cld_dir, 'R1', self.chunk_file_name)), luigi.LocalTarget(djoin(gp.enhd_cld_dir, 'R2', self.chunk_file_name))

    def run(self):
        concat_annotations(djoin(gp.orig_map_dir, self.chunk_file_name), gp.enhd_cld_dir, False)


class Generate_Enhanced_Chunks(luigi.Task):
    """Convert reference sequence (and position) information to enhanced read clouds."""

    def __init__(self, *args, **kwargs):
        super(Generate_Enhanced_Chunks, self).__init__(*args, **kwargs)
        self.enhd_file_prefix = djoin(gp.enhd_cld_dir, 'bsort')

    def requires(self):
        for i in range(gp.num_chunks):
            yield Single_Enhanced_Chunk(chunk_num = i)

    def output(self):
        return [luigi.LocalTarget(self.enhd_file_prefix + '.R1.csv'), luigi.LocalTarget(self.enhd_file_prefix + '.R2.csv')]

    def run(self):
        for i, j in enumerate(['R1','R2']):
            enhd_file_path = '.'.join([self.enhd_file_prefix, j, 'csv'])
            with open(enhd_file_path, 'w') as outf:
                for chunk_file in listdir(djoin(gp.enhd_cld_dir, j)):
                    with open(djoin(gp.enhd_cld_dir, j, chunk_file)) as f: 
                        for l in f:
                            outf.write(l)


class Single_Enhanced_FastQ(luigi.Task):
    enh_csv = luigi.Parameter()
    direction = luigi.Parameter()

    def requires(self):
        return Generate_Enhanced_Chunks()

    def output(self):
        return luigi.LocalTarget(djoin(gp.enhd_cld_dir, '.'.join([gp.edit_prefix, self.direction, 'fastq'])))

    def run(self): 
        fastq_enhance(djoin(gp.read_dir, '.'.join([gp.base_prefix + '_spc_bsort', self.direction, 'fastq'])), self.enh_csv, gp.enhd_cld_dir, gp.edit_prefix) 


class Enhance_Original_FastQs(luigi.Task):
    """Create read files with enhanced read clouds."""

    def requires(self):
        for i in ['R1','R2']:
            yield Single_Sort_FastQ(in_dir = gp.enhd_cld_dir, prefix = gp.edit_prefix, direction = i, enhance_mode = True)

    def output(self):
        return [luigi.LocalTarget(djoin(gp.read_dir, gp.final_prefix + '.R1.fastq')), luigi.LocalTarget(djoin(gp.read_dir, gp.final_prefix + '.R2.fastq'))]


class cloudSPAdes(luigi.Task):
    """Run cloudSPAdes on the bwt-enhanced read clouds."""
    memory = luigi.Parameter()

    def requires(self):
        return Enhance_Original_FastQs()

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, gp.final_prefix + '.scaffolds.fasta'))

    def run(self):
        subprocess.run(['/home/lam4003/bin/spades/assembler/spades.py', '--meta', '--only-assembler', '--gemcode1-1', self.input()[0].path, '--gemcode1-2', self.input()[1].path, '--search-distance', '0', '--size-cutoff', '6', '-t', gp.num_threads, '-m', self.memory, '-o', gp.cldspades_dir])
        shutil.copy(djoin(gp.cldspades_dir, 'scaffolds.fasta'), self.output().path)


class Single_FastQ_to_Table(luigi.Task):
    direction = luigi.Parameter()

    def requires(self):
        return Enhance_Original_FastQs()

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, '.'.join([gp.final_prefix, self.direction, 'csv'])))

    def run(self): 
        fastq_to_table(djoin(gp.read_dir, '.'.join([gp.final_prefix, self.direction, 'fastq'])), 
                       djoin(gp.analyses_dir, '.'.join([gp.gstd_prefix, self.direction, 'csv'])), gp.analyses_dir)


class Summarize_FastQ_Statistics(luigi.Task):
    """Generate size, purity, and entropy summary statistics from enhanced reads."""

    def requires(self):
        # Match predicted read clouds to actual (reference sequence) read clouds
        fastq_tbls = []
        for i in ['R1','R2']:
            fastq_tbls.append(Single_FastQ_to_Table(direction = i))
        return fastq_tbls

    def output(self):
        return luigi.LocalTarget(djoin(gp.analyses_dir, gp.final_prefix + '.statistics.csv'))

    def run(self):
        # Generate read cloud quality statistics tables 
        generate_summaries(self.input()[0].path, self.input()[1].path, gp.analyses_dir, None)


class de_Novo_Assembly(luigi.WrapperTask):
    """Top-level function calling cloudSPAdes and enhancement summaries."""
    read_dir = luigi.Parameter()
    prefix = luigi.Parameter()
    master_dir = luigi.Parameter()
    gold_standard = luigi.Parameter()
    whitelist = luigi.Parameter()
    num_threads = luigi.Parameter()
    num_chunks = luigi.IntParameter()

    def requires(self):
        gp.set_params(self.read_dir, self.prefix, self.master_dir, self.gold_standard, self.whitelist, self.num_threads, self.num_chunks)
        yield cloudSPAdes()
        yield Summarize_FastQ_Statistics()


class Global_Parameters:
    """Store parameters, prefixes, and target paths. Adapted from https://stackoverflow.com/questions/51152485/luigi-global-variables."""

    def __init__(self):
        self.read_dir = None
        self.work_dir = None
        self.base_prefix = None
        self.edit_prefix = None
        self.final_prefix = None
        self.gstd_prefix = None
        self.whitelist = None
        self.num_threads = None
        self.bin_dir = None
        self.aln_dir = None
        self.orig_map_dir = None
        self.enhd_cld_dir = None
        self.cldspades_dir = None
        self.analyses_dir = None
        self.num_chunks = None

    def set_params(self, read_dir, prefix, master_dir, gold_standard, whitelist, num_threads, num_chunks):
        self.read_dir = read_dir
        self.base_prefix = prefix
        self.sort_prefix = prefix + '_bsort'
        self.edit_prefix = prefix + '_ema'
        self.work_dir = check_make(master_dir, self.edit_prefix)
        self.final_prefix = self.edit_prefix + '_bsort'
        self.gstd_prefix = gold_standard
        self.whitelist = whitelist
        self.num_threads = num_threads
        self.bin_dir = check_make(self.work_dir, 'bins')
        self.aln_dir = check_make(self.work_dir, 'alns')
        self.orig_map_dir = check_make(self.work_dir, 'original_mapping')
        self.enhd_cld_dir = check_make(self.work_dir, 'enhanced_clouds')
        check_make(self.enhd_cld_dir, 'R1')
        check_make(self.enhd_cld_dir, 'R2')
        self.cldspades_dir = check_make(self.work_dir, 'cloudSPAdes')
        self.analyses_dir = check_make(master_dir, prefix + '_analyses')
        self.num_chunks = num_chunks


gp = Global_Parameters()

if __name__ == '__main__':
    luigi.run(main_cls_task = de_Novo_Assembly)
