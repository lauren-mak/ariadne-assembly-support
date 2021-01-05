import csv
import luigi
from os.path import basename, dirname, join, exists
from os import makedirs
import pandas as pd
from random import randint, shuffle
import subprocess
import time


# PYTHONPATH='.' luigi --module simulation Generate_Strain_Datasets --workers 10
def check_make(curr_dir, sub):
    outdir = join(curr_dir, sub)
    if not exists(curr_dir):
        makedirs(curr_dir)
    return outdir


class Single_Sort_FastQ(luigi.Task):
    # direction

    def output(self):
        # bwa_sorted.R{1,2}.fastq

    def run(self): 
        # python3 assembly_support.py barcode_sorter tests/mock5_10x.R1.fastq tests/ --cloud-sizes


class Sort_FastQs_by_Barcode(luigi.Task):
    """Sort raw barcoded reads by barcode."""

    def output(self): 
        # return {'output1' : luigi.LocalTarget('task_a_out1'),
        #         'output2' : luigi.LocalTarget('task_a_out2')}
        # bwa_sorted.R{1,2}.fastq
        pass

    def run(self):
        for i in [1,2]:
            # yield Single_Sort_FastQ(i)
            pass


class Map_Original_Clouds(luigi.Task):
    """Map barcoded reads to the metagenome and add barcode information."""
    # num_threads
    # memory

    def requires(self):
        # return Sort_by_Barcode()
        pass

    def output(self):
        # bwa_sorted.final.csv
        pass

    def run(self):
        # bwa mem -M ${REFERENCE} ${READIR}.R1.fastq ${READIR}.R2.fastq | samtools view -hb -f 0 -F 256 - | samtools sort -o ${READIR}\_psort.bam -
        # samtools index ${READIR}\_psort.bam
        # python /athena/ihlab/scratch/lam4003/scripts/testing_support.py bam_to_annotate mock5_10x_new_sorted.bam mock5_ref_ids.csv > bam_to_annotate.out 2>&1 &
        # python /athena/ihlab/scratch/lam4003/scripts/testing_support.py add_barcodes /athena/ihlab/scratch/lam4003/mock_original_reads/mock5_10x.R1.fastq mock5_10x_new_sorted.csv > bwa_new_add_barcodes.out 2>&1 &
        pass


class Subdivide_Original_Clouds(luigi.Task):
    # outdir = original_mapping

    def requires(self): 
        # return Map_Original_Clouds()
        pass

    def output():
        # return {'output1' : luigi.LocalTarget('task_a_out1'),
        # 'output2' : luigi.LocalTarget('task_a_out2')}
        # original_mapping/bwa_sorted.{i}.csv
        pass 

    def run(self):
        # python /athena/ihlab/scratch/lam4003/scripts/testing_support.py subdiv_annotations mock5_10x_new_sorted.final.csv > bwa_new_subdiv_annotations.out 2>&1 &
        pass


class Single_Enhanced_Chunk(luigi.Task):
    # chunk_num
    # direction
    # indir = original_mapping
    # outdir = enhanced_clouds

    def requires(self):
        # return Subdivide_Original_Clouds()
        pass

    def output(self):
        # enhanced_clouds/R{1,2}/bwa_sorted.{i}.R{1+2}.csv
        pass

    def run(self):
        # python /athena/ihlab/scratch/lam4003/scripts/testing_support.py concat_annotations mock5_10x_new_sorted.${i}.csv &
        pass


class Generate_Enhanced_Chunks(luigi.WrapperTask):
    """Convert reference sequence (and position) information to enhanced read clouds."""

    def requires(self):
        for i in range(100):
            # Single_Enhanced_Chunk(i)
        pass


class Single_Enhanced_FastQ(luigi.Task):
    # direction

    def output(self):
        # bwa_final.R{1,2}.fastq

    def run(self): 
        # python /athena/ihlab/scratch/lam4003/scripts/testing_support.py fastq_enhance /athena/ihlab/scratch/lam4003/mock_original_reads/mock5_10x.R1.fastq ./ --enh_dir R1/ > bwa_enhance_R1.out 2>&1 &


class Enhance_Original_FastQs(luigi.Task):
    """Create read files with enhanced read clouds."""

    def requires(self):
        # Generate_Enhanced_Chunks()
        pass

    def output(self):
        # bwa_final.R{1,2}.fastq (each of the above is a separate task)

    def run(self): 
        for i in [1,2]:
            # yield Single_Enhanced_FastQ(i)
            pass


class cloudSPAdes(luigi.Task):
    """Run cloudSPAdes on the bwa-enhanced read clouds."""
    # cs_out_dir
    # num_threads
    # memory

    def requires(self):
        # Enhance_Original_FastQs()
        pass

    def output(self):
        # x_analyses/bwa.scaffolds.fasta

    def run(self):
        # /home/lam4003/bin/spades/assembler/spades.py --meta --only-assembler --gemcode1-1 ${READIR}.R1.fastq --gemcode1-2 ${READIR}.R2.fastq --search-distance ${3} --size-cutoff 6 -t 30 -m 150 -o ${OUTDIR}
        with open(ms_out, 'w') as mout: 
            subprocess.Popen(['/Users/laurenmak/Programs/msdir/ms', str(num_total_strains), '1', '-s', self.num_variants, '-t', '0.000001'], stdout=mout)
        # cp scaffolds.fasta to results_directory


# TODO add fastq_to_table() and generate_summaries() tasks for purity/entropy analyses


class Targets:
    """Store parameters, prefixes, and target paths."""

    def __init__(self):
        """Initialize paths and variables."""
        self.read_dir = None
        self.work_dir = None
        self.prefix = None
        self.orig_map_dir = None
        self.enhd_cld_dir = None
        self.cldspades_dir = None
        self.analyses_dir = None

    def set_read_dir(self, read_dir):
        self.read_dir = read_dir

    def set_prefix(self, prefix):
        self.prefix = prefix

    def make_work_dir(self, master_dir):
        self.work_dir = check_make(join(master_dir, self.prefix + '_bwa'))

    def update_bsort(self):
        self.prefix += '_bsort'

    def update_psort(self):
        self.prefix += '_psort'

    def set_prefix_final(self, prefix):
        self.prefix = prefix + '_final'

    def make_orig_map_dir(self):
        self.orig_map_dir = check_make(join(self.work_dir, 'original_mapping'))

    def make_enhd_cld_dir(self):
        self.enhd_cld_dir = check_make(join(self.work_dir, 'enhanced_clouds'))

    def make_cldspades_dir(self):
        self.cldspades_dir = check_make(join(self.work_dir, 'cloudSPAdes'))

    def make_analyses_dir(self, master_dir):
        self.analyses_dir = check_make(join(self.master_dir, self.prefix + '_analysis'))


params = Targets()

if __name__ == '__main__':
    luigi.run(main_cls_task = cloudSPAdes)
