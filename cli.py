
import click
import assembly_support
import evaluate_support 

@click.group()
def main():
    pass


"""
General FastQ Preprocessing Functions
"""

@main.command('barcode_sorter')
@click.argument('in_fq') # FastQ file with original read clouds
@click.argument('outdir') 
@click.option('--cloud-sizes', flag_value = True)
def barcode_sorter(in_fq, outdir, cloud_sizes):
    """Sorts FastQ reads based on their barcode and cloud number. Tags barcoded reads with '-1' if not already there. Must be used at the end of any enhancement workflow and/or before usage in cloudSPAdes."""
    assembly_support.barcode_sorter(in_fq, outdir, cloud_sizes)


"""
Support Functions for Ariadne
"""


@main.command('complete_reads')
@click.argument('original_fq') # FastQ file with original read clouds
@click.argument('enhanced_fq') # FastQ file with enhanced read clouds
@click.argument('full_fq') # Output FastQ file with missing unbarcoded reads
def complete_reads(original_fq, enhanced_fq, full_fq):
    """Adds reads missing from the enhanced FastQs based on the total set of reads in the original FastQs."""
    assembly_support.complete_reads(original_fq, enhanced_fq, full_fq)


"""
FastQ Preprocessing Functions for Other Tools
"""


@main.command('bam_to_annotate')
@click.argument('bam') # Completed BAM file
@click.argument('id_csv') # Acccesion ID to reference sequence name  
@click.argument('outdir') 
@click.option('--est-fragments', flag_value=True)
def bam_to_annotate(bam, id_csv, outdir, est_fragments):
    """Generates read cloud enhancement information from BAM files. Output is [read, direction, reference name, (insert start position), (barcode)]."""
    assembly_support.bam_to_annotate(bam, id_csv, outdir, est_fragments)


@main.command('add_barcodes')
@click.argument('fastq')
@click.argument('annot_csv')
@click.argument('outdir') 
def add_barcodes(fastq, annot_csv, outdir):
    """Re-adds FastQ barcodes to bwa-based annotation file generated by bam_to_annotate. Output is [read, direction, reference name, (insert start position), barcode]. bwa-specific."""
    assembly_support.add_barcodes(fastq, annot_csv, outdir)


@main.command('subdiv_annotations')
@click.argument('infile') # Aggregated read information 
@click.argument('num_chunks')
@click.argument('outdir') 
def subdiv_annotations(infile, num_chunks, outdir):
    """Subsets total set of read clouds into smaller sets of read clouds so that the next step (concat_annotations) doesn't take forever. bwa-specific."""
    assembly_support.subdiv_annotations(infile, num_chunks, outdir)


@main.command('concat_annotations')
@click.argument('infile') # Aggregated read information 
@click.argument('outdir') 
@click.option('--est-fragments', flag_value=True)
def concat_annotations(infile, outdir, est_fragments):
    """Aggregates reference mapping information in each read cloud and generates enhanced groupings. Output is [read, enh_num, (ref-enh_num)]."""
    assembly_support.concat_annotations(infile, outdir, est_fragments)


@main.command('fastq_enhance')
@click.argument('in_fq')
@click.argument('enh_csv')
@click.argument('outdir')
@click.argument('tool')
def fastq_enhance(in_fq, enh_csv, outdir, tool):
    """Replaces original cloud number with tool-enhanced number."""
    assembly_support.fastq_enhance(in_fq, enh_csv, outdir, tool)


@main.command('fastq_to_table')
@click.argument('fastq') # FastQ file with original read clouds
@click.argument('mapping_tbl') # Table of read name,barcode,direction,ref_seq
@click.argument('outdir') # Analysis output directory
def fastq_to_table(fastq, mapping_tbl, outdir):
    """Generates tables pairing enhanced cloud numbers with the actual reference sequence. Output is [barcode, enhanced_num, reference name]"""
    evaluate_support.fastq_to_table(fastq, mapping_tbl, outdir)


@main.command('generate_summaries')
@click.argument('forward_tbl') # Forward FastQ table
@click.argument('reverse_tbl') # Reverse FastQ table 
@click.argument('outdir') # Analysis output directory
@click.option('-i', '--id_csv') # Acccesion ID to reference sequence name  
def generate_summaries(forward_tbl, reverse_tbl, outdir, id_csv): 
    """Calculate summary statistics for each and across all enhanced read clouds using the matched enhanced-actual information above. Output are [barcode, enhanced_num, size, purity, entropy] and [barcode, enhanced_num, species list]"""
    evaluate_support.generate_summaries(forward_tbl, reverse_tbl, outdir, id_csv)


@main.command('evaluate_clouds')
@click.option('--distances', '-d', help='Ariadne search distances to compare') # e.g.: 0,10000,15000,20000
@click.option('--prefixes', '-p', help='CSV prefixes') # e.g.: mock5_10x,5000_full,10000_full
@click.argument('outdir') 
def evaluate_clouds(distances, prefixes, outdir):
    """Make summary graphs of the summary information table from generate_summaries()."""
    evaluate_support.evaluate_clouds(distances, prefixes, outdir)


@main.command('pairwise_graph_align')
@click.argument('fastg') # SPAdes-format assembly graph
@click.argument('fasta') # FastA version of reads (names and sequences) 
@click.argument('outdir') # Analysis output directory
@click.option('--depth', '-d', type = int, default = 6) # Number of edges deep to search
@click.option('--fragment-mode', '-f', flag_value = True)
@click.option('--aligned-only', '-a', flag_value = True)
def pairwise_graph_align(fastg, fasta, outdir, depth, fragment_mode, aligned_only):
    """Identify read cloud assembly graph alignments and pairwise differences."""
    evaluate_support.pairwise_graph_align(fastg, fasta, outdir, depth, fragment_mode, aligned_only)


@main.command('make_barcode_whitelist')
@click.argument('in_fq') # FastQ file with original read clouds
@click.argument('outdir') 
def make_barcode_whitelist(in_fq, outdir):
    """Extract linked-read barcodes into a custom whitelist from FastQ. For the EMA tool."""
    assembly_support.make_barcode_whitelist(in_fq, outdir)


@main.command('interleave_fastqs')
@click.argument('in_fq_prefix') # Read directory plus sample name
@click.argument('outdir') 
@click.option('--ema', '-ema', flag_value = True)
def interleave_fastqs(in_fq_prefix, outdir, ema):
    """Interleave paired-end FastQs."""
    assembly_support.interleave_fastqs(in_fq_prefix, outdir, ema)


@main.command('tidy_quast_report')
@click.argument('prefix')
@click.argument('outdir') 
def tidy_quast_report(prefix, outdir):
    evaluate_support.tidy_quast_report(prefix, outdir)


@main.command('graph_assembly_stats')
@click.argument('base_dir')
def tidy_quast_report(base_dir):
    evaluate_support.graph_assembly_stats(base_dir)


@main.command('graph_cloud_stats')
@click.argument('base_dir')
@click.argument('num_clouds')
@click.argument('param')
@click.option('--scale', '-s', flag_value = True)
def tidy_quast_report(base_dir, num_clouds, param, scale):
    evaluate_support.graph_cloud_stats(base_dir, num_clouds, param, scale)


if __name__ == '__main__':
    main()


