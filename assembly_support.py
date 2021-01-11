
from csv import writer
from datetime import datetime
from math import floor
from os import listdir
from os.path import basename, isfile, join
from sys import stderr

# FastQ/BAM data manipulation
import pysam

# Pandas and numpy for data manipulation
import numpy as np
import pandas as pd


def logger(msg):
    print(f'{datetime.now()}: {msg}', file = stderr)


def prefix(fn):
    return basename(fn).split('.')[0]


"""
General FastQ Preprocessing Functions
"""


# @D00547:847:HYHNTBCXX:1:1101:10000:21465 BX:Z:CGCCGAAGTCACGCAC-1
def barcode_sorter(in_fq, outdir, cloud_sizes):
    """Sorts FastQ reads based on their barcode and cloud number. Tags barcoded reads with '-1' if not already there. Must be used at the end of any enhancement workflow and/or before usage in cloudSPAdes."""
    bc2reads = {} # {barcode,{cloud_num,[read_lines]}}
    bc2reads['NA'] = {}
    bc2reads['NA']['NA'] = []
    curr_bc = ''
    curr_enh = ''
    with open(in_fq, 'r') as fq:
        for i, line in enumerate(fq): 
            output = line
            if (i % 4) == 0:
                bc = ''
                enh = ''
                if 'BX:Z:' in line:
                    info = line.strip()[1:].replace(' BX:Z:','-').split('-')
                    bc = info[1]
                    if len(info) == 3:
                        enh = info[2]
                        cloud_numbered = True
                    else:
                        enh = '1'
                        cloud_numbered = False
                    if bc not in bc2reads: 
                        bc2reads[bc] = {}
                    if enh not in bc2reads[bc]:
                        bc2reads[bc][enh] = []
                    curr_bc = bc
                    curr_enh = enh
                else: # This read has not been barcoded. 
                    curr_bc = 'NA'
                    curr_enh = 'NA'
                    cloud_numbered = True
                if not cloud_numbered:
                    output = line.strip() + '-1\n'
            bc2reads[curr_bc][curr_enh].append(output)
            if (i % 10000000) == 0: 
                logger(f'Finished processing line {i}')

    in_fq_parts = basename(in_fq).split('.')
    out_fq = join(outdir, '.'.join([in_fq_parts[0] + '_bsort', in_fq_parts[1], in_fq_parts[2]]))
    size_lst = []
    with open(out_fq, 'w') as of:
        for b in bc2reads: 
            for e in bc2reads[b]: 
                size_lst.append([b, e, len(bc2reads[b][e]) / 4]) # [[barcode,cloud_num,num_reads]]
                for line in bc2reads[b][e]: 
                    of.write(line)
    if cloud_sizes:
        with open(join(outdir, '.'.join([in_fq_parts[0] + '_sizes', in_fq_parts[1], 'csv'])), 'w') as sf:
            sw = writer(sf)
            sw.writerows(size_lst)


"""
Support Functions for Ariadne
"""


def get_read_names(fq_file):
    logger('Importing read names from the file ' + fq_file)
    names = {}
    with open(fq_file, 'r') as fq:
        for i, line in enumerate(fq): 
            if (i % 4) == 0:
                full_name = line.strip()
                names[full_name.split(' ')[0]] = full_name
                if (i % 10000000) == 0:
                    logger('Processed ' + str(i) + ' lines')            
    logger('There are ' + str(i) + ' lines in total in ' + fq_file)
    return names


def get_read_full(fq_file):
    logger('Importing reads from the file ' + fq_file)
    reads = {}
    with open(fq_file, 'r') as fq:
        for i, line in enumerate(fq): 
            if (i % 4) == 0:
                curr_name = line.strip().split(' ')[0]
                reads[curr_name] = []
                if (i % 10000000) == 0:
                    logger('Processed ' + str(i) + ' lines')            
            reads[curr_name].append(line.strip())
    logger('There are ' + str(i) + ' lines in total in ' + fq_file)
    return reads 


def complete_reads(original_fq, enhanced_fq, full_fq):
    """Adds reads missing from the enhanced FastQs based on the total set of reads in the original FastQs."""
    enh_names = get_read_names(enhanced_fq)
    all_reads = get_read_full(original_fq)
    num_processed = 0
    num_missing = 0
    with open(full_fq, 'w') as ff:
        for r in all_reads:
            r_info = all_reads[r]
            if r in enh_names:
                r_name = enh_names[r]
            else:
                r_name = r_info[0].replace('-1', '-0') # Indicates that the read was not deconvolved. 
                num_missing += 1
            ff.write(r_name + '\n' + r_info[1] + '\n' + r_info[2] + '\n' + r_info[3] + '\n')
            num_processed += 1
            if (num_processed % 1000000) == 0:
                logger('Processed ' + str(num_processed) + ' reads')
    logger('There are ' + str(num_missing) + ' reads from ' + original_fq + ' that were not included')


"""
FastQ Preprocessing Functions for Other Tools
"""


def decimalToBinary(n):  
    return bin(int(n)).replace("0b", "")  


def bam_to_annotate(bam, id_csv, outdir, est_fragments):
    """Generates read cloud enhancement information from BAM files. Output is [read, direction, reference name, (insert start position), (barcode)]."""
    inbam = pysam.AlignmentFile(bam, 'rb')
    id2seq = {}
    id2seq['None'] = 'None'
    with open(id_csv, 'r') as ic:
        for line in ic:
            info = line.strip().split(',')
            id2seq[info[0]] = info[1]

    # Ingest BAM file and extract mapped-to references and read-directions in pairs. 
    read_info_tbl = []
    inbam = pysam.AlignmentFile(bam, 'rb')
    for read in inbam.fetch(until_eof=True):
        read_info = [read.query_name]
        read_info.append(1 - int((decimalToBinary(read.flag))[-7])) # Direction: 0 = forward, 1 = reverse
        if read.reference_name and est_fragments:
            # Find the left-most coordinate, and round down to the nearest multiple of 100,000, the estimated largest fragment size. 
            # start_pos = '-' + str(100000 * floor(min(read.reference_start, read.next_reference_start)) / 100000)  
            read_info += [id2seq[read.reference_name], min(read.reference_start, read.next_reference_start)] # Mapped-to reference
        elif read.reference_name: 
            read_info += [id2seq[read.reference_name]]
        elif est_fragments:
            read_info += ['None', '-1']
        else:
            read_info += ['None']
        # if read.has_tag('BX'): # TODO correct this for EMA
        #     read_barcode = read.get_tag('BX')[2:-2]
        #     read_info_tbl.append([read_name, read_barcode, direction, read_aln_info])
        # else:
        read_info_tbl.append(read_info)
    col_names = ['Name', 'Direction', 'Reference']
    col_names.append('Start') if est_fragments else None
    pd.DataFrame(read_info_tbl).to_csv(join(outdir, prefix(bam) + '.csv'), index = False, header = col_names)    


def add_barcodes(fastq, annot_csv, outdir):
    """Re-adds FastQ barcodes to bwa-based annotation file generated by bam_to_annotate. Output is [read, direction, reference name, (insert start position), barcode]. bwa-specific."""
    read2bc = {}
    logger(f'Loading read cloud information from {fastq}')
    with open(fastq,'r') as fq:
        for i, line in enumerate(fq):
            if (i % 4) == 0:
                name_components = line.strip()[1:].split()
                if (len(name_components) > 1): # If the read has a barcode...
                    read2bc[name_components[0]] = name_components[1].split(':')[2].split('-')[0]
                else:
                    read2bc[name_components[0]] = 'None'
            if (i % 10000000) == 0:
                logger(f'Finished processing {i} lines')
    df = pd.read_csv(annot_csv)
    logger(f'Adding barcodes to the read mapping information')
    df['Barcode']= df['Name'].map(read2bc)
    logger(f'Sorting in order of barcodes, read names, and direction')
    df.sort_values(['Barcode', 'Reference'], axis = 0, ascending = [True, True], inplace = True, na_position ='last') 
    df.to_csv(join(outdir, prefix(annot_csv) + '.final.csv'), index = False, header = True)


def split_into_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]   # yields successive n-sized chunks of data


def subdiv_annotations(infile, num_chunks, outdir):
    """Subsets total set of read clouds into smaller sets of read clouds so that the next step (concat_annotations) doesn't take forever. bwa-specific."""
    logger(f'Started loading the annotations')
    df = pd.read_csv(infile, header = 0)
    logger(f'Finished loading the annotations')
    barcodes = df.Barcode.unique().tolist()
    # barcodes.remove('None')
    logger(f'{len(barcodes)} read clouds')
    chunk_size = (int)(len(barcodes) / int(num_chunks)) + 1
    barcode_sublists = list(split_into_chunks(barcodes, chunk_size))
    logger(f'{len(barcode_sublists)} sublists of approximately {chunk_size} read clouds each')
    for i, bl in enumerate(barcode_sublists):
        df_b = df[df.Barcode.isin(bl)]
        df_b.to_csv(join(outdir, '.'.join([prefix(infile), str(i), 'csv'])), index = False)
        if (i % 10) == 0:
            logger(f'Finished processing {i} barcode sublists')


def species_annot(df, barcodes): 
    """[read, direction, reference name, barcode]"""
    fw = []
    rv = []
    for i, b in enumerate(barcodes):
        df_b = df.loc[df['Barcode'] == b]
        references = df_b.Reference.unique().tolist()
        if len(references) > 1:
            for j in range(len(df_b)):
                read_ref = references.index(df_b.iloc[j,2]) # Index/enhanced cloud of the read's reference
                if df_b.iloc[j,1]: 
                    rv.append([df_b.iloc[j,0], read_ref])
                else:
                    fw.append([df_b.iloc[j,0], read_ref])
        if (i % 10000) == 0:
            logger(f'{i} original read clouds processed')
    return fw, rv


def fragment_annot(df, barcodes):
    """[read, direction, reference name, insert start position, barcode]"""
    fw = []
    rv = []
    for i, b in enumerate(barcodes):
        df_b = df.loc[df['Barcode'] == b]
        cloud_idx = 0
        for j in df_b.Reference.unique().tolist():
            ref_idx = cloud_idx
            df_r = df_b.loc[df_b['Reference'] == j]
            start_pos = df_r.Start.tolist()[0] # Leftmost start on reference j. Should be pre-sorted from SAM
            for k in range(len(df_r)):
                if df_r.iloc[k,3] > start_pos + 200000: # Read start outside 2 x max. est. fragment size
                    ref_idx += 1
                    start_pos = df_r.iloc[k,3]
                read_info = [df_r.iloc[k,0], ref_idx, j + '-' + str(ref_idx)] # ref-enh_num
                rv.append(read_info) if df_r.iloc[k,1] else fw.append(read_info)
            cloud_idx = ref_idx + 1
    return fw, rv


def concat_annotations(infile, outdir, est_fragments):
    """Aggregates reference mapping information in each read cloud and generates enhanced groupings. Output is [read, enh_num, (ref-enh_num)]."""
    df = pd.read_csv(infile)
    logger(f'Finished loading the annotations')
    barcodes = df.Barcode.unique()
    logger(f'{len(barcodes)} read clouds')
    fw, rv = fragment_annot(df, barcodes) if est_fragments else species_annot(df, barcodes)
    pd.DataFrame(fw).to_csv(join(outdir, 'R1', basename(infile)), index = False, header = False)    
    pd.DataFrame(rv).to_csv(join(outdir, 'R2', basename(infile)), index = False, header = False)    


def load_from_dir(enh_dir):
    """TODO: Probably deprecate in the future since concatenation function is in the pipeline now."""
    name_enh_dict = {}
    enh_list = [f for f in listdir(enh_dir) if isfile(join(enh_dir, f))]
    for i, f in enumerate(enh_list):
        name_enh_dict.update(load_from_csv(join(enh_dir,f)))
        if (i % 20 == 0):
            logger(f'Loaded {i} enhanced CSVs')            
    return name_enh_dict


def load_from_csv(enh_csv):
    name_enh_dict = {}
    with open(enh_csv,'r') as ef:
        for line in ef:
            read_info = line.strip().split(',')
            name_enh_dict[read_info[0]] = read_info[1]
    return name_enh_dict


def fastq_enhance(in_fq, enh_csv, outdir, prefix):
    """Replaces original cloud number with tool-enhanced number."""
    enhanced_dict = load_from_csv(enh_csv)

    # Replace the original groupings with enhanced groupings
    in_fq_parts = basename(in_fq).split('.')
    out_fq = join(outdir, '.'.join([prefix, in_fq_parts[1], in_fq_parts[2]]))
    with open(out_fq, 'w') as of:
        with open(in_fq, 'r') as fq: 
            for i, line in enumerate(fq):
                if (i % 4) == 0: # Read name
                    out_read_name = line.strip()
                    name_components = out_read_name.split() 
                    if len(name_components) > 1: # If read is i) in a cloud, replace it with enhanced cloud numbers
                        if name_components[0][1:] in enhanced_dict: # Get rid of this post-testing
                            enh_num = enhanced_dict[name_components[0][1:]] # Enhanced cloud numbers
                            out_read_name = out_read_name[:-1] + str(enh_num) # Tag it on
                    of.write(out_read_name + '\n')
                else:
                    of.write(line)
                if (i % 10000000 == 0):
                    logger(f'Finished processing {i} lines from {in_fq}')            

