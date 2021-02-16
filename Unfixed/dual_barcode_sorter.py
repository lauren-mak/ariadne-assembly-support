def fastq_ingest(in_fq):
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
                logger(f'Finished processing line {i} of {in_fq}')
    return bc2reads

def dual_barcode_sorter(forward_fq, reverse_fq, outdir):
    """Sorts FastQ reads based on their barcode and cloud number as a paired set. Tags barcoded reads with '-1' if not already there."""
    bc2reads = fastq_ingest(forward_fq)
    fw_parts = basename(forward_fq).split('.')
    fw_out = join(outdir, '.'.join([fw_parts[0] + '_bsort', fw_parts[1], fw_parts[2]]))
    with open(fw_out, 'w') as of:
        for b in bc2reads: 
            for e in bc2reads[b]: 
                for line in bc2reads[b][e]: 
                    of.write(line)
    bc2reads.clear()