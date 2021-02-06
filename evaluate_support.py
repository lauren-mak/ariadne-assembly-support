
from csv import reader, writer
from collections import Counter
from datetime import datetime
from os.path import basename, join, exists
import subprocess
from sys import stderr

# Pandas and numpy for data manipulation
import numpy as np
import pandas as pd
from math import log2, fabs
from sklearn.metrics.pairwise import manhattan_distances

# Matplotlib and colormaps for plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm


def logger(msg):
    print(f'{datetime.now()}: {msg}', file = stderr)


def fastq_to_table(fastq, mapping_tbl, outdir):
    """Generates tables pairing enhanced cloud numbers with the actual reference sequence. Output is [barcode, enhanced_num, reference name]"""
    logger(f'Loading read cloud information from {fastq}')
    full_name_list = []
    with open(fastq,'r') as fq:
        for i, line in enumerate(fq):
            if (i % 4) == 0 and 'BX:Z:' in line:
                full_name_list.append(line.strip()[1:].replace(' BX:Z:','-').split('-')) # [[name,letters,cloud_num]]

    logger(f'Loading reference sequence information from {mapping_tbl}')
    read_seq_map = {}
    with open(mapping_tbl, 'r') as mf:
        for line in mf:
            read_to_seq = line.strip().split(',')
            read_seq_map[read_to_seq[0]] = read_to_seq[2]

    logger(f'Adding reference sequence information to reads')
    for i, read_name in enumerate(full_name_list):
        read_name.append(read_seq_map[read_name[0]])
        if (i % 1000000) == 0:
            logger(f'{i} reads have been processed')
    full_name_df = pd.DataFrame(full_name_list, columns = ['Read_Name', 'Barcode', 'Cloud_Num', 'Ref_Seq'])
    del full_name_df['Read_Name']
    prefix = basename(fastq).split('.')
    full_name_df.to_csv(outdir + '/' + '.'.join(prefix[0:2]) + '.csv', header=False, index=False)


def load_tbl(tbl):
    out = []
    with open(tbl, 'r') as tf:
        for line in tf:
            out.append(line.strip().split(','))
    return out


def count_species(cloud_refs, sl):
    species_ct_lst = []
    for s in sl:
        species_ct_lst.append(cloud_refs.count(s))
    species_ct_lst.append(len(cloud_refs) - sum(species_ct_lst))
    return species_ct_lst


def evaluate_cloud(species_cts):
    size = sum(species_cts)
    entropy = 0
    for ct in species_cts:
        p = ct / size # p < 1 unless count = size i.e.: a pure cluster.
        if p > 0:
            entropy += p * log2(p)
    return size, max(species_cts) / size, fabs(entropy) # entropy < 0 unless count = size i.e.: a pure cluster. 


def merge_clouds_species(forward_lst, reverse_lst, sl):
    summary_info = [] # [letters, number, size, purity, entropy]
    species_info = [] # [species_counts...]
    cloud_counts = []
    num_enhanced_clouds = 0 
    curr_barcode = forward_lst[0][0]
    cloud_ref_dict = {}
    for i in range(len(forward_lst)):
        barcode = forward_lst[i][0]
        if barcode != curr_barcode:
            for c in cloud_ref_dict:
                species_cts = count_species(cloud_ref_dict[c], sl)
                species_info.append([curr_barcode, c] + species_cts)
                size, purity, entropy = evaluate_cloud(species_cts)
                summary_info.append([curr_barcode, c, size, purity, entropy])
                cloud_counts.append(size)
            if (num_enhanced_clouds % 100000) == 0:
                logger(f'Processed {num_enhanced_clouds} enhanced clouds')
                num_enhanced_clouds += 1
            curr_barcode = barcode
            cloud_ref_dict = {}
        cloud_num = forward_lst[i][1] 
        if cloud_num not in cloud_ref_dict:
            cloud_ref_dict[cloud_num] = []
        cloud_ref_dict[cloud_num].append(forward_lst[i][2])
        cloud_ref_dict[cloud_num].append(reverse_lst[i][2])
    for c in cloud_ref_dict:
        species_cts = count_species(cloud_ref_dict[c], sl)
        species_info.append([curr_barcode, c] + species_cts)
        size, purity, entropy = evaluate_cloud(species_cts)
        summary_info.append([curr_barcode, c, size, purity, entropy])
        cloud_counts.append(size)
    return species_info, summary_info, cloud_counts


def merge_clouds_fragments(forward_lst, reverse_lst):
    summary_info = [] # [letters, number, size, purity, entropy]
    cloud_counts = []
    num_enhanced_clouds = 0 
    curr_barcode = forward_lst[0][0]
    cloud_ref_dict = {}
    for i in range(len(forward_lst)):
        barcode = forward_lst[i][0]
        if barcode != curr_barcode:
            for c in cloud_ref_dict:
                frag_cts = Counter(cloud_ref_dict[c]).values()
                size, purity, entropy = evaluate_cloud(frag_cts)
                summary_info.append([curr_barcode, c, size, purity, entropy])
                cloud_counts.append(size)
            if (num_enhanced_clouds % 100000) == 0:
                logger(f'Processed {num_enhanced_clouds} enhanced clouds')
                num_enhanced_clouds += 1
            curr_barcode = barcode
            cloud_ref_dict = {}
        cloud_num = forward_lst[i][1] 
        if cloud_num not in cloud_ref_dict:
            cloud_ref_dict[cloud_num] = []
        cloud_ref_dict[cloud_num].append(forward_lst[i][2])
        cloud_ref_dict[cloud_num].append(reverse_lst[i][2])
    for c in cloud_ref_dict:
        frag_cts = Counter(cloud_ref_dict[c]).values()
        size, purity, entropy = evaluate_cloud(frag_cts)
        summary_info.append([curr_barcode, c, size, purity, entropy])
        cloud_counts.append(size)
    return summary_info, cloud_counts


def dataset_summary(cloud_counts, prefix):
    num_unique_barcodes = len(cloud_counts)
    size_list = sorted(cloud_counts)
    size_cts = Counter(size_list)
    size_df = pd.DataFrame.from_dict(size_cts, orient='index')
    logger('Finished generating dataframe, plotting frequency information')
    plt.plot(size_df)
    plt.savefig(prefix + '.cumulative_counts.png', format = 'png', dpi = 1200)
    # size_df['Count_Freq'] = size_df.iloc[:,0]/num_unique_barcodes
    # size_df['Cum_Freq'] = size_df.iloc[:,0].cumsum()/num_unique_barcodes

    middle_index = int(len(size_list)/2)
    size_median = size_list[middle_index]
    if len(size_list) % 2 == 0:
        size_median = (size_list[middle_index - 1] + size_list[middle_index])/2
    summary_info = [['Num_Clouds', num_unique_barcodes], 
                    ['Min_Num_Reads', size_list[0]], 
                    ['Max_Num_Reads', size_list[-1]],
                    ['Mean_Num_Reads', sum(size_list) / num_unique_barcodes],
                    ['StDev_Num_Reads', np.std(size_list)],
                    ['Med_Num_Reads', size_median]]
    with open(prefix + '.summary.csv', 'w') as of: # Barcode,Cloud_Num,Species_0,Species_1,Species_2...
        ow = writer(of)
        ow.writerows(summary_info)
    # summary_info = [pd.Series([num_unique_barcodes, 'Num_Clouds', 'NA'], index=size_df.columns),
    #             pd.Series([size_list[0], 'Min_Num_Reads', 'NA'], index=size_df.columns),
    #             pd.Series([size_list[-1], 'Max_Num_Reads', 'NA'], index=size_df.columns),
    #             pd.Series([sum(size_list) / num_unique_barcodes, 'Mean_Num_Reads', 'NA'], index=size_df.columns),
    #             pd.Series([size_median, 'Med_Num_Reads', 'NA'], index=size_df.columns)]
    # size_df = size_df.append(summary_info, ignore_index=True)
    # size_df.to_csv(outdir + '/' + prefix + '.dataset.csv', header = False)


def generate_summaries(forward_tbl, reverse_tbl, outdir, id_csv): 
    """Calculate summary statistics for each and across all enhanced read clouds using the matched enhanced-actual information above. Output are [barcode, enhanced_num, size, purity, entropy] and [barcode, enhanced_num, species list]"""
    logger('Loading forward FastQ information')
    forward_lst = load_tbl(forward_tbl)
    logger('Loading reverse FastQ information')
    reverse_lst = load_tbl(reverse_tbl)

    prefix = join(outdir, basename(forward_tbl).split('.')[0])
    if id_csv:
        species_df = pd.read_csv(id_csv, header = None)
        species_lst = list(species_df.iloc[:,1])
        species_info, summary_info, cloud_counts = merge_clouds_species(forward_lst, reverse_lst, species_lst)
        with open(prefix + '.species.csv', 'w') as of: # Barcode,Cloud_Num,Species_0,Species_1,Species_2...
            ow = writer(of)
            ow.writerow(['Barcode', 'Cloud_Num'] + species_lst + ['None'])
            ow.writerows(species_info)
    else:
        summary_info, cloud_counts = merge_clouds_fragments(forward_lst, reverse_lst)

    dataset_summary(cloud_counts, prefix)
    with open(prefix + '.statistics.csv', 'w') as of:
        ow = writer(of)
        ow.writerow(['Barcode', 'Cloud_Num', 'Size', 'Purity', 'Entropy'])
        ow.writerows(summary_info)


def load_dataframes(file_prefixes, outdir):
    all_df_list = []
    size_filter = 2 # Exclude all read clouds that are composed of a single pair of reads 
    for fp in file_prefixes:
        df = pd.read_csv('{}/{}.statistics.csv'.format(outdir,fp), header = 0)
        filtered_df = df[df['Size'] > size_filter]
        all_df_list.append(filtered_df)
    return all_df_list


def make_main_graphs(all_df_list, distances, param_name, outdir):
    # Matplotlib settings needed to create complex bar plots
    mpl.rcParams['font.size'] = 6
    mpl.rcParams['figure.dpi'] = 250
    scale_y = 1e3

    bins = np.linspace(0, 1, 20)
    purity_list = []
    hist_list = []
    for df in all_df_list:
        curr_lst = df[param_name]
        purity_list.append(curr_lst)
        curr_hist, bins = np.histogram(curr_lst, bins=bins)
        curr_freqs = np.divide(curr_hist, len(curr_lst))
        hist_list.append(curr_freqs)

    fig = plt.figure()
    # norm = colors.Normalize(vmin = 0, vmax = len(all_df_list) - 1)
    # normed_distances = list(map(lambda x: norm(x), range(len(all_df_list))))
    # cmap = cm.RdPu_r(normed_distances, bytes=False) # Luminance increases monotonically
    cmap = cm.Set2(np.linspace(0, 1, len(all_df_list)))

    plt.hist(purity_list, bins, label = distances, color = cmap)
    plt.legend(prop={'size': 6}, title = 'Search\nDistance')
    plt.xlabel(param_name)
    plt.ylabel('Num. Read Clouds')
    fig.savefig('{}/{}.png'.format(outdir, param_name), format = 'png', dpi = 1200)
    plt.clf()

    num_rows = int(len(distances)/3 + 1)
    fig, axs = plt.subplots(num_rows, 3)
    x = []
    for i in range(num_rows):
        x += [i]*3
    y = list(range(3)) * num_rows
    for i in range(len(hist_list)):
        axs[x[i],y[i]].plot(bins[:-1], hist_list[i], color = cmap[i])
        axs[x[i],y[i]].set_title(f'Search Distance = {distances[i]}')
    for ax in axs.flat:
        ax.set(xlabel='Purity', ylabel='Prop. Read Clouds')
    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    fig.savefig('{}/{}_scaled.png'.format(outdir, param_name), format = 'png', dpi = 1200)


def make_secondary_graphs(param_list, distances, param_name, outdir): 
    # Matplotlib settings needed to create complex bar plots
    mpl.rcParams['font.size'] = 6
    mpl.rcParams['figure.dpi'] = 250
    scale_y = 1e3

    fig = plt.figure()
    # norm = colors.Normalize(vmin = 0, vmax = len(param_list) - 1)
    # normed_distances = list(map(lambda x: norm(x), range(len(param_list))))
    # cmap = cm.RdPu_r(normed_distances, bytes=False) # Luminance increases monotonically
    cmap = cm.Set2(np.linspace(0, 1, len(param_list)))

    param_std = list(map(lambda lst: np.std(lst), param_list))
    # y_error = [np.subtract(param_list, param_std), np.add(param_list, param_std)]
    
    plt.bar(distances, param_list, color = cmap) # yerr = y_error, 
    plt.xlabel('Search Distance')
    plt.ylabel(param_name)
    param2strng = param_name.replace(' ', '_').replace('.', '')
    fig.savefig('{}/{}.png'.format(outdir, param2strng), format = 'png', dpi = 1200)
    df = pd.DataFrame(list(zip(distances,param_list)), columns=['Search_Distances', param2strng])
    df.to_csv('{}/{}.tbl'.format(outdir, param2strng))


def evaluate_clouds(distances, prefixes, outdir):
    """Make summary graphs of the summary information table from generate_summaries()."""
    start_time = datetime.now()
    logger(f'Generating summary statistics for search distances {distances} of the dataset {"_".join(outdir.split("_")[:-1])}')

    search_distances = distances.split(',')
    file_prefixes = prefixes.split(',')
    all_df_list = load_dataframes(file_prefixes, outdir) 
    # Load alignment info tables for each search distance and filter for clouds of size < 2
    logger('Loaded the search distance dataframes ' + distances)

    # Comparing the purity and Shannon entropy distributions for each search distances
    make_main_graphs(all_df_list, search_distances, 'Purity', outdir)
    make_main_graphs(all_df_list, search_distances, 'Entropy', outdir)
    logger('Finished the main purity and entropy comparison graphs')

    # How many read clouds are there in total? 
    # num_read_clouds = list(map(lambda df: len(df.index), all_df_list))
    # make_secondary_graphs(num_read_clouds, search_distances, 'Num. Read Clouds', outdir)
    # What is the average size of a read-cloud?
    avg_cloud_size = list(map(lambda df: np.mean(df['Size']), all_df_list))
    make_secondary_graphs(avg_cloud_size, search_distances, 'Avg. Cloud Size', outdir)
    # # How many reads were deconvolved?
    # num_reads_decon = list(map(lambda df: df['Size'].sum(), all_df_list[1:]))
    # make_secondary_graphs(num_reads_decon, search_distances[1:], 'Num. Read Deconv.', outdir)
    # What is the average purity of the read clouds?
    avg_cloud_purity = list(map(lambda df: np.mean(df['Purity']), all_df_list))
    make_secondary_graphs(avg_cloud_purity, search_distances, 'Avg. Cloud Purity.', outdir)
    # What is the average purity of the read clouds?
    avg_cloud_entropy = list(map(lambda df: np.mean(df['Entropy']), all_df_list))
    make_secondary_graphs(avg_cloud_entropy, search_distances, 'Avg. Cloud Entropy.', outdir)
    logger('Finished the accessory graphs')


def map_read_clouds(fg, prefix, fa, seq_to_cld_num):
    if not exists(prefix + '.tsv'):
        subprocess.run(['/Users/laurenmak/Programs/Bandage_Mac_v0_8_1/Bandage.app/Contents/MacOS/Bandage', 'querypaths', fg, fa, prefix])

    # Match Path information to read cloud (fragment) number
    bandage_tbl = pd.read_csv(prefix + '.tsv', header = 0, sep = '\t')
    clean_lst = [s.replace('(', '').replace(')', '').replace('+', '').replace('-', '').split() for s in bandage_tbl['Path'].tolist()] # Clean and split Path info into [[start,node_name,end]]
    name_list = bandage_tbl['Query'].tolist() # Extract only alignment into
    full_read_seqs = seq_to_cld_num.keys()
    for i, n in enumerate(name_list):
        clean_lst[i] = [seq_to_cld_num[n]] + clean_lst[i] # [[cld_num,start,node_name,end]]
    return pd.DataFrame(clean_lst, columns = ['Cloud_Num', 'Start', 'Node', 'End'])


def fragment_comparison(clean_df, node_to_node, depth, prefix):
    aln_nodes = clean_df.Node.unique().tolist()
    fragments = clean_df.Cloud_Num.unique().tolist()
    frag_df_lst = []
    aln_df = pd.DataFrame(0, index = fragments, columns = aln_nodes)
    ctg_dct = {}
    for i, f in enumerate(fragments):
        fragment_df = clean_df.loc[clean_df['Cloud_Num'] == f]
        frag_lst = []
        frag_aln_nodes = fragment_df.Node.unique().tolist()
        frag_ctg_nodes = []
        for j, a in enumerate(aln_nodes):
            if a in frag_aln_nodes:
                node_df = fragment_df.loc[fragment_df['Node'] == a]
                min_start = int(min(node_df['Start']))
                max_end = int(max(node_df['End']))
                contiguous_edges = depth_based_search(a, depth, node_to_node, [])
                frag_lst.append([a, len(node_df), min_start, max_end, max_end - min_start + 1, len(contiguous_edges), '/'.join(contiguous_edges)])
                frag_ctg_nodes += contiguous_edges
                aln_df.iloc[i,j] = len(node_df)
        ctg_dct[f] = dict(Counter(frag_ctg_nodes))
        frag_df_lst.append(pd.DataFrame(frag_lst, index = [f] * len(frag_lst), columns = ['Node', 'Num_Reads', 'Start', 'End', 'Distance', 'Num_Contig_Nodes', 'Contiguous_Nodes']))
        logger(f'{f}th read cloud: {len(fragment_df)} reads, {len(frag_aln_nodes)} unique aligned edges, {len(ctg_dct[f])} unique contiguous edges')
    pd.concat(frag_df_lst).to_csv('.'.join([prefix, 'edges', 'csv']))
    return aln_df, ctg_dct


def edge_comparison(clean_df, node_to_node, depth, prefix):
    aln_nodes = clean_df.Node.unique().tolist()
    fragments = clean_df.Cloud_Num.unique().tolist()
    edge_lst = []
    frag_name_lst = []
    aln_ctg_df_lst = []
    for i, f in enumerate(fragments):
        fragment_df = clean_df.loc[clean_df['Cloud_Num'] == f]
        frag_aln_nodes = fragment_df.Node.unique().tolist()
        for j, a in enumerate(aln_nodes):
            if a in frag_aln_nodes:
                node_df = fragment_df.loc[fragment_df['Node'] == a]
                min_start = int(min(node_df['Start']))
                max_end = int(max(node_df['End']))
                contiguous_edges = depth_based_search(a, depth, node_to_node, [])
                edge_name = '-'.join([f,a])
                edge_lst.append([len(node_df), min_start, max_end, max_end - min_start + 1, len(contiguous_edges), '/'.join(contiguous_edges)])
                frag_name_lst.append(edge_name)
                aln_ctg_dct = dict(Counter(contiguous_edges)) # Contiguous edges only
                if a in aln_ctg_dct:
                    aln_ctg_dct[a] += len(node_df)
                else:
                    aln_ctg_dct[a] = len(node_df)
                # aln_ctg_df_lst.append(pd.DataFrame.from_dict(aln_ctg_dct))
                aln_ctg_df_lst.append(pd.DataFrame([aln_ctg_dct.values()], index = [edge_name], columns = aln_ctg_dct.keys()))
        logger(f'{f}th read cloud: {len(fragment_df)} reads, {len(frag_aln_nodes)} unique aligned edges')
    pd.DataFrame(edge_lst, index = frag_name_lst, columns = ['Num_Reads', 'Start', 'End', 'Distance', 'Num_Contig_Nodes', 'Contiguous_Nodes']).to_csv('.'.join([prefix, 'edges', 'csv']))
    aln_df = pd.concat(aln_ctg_df_lst, axis = 0, sort = True)
    aln_df.fillna(value = 0, inplace = True)
    return aln_df.astype(int) 


def depth_based_search(a, d, node_to_node, existing_edges):
    connected_edges = node_to_node[a]
    if d is 0 or len(connected_edges) is 0: return [] 
    contiguous_edges = existing_edges + connected_edges
    # print(f'{a} {d} {contiguous_edges}')
    for b in connected_edges:
        if b not in existing_edges:
            contiguous_edges += depth_based_search(b, d - 1, node_to_node, contiguous_edges)
    cleaned_edge_lst = list(set(contiguous_edges))
    # print(f'{a} {d} {cleaned_edge_lst}')
    return cleaned_edge_lst


def pairwise_graph_align(fastg, fasta_prefix, outdir, depth, fragment_mode, aligned_only):
    """Identify read cloud assembly graph alignments and pairwise differences."""

    # Match the read sequence to the read cloud number
    seq_to_cld_num = {} # Sequence instead of read name because paired-end file
    with open(fasta_prefix + '.R1.fasta', 'r') as f:
        for i, l in enumerate(f): 
            if (i % 2 == 0):  # >MN00867:6:000H2J7M3:1:22110:23088:18151 BX:Z:GGATTTATGTTTGAATGG-4
                seq_to_cld_num[l.split()[0][1:]] = l.strip()[-1] # {seq:cld_num}
    logger(f'Finished loading read cloud numbers from {fasta_prefix + ".R1.fasta"}')

    # Map all reads to assembly graph
    prefix = join(outdir, basename(fastg).split('.')[0])
    r1_clean_df = map_read_clouds(fastg, prefix + '.R1', fasta_prefix + '.R1.fasta', seq_to_cld_num)
    r2_clean_df = map_read_clouds(fastg, prefix + '.R2', fasta_prefix + '.R2.fasta', seq_to_cld_num)
    clean_df = pd.concat([r1_clean_df, r2_clean_df])
    logger(f'There are {len(clean_df)} read-graph alignments in total')

    # Match overlapping nodes to each other
    node_to_node = {}
    with open(fastg, 'r') as fg:
        for i, l in enumerate(fg): 
            if ':' in l: # Two nodes that overlap
                edge_ids = l.split('_') # Node name at indices 1, 6
                if edge_ids[1] not in node_to_node:
                    node_to_node[edge_ids[1]] = []
                node_to_node[edge_ids[1]].append(edge_ids[6])
    logger(f'Finished extracting {len(node_to_node)} edge overlaps from {fastg}')

    # For each (fragment) read cloud, count the number of aligned and contiguous nodes. Output the summary information for each fragment.
    fragments = clean_df.Cloud_Num.unique().tolist()
    if fragment_mode:
        aln_df, ctg_dct = fragment_comparison(clean_df, node_to_node, depth, prefix)
        if aligned_only:
            # Pairwise comparison between i) aligned...
            logger(f'Pairwise comparison between aligned edges')
            aln_df.to_csv('.'.join([prefix, 'aln', 'csv']))
            aln_dist = pd.DataFrame(manhattan_distances(aln_df), index = aln_df.index.values, columns = aln_df.index.values)
            scaled_aln_dist = aln_dist
            for i in range(len(fragments)):
                for j in range(len(fragments)):
                    scaled_aln_dist.iloc[i,j] = aln_dist.iloc[i,j] / ( sum(aln_df.iloc[i,:]) + sum(aln_df.iloc[j,:]) )
                # scaled_aln_dist.iloc[i,:] = aln_dist.iloc[i,:] / sum(aln_df.iloc[i,:])
            scaled_aln_dist.to_csv('.'.join([prefix, 'pairwise_aln', 'csv']))
        ctg_df = pd.DataFrame.from_dict(ctg_dct, orient = 'index') 
        ctg_df.fillna(value = 0, inplace = True)
        ctg_nodes = ctg_df.columns.values.tolist()
        for i in set(aln_nodes) & set(ctg_nodes):
            aln_df[i] += ctg_df[i]
            del ctg_df[i]
        aln_df = pd.concat([aln_df, ctg_df], axis = 1, sort = True).astype(int) 
    else:
        aln_df = edge_comparison(clean_df, node_to_node, depth, prefix)

    # ...and ii) aligned + contiguous nodes.
    logger(f'Pairwise comparison between all (aligned + contiguous) edges')
    aln_df.to_csv('.'.join([prefix, 'all', 'csv']))
    aln_ctg_dist = pd.DataFrame(manhattan_distances(aln_df), index = aln_df.index.values, columns = aln_df.index.values)
    scaled_aln_ctg_dist = aln_ctg_dist
    for i in range(len(aln_ctg_dist)):
        for j in range(len(aln_ctg_dist)):
            scaled_aln_ctg_dist.iloc[i,j] = aln_ctg_dist.iloc[i,j] / ( sum(aln_df.iloc[i,:]) + sum(aln_df.iloc[j,:]) )
        # scaled_aln_ctg_dist.iloc[i,:] = aln_ctg_dist.iloc[i,:] / sum(aln_df.iloc[i,:])
    scaled_aln_ctg_dist.to_csv('.'.join([prefix, 'pairwise_all', 'csv']))

