
from csv import reader, writer
from collections import Counter
from datetime import datetime
from os.path import basename, join, exists
import subprocess
from sys import stderr

# Pandas and numpy for data manipulation
import altair as alt
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
            read_seq_map[read_to_seq[0]] = read_to_seq[1]

    logger(f'Adding reference sequence information to reads')
    for i, read_name in enumerate(full_name_list):
        read_name.append(read_seq_map[read_name[0]] if read_name[0] in read_seq_map else -1)
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
    bc_status = {}
    gstd_dct = {}
    for i in range(len(forward_lst)):
        barcode = forward_lst[i][0]
        if barcode != curr_barcode:
            for c in cloud_ref_dict:
                frag_cts = Counter(cloud_ref_dict[c])
                size, purity, entropy = evaluate_cloud(frag_cts.values())
                entire_barcode = curr_barcode + '-' + c
                if purity < 1: 
                    bc_status[entire_barcode] = 'Under'
                for f in frag_cts.keys():
                    if f not in gstd_dct:
                        gstd_dct[f] = [] 
                    gstd_dct[f].append(entire_barcode)
                summary_info.append([entire_barcode, size, purity, entropy])
                cloud_counts.append(size)
            for f in gstd_dct:
                if len(gstd_dct[f]) == 1:
                    if gstd_dct[f][0] not in bc_status:
                        bc_status[gstd_dct[f][0]] = 'Accurate'
                else:
                    for b in gstd_dct[f]:
                        if b not in bc_status:
                            bc_status[b] = 'Over'
            if (num_enhanced_clouds % 100000) == 0:
                logger(f'Processed {num_enhanced_clouds} enhanced clouds')
                num_enhanced_clouds += 1
            curr_barcode = barcode
            cloud_ref_dict = {}
            gstd_dct = {}
        cloud_num = forward_lst[i][1] 
        if cloud_num not in cloud_ref_dict:
            cloud_ref_dict[cloud_num] = []
        cloud_ref_dict[cloud_num].append(forward_lst[i][2])
        cloud_ref_dict[cloud_num].append(reverse_lst[i][2])
    for c in cloud_ref_dict:
        frag_cts = Counter(cloud_ref_dict[c])
        size, purity, entropy = evaluate_cloud(frag_cts.values())
        summary_info.append([curr_barcode + '-' + c, size, purity, entropy])
        cloud_counts.append(size)
        entire_barcode = curr_barcode + '-' + c
        if purity < 1: 
            bc_status[entire_barcode] = 'Under'
        for f in frag_cts.keys():
            if f not in gstd_dct:
                gstd_dct[f] = [] 
            gstd_dct[f].append(entire_barcode)
    for f in gstd_dct:
        if len(gstd_dct[f]) == 1:
            if gstd_dct[f][0] not in bc_status:
                bc_status[gstd_dct[f][0]] = 'Accurate'
        else:
            for b in gstd_dct[f]:
                if b not in bc_status:
                    bc_status[b] = 'Over'
    return summary_info, cloud_counts, bc_status


def dataset_summary(cloud_counts, prefix):
    num_unique_barcodes = len(cloud_counts)
    size_list = sorted(cloud_counts)
    size_cts = Counter(size_list)
    size_df = pd.DataFrame.from_dict(size_cts, orient='index')
    logger('Finished generating dataframe, plotting frequency information')
    plt.plot(size_df)
    plt.savefig(prefix + '.cumulative_counts.png', format = 'png', dpi = 1200)

    middle_index = int(len(size_list)/2)
    size_median = size_list[middle_index]
    if len(size_list) % 2 == 0:
        size_median = (size_list[middle_index - 1] + size_list[middle_index])/2
    ct_summary = [['Num_Clouds', num_unique_barcodes], 
                    ['Min_Num_Reads', size_list[0]], 
                    ['Max_Num_Reads', size_list[-1]],
                    ['Mean_Num_Reads', sum(size_list) / num_unique_barcodes],
                    ['StDev_Num_Reads', np.std(size_list)],
                    ['Med_Num_Reads', size_median]]
    with open(prefix + '.summary.csv', 'w') as of: # Barcode,Cloud_Num,Species_0,Species_1,Species_2...
        ow = writer(of)
        ow.writerows(ct_summary)


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
        summary_tbl = pd.DataFrame(summary_info, columns = ['Barcode-Cloud_Num', 'Size', 'Purity', 'Entropy'])
    else:
        summary_info, cloud_counts, bc_status = merge_clouds_fragments(forward_lst, reverse_lst)
        summary_tbl = pd.DataFrame(summary_info, columns = ['Barcode-Cloud_Num', 'Size', 'Purity', 'Entropy'])
        summary_tbl['Status'] = summary_tbl['Barcode-Cloud_Num'].map(bc_status)
    
    dataset_summary(cloud_counts, prefix)
    summary_tbl.to_csv(prefix + '.statistics.csv', index = False)


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
        curr_hist, bins = np.histogram(curr_lst, bins = bins)
        curr_freqs = np.divide(curr_hist, len(curr_lst))
        hist_list.append(curr_freqs)

    fig = plt.figure()
    cmap = cm.Set2(np.linspace(0, 1, len(all_df_list)))

    plt.hist(purity_list, bins, label = distances, color = cmap)
    plt.legend(prop={'size': 6}, title = 'Search\nDistance')
    plt.xlabel(param_name)
    plt.ylabel('Num. Read Clouds')
    fig.savefig('{}/{}.png'.format(outdir, param_name), format = 'png', dpi = 1200)
    plt.clf()

    graphs_per_row = 4
    num_rows = int(len(distances)/graphs_per_row + 1)
    fig = plt.figure() # plt.subplots(num_rows, graphs_per_row, constrained_layout = True)
    x = []
    for i in range(num_rows):
        x += [i] * graphs_per_row
    y = list(range(graphs_per_row)) * num_rows
    for i in range(len(hist_list)):
        ax = plt.subplot2grid((num_rows, graphs_per_row), (x[i], y[i])) #, constrained_layout = True)
        ax.plot(bins[:-1], hist_list[i], color = cmap[i])
        ax.set_title(f'Search Distance = {distances[i]}')
        ax.set(xlabel=param_name, ylabel='Prop. Read Clouds')
        # axs[x[i],y[i]].plot(bins[:-1], hist_list[i], color = cmap[i])
        # axs[x[i],y[i]].set_title(f'Search Distance = {distances[i]}')
    # for ax in axs.flat:
    #     ax.set(xlabel='Purity', ylabel='Prop. Read Clouds')
    # Hide x labels and tick labels for everything but the left- and right-most plots. Removed because different axis scales.
    # for ax in axs.flat:
    #     ax.label_outer()
    fig.tight_layout() 
    fig.savefig('{}/{}_scaled.png'.format(outdir, param_name), format = 'png', dpi = 1200)


def make_secondary_graphs(param_list, distances, param_name, outdir): # param_stdev
    # Matplotlib settings needed to create complex bar plots
    mpl.rcParams['font.size'] = 6
    mpl.rcParams['figure.dpi'] = 250
    scale_y = 1e3

    fig = plt.figure()
    cmap = cm.Set2(np.linspace(0, 1, len(param_list)))

    # y_error = [np.subtract(param_list, param_stdev), np.add(param_list, param_stdev)]
    plt.bar(distances, param_list, color = cmap) # yerr = y_error, 
    plt.xlabel('Search Distance')
    plt.ylabel(param_name)
    param2strng = param_name.replace(' ', '_').replace('.', '')
    fig.savefig('{}/{}.png'.format(outdir, param2strng), format = 'png', dpi = 1200)


def make_size_graphs(all_df_list, distances, param_name, outdir):
    # Matplotlib settings needed to create complex bar plots
    mpl.rcParams['font.size'] = 6
    mpl.rcParams['figure.dpi'] = 250
    scale_y = 1e3

    bins = np.linspace(0, 150, 20)
    purity_list = []
    for df in all_df_list:
        purity_list.append(df[param_name])

    graphs_per_row = 4
    num_rows = int(len(distances)/graphs_per_row + 1)
    fig = plt.figure() # plt.subplots(num_rows, graphs_per_row, constrained_layout = True)
    cmap = cm.Set2(np.linspace(0, 1, len(purity_list)))
    x = []
    for i in range(num_rows):
        x += [i] * graphs_per_row
    y = list(range(graphs_per_row)) * num_rows
    for i, p in enumerate(purity_list):
        ax = plt.subplot2grid((num_rows, graphs_per_row), (x[i], y[i])) #, constrained_layout = True)
        ax.hist(p, bins, color = cmap[i])
        # ax.plot(bins[:-1], hist_list[i], color = cmap[i])
        ax.set_title(f'Search Distance = {distances[i]}')
        ax.set(xlabel=param_name, ylabel='Prop. Read Clouds')
    fig.tight_layout() 
    fig.savefig('{}/{}.png'.format(outdir, param_name), format = 'png', dpi = 1200)


def make_status_graphs(param_df, distances, param_name, outdir): # param_stdev
    # Matplotlib settings needed to create complex bar plots
    mpl.rcParams['font.size'] = 6
    mpl.rcParams['figure.dpi'] = 250
    scale_y = 1e3

    # fig = plt.figure()
    cmap = cm.Set2(np.linspace(0, 1, param_df.shape[1]))

    ax = param_df.plot(kind = 'bar', color = cmap)
    fig = ax.get_figure()
    # y_error = [np.subtract(param_list, param_stdev), np.add(param_list, param_stdev)]
    # plt.bar(distances, param_df, color = cmap) # yerr = y_error, 
    ax.set_xlabel('Search Distance')
    ax.set_ylabel(param_name)
    param2strng = param_name.replace(' ', '_').replace('.', '')
    fig.savefig('{}/{}.png'.format(outdir, param2strng), format = 'png', dpi = 1200)


def evaluate_clouds(distances, prefixes, outdir):
    """Make summary graphs of the summary information table from generate_summaries()."""
    start_time = datetime.now()
    logger(f'Generating summary statistics for search distances {distances} of the dataset {"_".join(outdir.split("_")[:-1])}')

    # Load alignment info tables for each search distance and filter for clouds of size < 2
    search_distances = distances.split(',')
    file_prefixes = prefixes.split(',')
    all_df_list = load_dataframes(file_prefixes, outdir) 
    logger('Loaded the search distance dataframes ' + distances)

    # Comparing the purity Shannon entropy, and size distributions for each search distances
    avg_param_lst = []
    std_param_lst = []
    for i, param in enumerate(['Purity', 'Entropy']):
        make_main_graphs(all_df_list, search_distances, param, outdir)
        avg_param_lst.append(list(map(lambda df: np.mean(df[param]), all_df_list)))
        std_param_lst.append(list(map(lambda df: np.std(df[param]), all_df_list)))
        make_secondary_graphs(avg_param_lst[i], search_distances, 'Avg. Cloud ' + param, outdir) # std_param_lst[i]
        logger('Finished ' + param + ' comparison graphs')

    # Comparing the size of the deconvolved read clouds.
    make_size_graphs(all_df_list, search_distances, 'Size', outdir)
    avg_param_lst.append(list(map(lambda df: np.mean(df['Size']), all_df_list)))
    std_param_lst.append(list(map(lambda df: np.std(df['Size']), all_df_list)))
    make_secondary_graphs(avg_param_lst[-1], search_distances, 'Avg. Cloud Size', outdir) # std_param_lst[i]
    logger('Finished Size comparison graphs')

    # Comparing the deconvolution status (i.e.: is it over, under, or accurately deconvolved) for each search distance
    status_df_list = []
    for i, df in enumerate(all_df_list):
        status_df_list.append(pd.DataFrame(dict(Counter(df['Status'])), index = [search_distances[i]]))
    status_df = pd.concat(status_df_list).fillna(0).T
    normed_df = status_df.div(status_df.sum(axis = 0), axis = 1)
    make_status_graphs(status_df, search_distances, 'Cloud Deconv. Status', outdir)
    make_status_graphs(normed_df, search_distances, 'Normalized Deconv. Status', outdir)

    # TODO How many reads were deconvolved?

    avg_stats_df = pd.DataFrame(list(zip(search_distances, avg_param_lst[0], std_param_lst[0], avg_param_lst[1], std_param_lst[1], \
        avg_param_lst[2], std_param_lst[2])), 
        columns = ['Search_Distances', 'Avg. Cloud Purity', 'Std. Dev.', 'Avg. Cloud Entropy', 'Std. Dev.', 'Avg. Cloud Size', \
        'Std. Dev.'])
    avg_stats_df.set_index('Search_Distances', inplace = True)
    df = pd.concat([avg_stats_df, status_df.T, normed_df.T], axis = 1)
    df.round(decimals = 2).to_csv('{}/{}.tbl'.format(outdir, 'Avg_Summary_Stats'))


def graph_assembly_stats(base_dir):
    """Make gridded comparison of metaQUAST assembly statistics."""

    na50_tbl = pd.read_csv('NA50.csv', header = 0)
    # Load metaQUAST tables wth the four relevant assembly statistics
    mq_tbls = []
    for dp in ['mock5_10x', 'mock6_lsq', 'mock20_10x_100m', 'mock20_tsq_100m']: 
        tbls = []
        dfs = []
        na50_sub = na50_tbl.loc[na50_tbl['Dataset'] == dp]
        renamed_df = na50_sub.rename(columns = {'Fragments': 'Reference', '5000': 'Ariadne'})
        if dp is 'mock20_tsq_100m':
            tmp_df = renamed_df.rename(columns = {'10000': '1000', '15000': '2000', '20000': '4000'}).astype('str')
            clean_df = tmp_df[['Illumina', 'No_Deconv', 'Reference', 'Ariadne', '1000', '2000', '4000']].astype('str')
        else:
            clean_df = renamed_df[['Illumina', 'No_Deconv', 'Reference', 'Ariadne', '10000', '15000', '20000']].astype('str')
        clean_df.replace(['-'], ['0'], inplace = True)
        aln_only = clean_df.drop(clean_df.loc[clean_df.index=='not_aligned'].index)
        tbls.append(aln_only.astype('int64')) 
        dfs.append(pd.DataFrame(aln_only.astype('int64')))
        full_prefix = join(base_dir, dp + '_analyses', 'metaQUAST', 'summary', 'TSV')
        for f in ['Largest_alignment.tsv', 'Misassembled_contigs_length.tsv', 'Total_length.tsv']:
            tmp_df = pd.read_csv(join(full_prefix, f), sep = '\t', header = 0, index_col = 0)
            renamed_df = tmp_df.rename(columns = {'Fragments': 'Reference', '5000': 'Ariadne'})
            if dp is 'mock20_tsq_100m':
                clean_df = renamed_df[['Illumina', 'No_Deconv', 'Reference', 'Ariadne', '1000', '2000', '4000']].astype('str')
            else:
                clean_df = renamed_df[['Illumina', 'No_Deconv', 'Reference', 'Ariadne', '10000', '15000', '20000']].astype('str')
            clean_df.replace(['-'], ['0'], inplace = True)
            aln_only = clean_df.drop(clean_df.loc[clean_df.index=='not_aligned'].index)
            tbls.append(aln_only.astype('int64')) 
            dfs.append(pd.DataFrame(aln_only.astype('int64')))
        pd.concat(tbls).to_csv(dp + '.csv')
        mq_tbls.append(tbls)

    # Scale Reference and Ariadne approaches by the No_Deconv column
    mq_dct = {'Dataset': [], 'Deconv_Method': [], 'NA50': [], 'Largest_Aln': [], 'Misassembly_Rate': []}
    for i, dp in enumerate(['MOCK5 10x', 'MOCK5 LoopSeq', 'MOCK20 10x', 'MOCK20 TELLSeq']): 
        num_species = len(mq_tbls[i][0])
        mq_dct['Dataset'].extend([dp] * num_species * 3)
        mq_dct['Deconv_Method'].extend(['No_Deconv'] * num_species + ['Reference'] * num_species + ['Ariadne'] * num_species)
        for j, f in enumerate(['NA50', 'Largest_Aln', 'Misassembly_Rate']):
            clean_df = mq_tbls[i][j]
            if f is 'Misassembly_Rate': # First divide the whole Misassembled_contigs_length by Total_length
                clean_df = mq_tbls[i][j].div(mq_tbls[i][j + 1])
            for d in ['No_Deconv', 'Reference', 'Ariadne']:
                mq_dct[f].extend(clean_df[d])
            #     for d in ['Reference', 'Ariadne']:
            #         mq_dct[f].extend(clean_df[d] / clean_df['No_Deconv'])
            # else:
            #     for d in ['Reference', 'Ariadne']:
            #         mq_dct[f].extend(clean_df[d] - clean_df['No_Deconv'])
    pd.DataFrame.from_dict(mq_dct).to_csv(join(base_dir, 'metaQUAST_all_three.csv'))


def graph_cloud_stats(base_dir, num_clouds, param, scale):
    """Make gridded comparison of read cloud summary statistics."""

    # Load read cloud statistics from generate_summaries()
    deconv_tbls = []
    dataset_prefixes = ['mock5_10x', 'mock6_lsq', 'mock20_10x_100m', 'mock20_tsq_100m']
    deconv_prefixes = ['no_deconv', '_frg_bsort', '5000']
    for i, dp in enumerate(['MOCK5 10x', 'MOCK5 LoopSeq', 'MOCK20 10x', 'MOCK20 TELLSeq']): 
        for j, d in enumerate(['No_Deconv', 'Reference', 'Ariadne']):
            full_prefix = join(base_dir, dataset_prefixes[i] + '_analyses') + '/'
            full_prefix += dataset_prefixes[i] if j == 1 else ''
            full_prefix += deconv_prefixes[j]
            tmp_df = pd.read_csv(full_prefix + '.statistics.csv', header = 0)[['Size', param]]
            tmp_df['Dataset'] = [dp] * len(tmp_df)
            tmp_df['Deconv_Method'] = [d] * len(tmp_df)
            downsampled_df = tmp_df.sample(n = int(num_clouds))
            if scale:
                downsampled_df['Size'] = downsampled_df['Size'] / max(downsampled_df['Size'])
            deconv_tbls.append(downsampled_df)
    deconv_df = pd.concat(deconv_tbls)
    deconv_df.to_csv(join(base_dir, param + '_Cloud_Stats.csv'))


def tidy_quast_report(prefix, outdir):
    """Trim full metaQUAST report to more easily interpretable and camera-ready versions."""
    report_df = pd.read_csv(join(outdir, 'report.tsv'), sep = '\t', header = 0, index_col = 0)
    kept_info = [ 'Genome fraction (%)', 'Duplication ratio', 'Largest alignment', 'Total aligned length',\
        'NA50', '# misassemblies', '# misassembled contigs', 'Misassembled contigs length', '# unaligned contigs',\
        'Unaligned length', '# N\'s per 100 kbp', '# mismatches per 100 kbp', '# contigs', 'Total length (>= 0 bp)' ]
    trimmed_report = report_df.loc[report_df.index.intersection(kept_info)]
    sorted_report = trimmed_report.reindex(kept_info, axis = 'index')
    sorted_report.to_csv(join(outdir, prefix + '.metaQUAST.csv'))
    sorted_report.to_latex(join(outdir, prefix + '.metaQUAST.tex'))


def map_read_clouds(fg, prefix, fa, seq_to_cld_num):
    if not exists(prefix + '.tsv'):
        subprocess.run(['/Users/laurenmak/Programs/Bandage_Mac_v0_8_1/Bandage.app/Contents/MacOS/Bandage', 'querypaths', fg, fa, prefix])

    # Match Path information to read cloud (fragment) number
    bandage_tbl = pd.read_csv(prefix + '.tsv', header = 0, sep = '\t')
    clean_lst = [] # [s.replace('(', '').replace(')', '').replace('+', '').replace('-', '').split() for s in bandage_tbl['Path'].tolist()] # Clean and split Path info into [[start,node_name,end]]
    for r in bandage_tbl['Path'].tolist():
        tmp = r.strip().split()
        cleaned_info = []
        if '(' not in tmp[0]:
            cleaned_info.append(0)
        else:
            cleaned_info.append(tmp[0].replace('(', '').replace(')', ''))
        cleaned_info.append(tmp[1].replace('+', '').replace('-', '').replace(',', ''))
        if '(' not in tmp[-1]:
            cleaned_info.append(0)
        else:
            cleaned_info.append(tmp[-1].replace('(', '').replace(')', ''))
        clean_lst.append(cleaned_info)
    name_list = bandage_tbl['Query'].tolist() # Extract only alignment into
    full_read_seqs = seq_to_cld_num.keys()
    for i, n in enumerate(name_list):
        clean_lst[i] = [seq_to_cld_num[n]] + clean_lst[i] # [[cld_num,start,node_name,end]]
        if len(clean_lst[i]) > 4: 
            clean_lst[i] = clean_lst[i][0:3] + [clean_lst[i][-1]]
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
                edge_name = '-'.join([str(f),str(a)])
                edge_lst.append([len(node_df), min_start, max_end, max_end - min_start + 1, len(contiguous_edges), str(contiguous_edges)])
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
                seq_to_cld_num[l.split()[0][1:]] = l.strip().split('-')[1] # {seq:cld_num}
    logger(f'Finished loading read cloud numbers from {fasta_prefix + ".R1.fasta"}')

    # Map all reads to assembly graph
    prefix = join(outdir, basename(fastg).split('.')[0])
    r1_clean_df = map_read_clouds(fastg, prefix + '.R1', fasta_prefix + '.R1.fasta', seq_to_cld_num).astype(int)
    r2_clean_df = map_read_clouds(fastg, prefix + '.R2', fasta_prefix + '.R2.fasta', seq_to_cld_num).astype(int)
    clean_df = pd.concat([r1_clean_df, r2_clean_df])
    logger(f'There are {len(clean_df)} read-graph alignments in total')

    # Match overlapping nodes to each other
    node_to_node = {}
    with open(fastg, 'r') as fg:
        for i, l in enumerate(fg): 
            if '>' in l: # Two nodes that overlap
                edge_ids = l.split('_') # Node name at indices 1, 6
                start_node = int(edge_ids[1])
                if start_node not in node_to_node:
                    node_to_node[start_node] = []
                if ':' in l: # Two nodes that overlap
                    node_to_node[start_node].append(int(edge_ids[6]))
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
        aln_nodes = clean_df.Node.unique().tolist()
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

