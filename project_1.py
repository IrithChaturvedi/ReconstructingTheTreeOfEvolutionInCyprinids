def intersect_sort_key(line):
    fields = line.strip().split("\t")
    fields[-1] = '\t'.join(fields[-1].split(';')) 
    return '\t'.join(fields)

def intersected_hits_tab_delimited(lines):
    sorted_lines = sorted((intersect_sort_key(line) for line in lines), key=lambda x: x.lower())
    return '\n'.join(sorted_lines) + '\n'
    
def return_exon_scaffolds(lines, ind):
    for line in lines:
        fields = line.strip().split("\t")
        if (fields[ind] == 'exon'): 
            return '\t'.join(fields) + '\n'
        else:
            return None

def return_gene_scaffolds(lines, ind):
    for line in lines:
        fields = line.strip().split("\t")
        if (fields[ind] == 'gene'): 
            return '\t'.join(fields) + '\n'
        else:
            return None
        
def min_max_coverage(lines):
    mini = 50
    maxim = 50
    for line in lines:
        fields = line.strip().split("\t")
        id = fields[3].split('_')
        if float(id[5]) > maxim:
            maxim = float(id[5])
        elif float(id[5]) < mini:
            mini = float(id[5])
    print(mini, maxim)

def make_coverage_dataframe_genes(lines, feature_index, ensembl_id_index):
    import pandas as pd
    cov = []
    num = 1
    for line in lines:
        fields = line.strip().split("\t")
        if fields[feature_index] == 'gene':
            id = fields[3].split('_')
            d = {"ensembl_id":fields[ensembl_id_index], "coverage":float(id[5]), "type":fields[feature_index], 'length':float(id[3])}
            num+=1
            cov.append(d)
    return pd.DataFrame(cov)

def make_coverage_dataframe_exons(lines, feature_index, ensembl_id_index):
    import pandas as pd
    cov = []
    num = 1
    for line in lines:
        fields = line.strip().split("\t")
        if fields[feature_index] == 'exon':
            id = fields[3].split('_')
            d = {"ensembl_id":fields[ensembl_id_index], "coverage":float(id[5]), "type":fields[feature_index], 'length':float(id[3])}
            num+=1
            cov.append(d)
    return pd.DataFrame(cov)

def gff3_sort_key(lines):
    for line in lines:
        fields = line.strip().split("\t")
        if ('#' not in fields[0]) and ((fields[2] == 'exon') or (fields[2] == 'gene')):
            return '\t'.join(fields) + '\n'
        else:
            return None
        
def eval_sort_key(lines, eval_index, evalue):
    for line in lines:
        fields = line.strip().split("\t")
        if (float(fields[eval_index]) <= evalue):
            return '\t'.join(fields) + '\n'
        else:
            return None
        
def graph_maker(d1, d2, title):
    import matplotlib.pyplot as plt
    fig, (ax3) = plt.subplots(1)
    fig.set_size_inches(10, 5)
    fig.set_dpi(100)
    ax3.scatter(d1['coverage'], d1['length'], color = 'blue', label = 'Genes')
    ax3.scatter(d2['coverage'], d2['length'], color = 'red', label = 'Exons')
    ax3.set_xlabel("Cov")
    ax3.set_ylabel("Len")
    ax3.set_title(title)
    ax3.legend()
    return plt.show()

def long_exons_or_genes_sort(lines, sort_feature):
    l_seq_max_dist = {}
    sorted_lines = []

    for line in lines:
        fields = line.strip().split("\t")
        if fields[6] == sort_feature:
            l_seq = fields[3]
            dist = int(fields[8]) - int(fields[7])

            if l_seq not in l_seq_max_dist:
                l_seq_max_dist[l_seq] = dist
            else:
                if dist > l_seq_max_dist[l_seq]:
                    l_seq_max_dist[l_seq] = dist

    # Use a list to store all 'exon' features meeting the condition
    selected_exons = []
    for line in lines:
        fields = line.strip().split("\t")
        if fields[6] == sort_feature:
            if (fields[3] in l_seq_max_dist) and (int(fields[8]) - int(fields[7]) == l_seq_max_dist[fields[3]]):
                selected_exons.append('\t'.join(fields) + '\n')

    # Add the selected 'exon' features to the output
    sorted_lines.extend(selected_exons)

    return sorted_lines

def exons_and_gene_not_in_range(lines, output_task1, feature_sort, feature_index):
    exons_in_range = set()

    # Collect exons from lines2 and store their ranges in a set
    for line2 in output_task1:
        fields2 = line2.strip().split("\t")
        if fields2[feature_index] == feature_sort:  # Ensure it's the required feature
            seq_name = fields2[3]
            exon_start = int(fields2[feature_index+1])
            exon_end = int(fields2[feature_index+2])
            exons_in_range.add((seq_name, exon_start, exon_end))

    output_lines = []

    # Process feature from lines, and add them to the output if not in the range
    for line in lines:
        fields = line.strip().split("\t")
        if fields[feature_index] == feature_sort:  # Ensure it's the feature
            seq_name = fields[3]
            exon_start = int(fields[feature_index+1])
            exon_end = int(fields[feature_index+2])

            # Check if the exon does not overlap with any required feature in lines2
            if all(
                exon_start >= end2 or exon_end <= start2
                for (_, start2, end2) in exons_in_range
            ):
                output_lines.append('\t'.join(fields) + '\n')

    return output_lines

def redundant_gene_exon_parse(lines, feature_sort, feature_index):
    output_1 = long_exons_or_genes_sort(lines, feature_sort)
    output_2 =  exons_and_gene_not_in_range(lines, output_1, feature_sort, feature_index)
    return output_1.extend(output_2)

def count_unique_gene_exon(lines):
    c_g = 0
    c_e = 0
    ens = []
    for line in lines:
        fields = line.strip().split("\t")
        if fields[6] == 'gene':
            if fields[12] not in ens:
                ens.append(fields[12])
                c_g+=1
        else:
            if fields[13] not in ens:
                ens.append(fields[13])
                c_e+=1
    print(c_g, c_e)

def unique_gene_and_exon_list(*l1):
    uniq = []
    def process_list(lst):
        for l in lst:
            elements = l.strip().split("\t")
            if elements[6] == 'gene':
                if elements[12] not in uniq:
                    uniq.append(elements[12])
            elif elements[6] == 'exon':
                if elements[17] not in uniq:
                    uniq.append(elements[17])

    for i in l1:
        process_list(i)

    return uniq


def return_ids(l):
    ids = []
    for i in l:
        elements = i.strip().split('\t')
        if elements[6] == 'gene':
            if elements[12] not in ids:
                ids.append(elements[12])
        elif elements[6] == 'exon':
            if elements[13] not in ids:
                ids.append(elements[17])
    return ids

def good_coverage_range(lines, high_cov, low_cov):
    for line in lines:
        fields = line.strip().split("\t")
        id = fields[3].split('_')
        if float(id[5]) <= high_cov and float(id[5]) > low_cov:
            return '\t'.join(fields) + '\n'
        else:
            return None