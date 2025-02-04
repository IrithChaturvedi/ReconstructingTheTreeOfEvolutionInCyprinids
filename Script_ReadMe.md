ReadMe.md

Have Pandas and Matplotlib installed in the working directory
supply all functions returning scaffolds to a desired file

FUNCTIONS:
All files are accessed by opening them in read ('r') mode,
and passing the .readlines() function as the parameter to the user defined function

-> Tackling Intersected Hit files:

project_1.intersected_hits_tab_delimited(lines)
    Description:
        For each sequence/scaffolds passed, it delimits the last descriptive field with TAB for easy reading and analysis
    Parameters: 
        lines - <file>.readlines()

project_1.return_exon_scaffolds(lines, ind)
    Description:
        Returns all scaffolds containing exon features at the specified feature index
    Parameters:
        lines - <file>.readlines()
        ind - index as which the features are specified, here 'exon'

project_1.return_gene_scaffolds(lines, ind)
    Description:
        Returns all scaffolds containing gene features at the specified feature index
    Parameters:
        lines - <file>.readlines()
        ind - index as which the features are specified, here 'gene'

project_1.min_max_coverage(lines)
    Description:
        Prints the maximum and minimum coverage value from all supplied scaffolds
    Parameters:
        lines - <file>.readlines()

project_1.make_coverage_dataframe_genes(lines, feature_index, ensembl_id_index)
    Description:
        Returns a dataframe with the following columns
            ensembl_id, scaffold coverage, feature type, scaffold length
                Ensembl_id - The ENSEMBL id of the scaffold feature aligned to
                Scaffold feature coverage - Coverage of the feature on the Scaffold
                Feature type - Type of feature aligned to the scaffold
                Scaffold length - Lentgh of the feature aligned to
        Parameters:
            lines - <file>.readlines()
            feature_index - index at which the features are present
            ensembl_id_index - index value at which the ensembl_ids are present

project_1.make_coverage_dataframe_exons(lines, feature_index, ensembl_id_index)
    Description:
        Returns a dataframe with the following columns
            ensembl_id, scaffold coverage, feature type, scaffold length
                Ensembl_id - The ENSEMBL id of the scaffold feature aligned to
                Scaffold feature coverage - Coverage of the feature on the Scaffold
                Feature type - Type of feature aligned to the scaffold
                Scaffold length - Lentgh of the feature aligned to
        Parameters:
            lines - <file>.readlines()
            feature_index - index at which the features are present
            ensembl_id_index - index value at which the ensembl_ids are present

project_1.gff3_sort_key(lines)
    Description:
        Sorts the .gff3 file to contain only gene and exon features for intersections
    Parameters:
            lines - <file.gff3>.readlines()

project_1.eval_sort_key(lines, eval_index, evalue)
    Description:
        Parses the blast file results to only contain those sequences having an e-value <= 1e-5
    Parameters:
            lines - <makeblast_file>.readlines()

project_1.graph_maker(d1, d2, title)
    # matplotlib required
    Description:
        Using the two dataframes formed using functions make_coverage_dataframe_exons and make_coverage_dataframe_genes
        to form a coverage vs length graph for gene vs exon sequences
    Parameters:
        d1 - dataframe with gene features
        d2 - dataframe with exon features
        title - desired title supplied for the graph

The following two functions will work in tandem as two tasks for removing exon and gene redundancies in intersected hits files
All steps to be done twice, once for each feature

project_1.long_exons_or_genes_sort(lines,sort_feature)
    Description:
        task 1 involves the creation of another file that contains the longest length alignments only for that feature completely. 
    Parameters:
        lines - initial intersected hits file - <file>.readlines() with good coverages.
        sort_feature - feature that needs to be sorted for only the longest alignments
    
project_1.exons_and_gene_not_in_range(lines, output_task1, feature_sort, feature_index)
    Description:
        task 2 involves the comparison of initial intersects file with the file created from task1 to 
        create a file that contains only those features that are non redundant - do not lie in the range of the longest features and
        are the longest feature alignments for that scaffold.
    Parameters:
        lines - initial intersected hits file - <file>.readlines() with good coverages
        output_task1 - file output from task1
        feature_index - index at which the feature is mentioned
        feature_sort - feature that needs to be sorted for only the longest alignments

The result from task2 is combined with the result from task1 to create a file only containing non-redundant scaffolds for that feature
To simplify this method, a Combined function is created :
    project_1.redundant_gene_exon_parse(lines, feature_sort, feature_index)
        Description: 
            Combines both redundancy removal functions mentioned above to give a complete output of all non redundant gene and exon features.
        lines - initial intersected hits file - <file>.readlines() with good coverages
        feature_index - index at which the feature is mentioned
        feature_sort - feature that needs to be sorted for only the longest alignments

project_1.count_unique_gene_exon(lines)
    Description:
        returns a count of the uniques gene and exon features.
    Parameters:
        lines - initial intersected hits file - <file>.readlines()

project_1.unique_gene_and_exon_list(*l1)
    Description:
        creates a list of all unique gene and exons features for all supplied intersected_hits files for n number of species.
    Paramters:
        *l1 - initial intersected hits files - <file>.readlines(), seperated by ','

project_1.return_ids(l)
    Description:
        creates a list of all unique gene and exons features for an intersected_hits file.
    Parameters:
        lines - initial intersected hits file - <file>.readlines()

project_1.good_coverage_range(lines, high_cov, low_cov)
    Description:
        Parses the blast file results to only contain those sequences in the good coverage range
    Parameters:
        lines - intersected hits file, preferably the one matching to lemmon sequences - <file>.readlines()
        high_cov - highest coverage in the chosen range
        low_cov - lowest coverage in the good range. 

