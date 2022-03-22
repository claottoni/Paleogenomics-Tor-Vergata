############################
Analysis of Metagenomic data
############################

Here we analyse shotgun sequence data from dental calculus. After filtering and collapsing the reads, classify the reads of each sample with Kraken2 and estimate the abudances with Bracken. 
After that, you can use R to run analyses and create charts from the abundance tables generated with Bracken.


************************************************
Prepare abundance tables for downstream analyses
************************************************


Parse the kraken/bracken outputs to abundance tables
****************************************************

Run the ``brackenToAbundanceTable_v2.R`` R script to merge the species aundances of all the samples in one table. The script needs as argument the path to the folder containing the bracken results. 
Move to the bracken results folder and type the following command: 
::
  
  brackenToAbundanceTable_v2.R . 

The script will generate the abundance tables, ``taxa_abundance_bracken_names_IDs.txt``, which contains the species names and the NCBI IDs. 


.. note::
  You can generate abundance tables for various datasets (with bracken files contained in dedicated folders) in different folders and rename the final tables based on the folder (i.e. dataset) name. Try the following commands in the the folder containing all the datasets-subfolders:
  ::
  
    for i in $(find -name "*names_IDs.txt" -type f | sort); do tag=$(basename $(dirname "$i")); cp $i ${FOLDER}/${tag}_abundance_bracken_names_IDs.txt; done

Then use the ``abundanceTablesMerger_v2.R`` script to merge all the tables to obtain one comprehensive table (``abundance_table_IDs.merged``) that you will use for downstream analyses. Remember to rename each file based on the specific dataset before pooling them: 
::

  abundanceTablesMerger_v2.R *.txt


Normalization of the abundance table for genome length
******************************************************

Since the Kraken database was built from complete Bacterial, Archaeal and Viral genomes, we must make a normalization for the genome leghts of each species. This normalization is important when we want to analyse relative taxa abundances to characterise full microbiomes. 
To do that will use a Python `script`_ that takes three arguments: 

  .. _script: https://github.com/claottoni/toolbox/tree/main/gL-nomalizer


1) The abundance table with names
2) The table of genome lengths (parsed from the NCBI genome browser_)
3) The name of the output file

  .. _browser: https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/

The script searches the **names** of the species in the file with the list of genome lengths of species in the NCBI and divides the abundance of each species by the length of the assembly available in the NCBI.
Run the normalization with ``gL-normalizer-lite_v2.py`` and type the command: 
::

  python gL-normalizer-lite_v2.py abundance_table_IDs.merged prokaryotes_viruses_organelles.table abundance_table_IDs.merged.norm

The output table reports the **Species**, the **NCBI IDs**, the **Assembly stage** of the genome in NCBI, the **Length in Megabases**, the kind of **Match**, exact, species (when a peculiar strain is not found in the list of lenghts), and genus (when a species is not found).


Retrieve full taxonomic ranks of the abundance table
****************************************************

In the next step, you will retrieve the full taxonomic ranks for the species of your abundance table. To do that you must have the program Taxaranks_ installed.

  .. _Taxaranks: https://github.com/linzhi2013/taxonomy_ranks

Run the following script, which incorporates Taxaranks, to attach the full taxonomy to every species ID. 
::

  getFullTaxaranks.sh -i abundance_table_IDs.merged.norm

The script generates the following files: 

1) a ``file.taxonomy`` - the full taxonomic ranks (from Kingdom) of the species entries.
2) a ``file.taxonomy.err`` - the species for which no taxonomic ranks could be found.
3) a ``file.taxonomy.final`` - full taxonomy ranks are included in the original abundance table. 

.. note::
  The ``file.taxonomy.err`` should be empty in principle. If not, then the ful taxonomic ranks could not be retrieved for the species listed. In that case, try to update the NCBI taxonomy database as reported in the Taxaranks website.

.. note::
  Depending on the composition of the Kraken DB that you used, you may want to remove unwanted taxa. To do that use ``grep``, for example to remove Eukaryota and Viruses:
  ::
  
    grep -v "Viruses\|Eukaryota\|Fungi" abundance_table_IDs.merged.norm.taxonomy.final > abundance_table_IDs.merged.norm.final.Archaea_Bacteria
  
  Then, you can generate a new table including only species names (without Eukaryota and Viruses)
  ::
      
    cut -f 20- abundance_table_IDs.merged.norm.final.Archaea_Bacteria > abundance_table_IDs.merged.norm.final.Archaea_Bacteria.species


*************
Sourcetracker
*************

First you must prepare the sourcetracker table, in which the NCBI taxonomic IDs and the abundance values (as integers) are used.
You can parse the table from the genome length normalized abudance table with the R script ``makeIntegers.R``, which will generate a ``filename.int`` 
::

  Rscript makeIntegers.R abundance_table_IDs.merged.norm

Then you can add the header typical of the sourcetracker format by tweaking the file: 
::

  sed '1s/^/#Constructed manually\n#OTU id\t/' abundance_table_IDs.merged.norm.int > abundance_table_IDs.merged.norm.int.sourcetracker

After preparing the abundance table, prepare the **mapping file** according to the required `layout`_. 
In the **SourceSink** column, mark as ``sink`` the input samples to be tested, and as ``source`` the reference microbiome samples (e.g. soil, skin, laboratory controls microbiomes).

  .. _layout: https://github.com/danknights/sourcetracker/blob/master/data/metadata.txt

`Sourcetracker`_ is a tool used to predict the source of microbial communities in a set of input samples (i.e., the sink samples). More info in the `publication`_. 

  .. _Sourcetracker: https://github.com/danknights/sourcetracker
  .. _publication: https://www.nature.com/articles/nmeth.1650

Run sourcetracker using the following command:
::

  Rscript sourcetracker_for_qiime.r -i abundance_table_IDs.merged.norm.int.sourcetracker -m mapping_file.txt -o output_directory -v

.. warning::
  The names of the samples in the mapping file (first column) must match those in the abundance table!  



**********************
Analyses with Phyloseq
**********************

.. highlight:: r


We will use the R package ``Phyloseq`` to explore the microbiome data in the normalized abundance table. Once imported in R, the table can be tweaked to focus on a subset of samples. 
Install Phyloseq from `Bioconductor`_: 

  .. _Bioconductor: https://bioconductor.org/packages/release/bioc/html/phyloseq.html


Importing abundance tables in Phyloseq
***************************************

Activate the following libraries: 
::

  library(phyloseq)
  library(ggplot2)
  library(gridExtra)
  library(vegan)

Phyloseq works with ``biom`` objects that store all the required information: abundance of taxa, sample metadata, taxonomy
First you can import the otu table with taxa abundances, and tweak it to remove unwanted columns (the first three resulting from the genome length normalization). 
::

  otu = as.matrix(read.delim("abundance_table_IDs.merged.norm", header=T, fill=T, row.names=1, sep="\t"))
  otu.final = (otu[,-c(1:3)])
  # make the matrix values as numeric
  class(otu.final) <- "numeric"


Then you can import the taxonomy file generated by Taxaranks and the ``getFullTaxaranks.sh`` script (above) and tweak it to remove unwanted columns (1=user_taxa, and 19=species, which is already reported in column 2)
::

  taxonomy = as.matrix(read.delim("abundance_table_IDs.merged.taxonomy", header=T, fill=T, row.names=1, sep="\t"))
  taxonomy.final = (taxonomy[,-c(1,19)])

Finally generate the biom object with phyloseq. 
::

  library(phyloseq)
  # assing otu and tax tables
  OTU = otu_table(otu.final, taxa_are_rows = TRUE)
  TAX = tax_table(taxonomy.final)
  my_biom = phyloseq(OTU, TAX)
  # import mapping file with metadata
  biom_metadata <- import_qiime_sample_data("map_file_full_dataset_Mar2022.txt")
  # merge data
  my_biom <- merge_phyloseq(my_biom, biom_metadata)

Here some commands to explore the biom object: 
::

  # check sample names:
  sample_data(my_biom)
  # check taxa:
  tax_table(my_biom)
  # check rank names:
  colnames(tax_table(my_biom))

If needed, it is convenient to rename the taxa columns more appropriately (instead of Rank1 etc.)
::

  colnames(tax_table(my_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  head(tax_table(my_biom))

Remove low-abundance taxa below a desired threshold (most commonly 0.02% is chosen). 
::

  minTotRelAbun = 0.0002
  keepTaxa = taxa_names(my_biom)[which((x / sum(x)) > minTotRelAbun)]
  my_biom_flt = prune_taxa(keepTaxa, my_biom)

Next, convert the absolute abundance  values in relative abundance. 
::

  my_biom_rel = transform_sample_counts(my_biom, function(x) x / sum(x))


Barplot of taxa abundances
**************************


Multidimensional Scaling (Principal Coordinate Analysis)
********************************************************

Differential taxonomic abundances with DESeq2
*********************************************

