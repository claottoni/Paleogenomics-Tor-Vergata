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

`Sourcetracker`_ is a tool used to predict the source of microbial communities in a set of input samples (i.e., the sink samples). More info in the `publication`_. 

  .. _Sourcetracker: https://github.com/danknights/sourcetracker
  .. _publication: https://www.nature.com/articles/nmeth.1650

First you must prepare the sourcetracker table, in which the NCBI taxonomic IDs and the abundance values (as integers) are used.
You can parse the table from the genome length normalized abudance table with the R script ``makeIntegers.R``, which will generate a ``filename.int`` 
::

  Rscript makeIntegers.R abundance_table_IDs.merged.norm

Then you can add the header typical of the sourcetracker format by tweaking the file: 
::

  sed '1s/^/#Constructed manually\n#OTU id\t/' abundance_table_IDs.merged.norm.int > abundance_table_IDs.merged.norm.int.sourcetracker

After preparing the abundance table, prepare the **mapping file** according to the required `layout`_. 
In the **SourceSink** column, mark as ``sink`` the input samples to be tested, and as ``source`` the reference microbiome samples (e.g. soil, skin, laboratory controls microbiomes).
Here is a `dataset`_ of microbiomes with metadata and a template mapping file to be used for the analysis in our lab.

  .. _layout: https://github.com/danknights/sourcetracker/blob/master/data/metadata.txt
  .. _dataset: https://docs.google.com/spreadsheets/d/18DmeiKz0nlG3TiqXciVZLiPUREh_BhlE/edit?usp=sharing&ouid=116430744637066257064&rtpof=true&sd=true

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
Install Phyloseq from `Bioconductor`_.

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

Then you can generate the biom object with phyloseq. 
::

  library(phyloseq)
  # assing otu and tax tables
  OTU = otu_table(otu.final, taxa_are_rows = TRUE)
  TAX = tax_table(taxonomy.final)
  my_biom = phyloseq(OTU, TAX)

As additional step, you can add to the biom file metadata that will be useful to describe the samples and select them in the downstream analyses (e.g. to work on subset of samples).
The metadata are stored in a mapping file, a tab-delimited text file where you have in the first column the names of your samples (matching those in your OTU abundance table), and the other columns as many metadata as you want (e.g. sample source, literature ID etc.). 
In the dataset of microbiomes mentioned about you can find metadata to be used for the analysis: `Database_metagenomics.xlsx`_

  .. _Database_metagenomics.xlsx: https://docs.google.com/spreadsheets/d/18DmeiKz0nlG3TiqXciVZLiPUREh_BhlE/edit?usp=sharing&ouid=116430744637066257064&rtpof=true&sd=true
 
The metadata can be incorporated in your biom file as follows:
::

  # import mapping file with metadata
  biom_metadata <- import_qiime_sample_data("map_file_full_dataset_Mar2022.txt")
  # merge data
  my_biom <- merge_phyloseq(my_biom, biom_metadata)

Finally, here some commands to explore the biom object: 
::

  # check sample names:
  sample_data(my_biom)
  # check taxa:
  tax_table(my_biom)
  # check rank names:
  colnames(tax_table(my_biom))

.. note::
  If needed, it is convenient to rename the taxa columns more appropriately (instead of Rank1 etc.)
  ::
  
    # check the labels used for rank names
    head(tax_table(my_biom))
    # change them if needed
    colnames(tax_table(my_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

Remove low-abundance taxa below a desired threshold (most commonly 0.02% is chosen). 
::

  minTotRelAbun = 0.0002
  x = taxa_sums(my_biom)
  keepTaxa = taxa_names(my_biom)[which((x / sum(x)) > minTotRelAbun)]
  my_biom_flt = prune_taxa(keepTaxa, my_biom)

Next, convert the absolute abundance  values in relative abundance. 
::

  my_biom_rel = transform_sample_counts(my_biom, function(x) x / sum(x))


Barplot of taxa abundances
**************************
We can generate barplots to display the relative abundance of each sample at different taxonomic ranks. 
To do that we will first merge the abundance of the taxa of our table (e.g. species) to a higher taxonomic rank (e.g. Phylum). 
::
  
  # merge taxa at the Phylum rank with tax_glom (in Phyloseq)
  my_biom_rel_phylum = tax_glom(my_biom_rel, "Phylum")
  # plot the relative abundances (you can use ggplot2 themes to improve the chart)
  plot_bar(my_biom_rel_phylum, fill = "Phylum") + theme(legend.position="bottom", axis.text.x = element_text(size = 4, angle = 90, hjust = 1))

.. warning::
  In the ``tax_glom`` command the label of the rank must match the one that is defined in your original biom file, as reported in the column names of the tax_table command (see above). 

Try to generate barplots at different taxonomic ranks (e.g. family). 
To save the plot in the pdf format you can use this command: 
::

  dev.print(pdf, 'filename.pdf)')

You can also focus your analysis on a subset of samples and use the metadata contained in your original mapping file. It would be convenient to add metadata the way it more convenient to select and compare your samples. Here we subsample for some groups of individuals as defined in the metadata column **Group2**.
::

  # subsample single group phyla
  subsample = subset_samples(my_biom_rel_phylum, Group2=="Portugal calculus")
  # subsample multiple groups phyla
  subgroup = c("Portugal calculus","Espinoso","NTC")
  subsample = subset_samples(my_biom_rel_phylum, Group2 %in% subgroup)
  # barplots in groups with facet_grid
  plot_bar(subsample, fill="phylum") + facet_grid(~Group2, scales= "free_x", space = "free_x") + theme(legend.position="bottom", axis.text.x = element_text(size = 6, angle = 90, hjust = 1))


.. _PCoA:

Principal Coordinate Analysis (Multidimensional Scaling)
********************************************************
The Principal Coordinate Analysis (PCoA), also referred to as metric Multidimensional Scaling (MDS) is a multivariate reduction method performed on distance (or dissimilarity) indexes. Here we will use the Bray-Curtis dissimilarity.
Read more about the MDS and so-called `ordination methods`_

  .. _ordination methods: https://ourcodingclub.github.io/tutorials/ordination/

To make the MDS you can focus on a subset of samples in your species abundance table and run the following commands: 
::
  
  # subsample multiple groups phyla
  subgroup = c("Portugal calculus","Espinoso","Prehistoric human","NTC","Soil","Skin","Gut")
  subsample = subset_samples(my_biom_rel, Group2 %in% subgroup)
  # calculate Bray-Curtis (BC) dissimilarities
  distBC = distance(subsample, method = "bray")
  # make PCoA ordination (MDS) with BC dissimilarities
  ordBC = ordinate(subsample, method = "PCoA", distance = distBC)
  # two-dimension plot PCoA (MDS) setting the colors of the points automatically as defined in the metadata "Group2".
  plot_ordination(subsample, ordBC, color = "Group2")
 
  # you could play with more plotting options and customizations with ggplot2
  plot_ordination(my_biom_rel, ordBC, color = "Group2") + geom_point(mapping = aes(size = xxx, shape = factor(yyy))) + ggtitle("PCoA: Bray-Curtis")

You can repeat the analysis at a higher taxonomic ranks (e.g. genus) and see the difference in the plot. 


.. _nMDS:

Non-Metric Multidimensional Scaling
***********************************
Unlike other ordination techniques that rely on (primarily Euclidean) distances, such as Principal Coordinates Analysis, nMDS uses rank orders (so not the abundances).
The use of ranks omits some of the issues associated with using absolute distance (e.g., sensitivity to transformation), and as a result is much more flexible technique that accepts a variety of types of data. (It’s also where the “non-metric” part of the name comes from.) 
To begin, NMDS requires a distance matrix, or a matrix of dissimilarities. Raw Euclidean distances are not ideal for this purpose: they’re sensitive to total abundances. Consequently, ecologists use the Bray-Curtis dissimilarity calculation.
NMDS arranges points to maximize rank-order correlation between real-world distance and ordination space distance. 
You can read more on the nMDS in the blogs curated by ecologists: `link1`_ `link2`_

  .. _link1: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
  .. _link2: https://archetypalecology.wordpress.com/2018/02/18/non-metric-multidimensional-scaling-nmds-what-how/

The technique uses a *trial and error* to find the best positioning in dimensional space. Goodness-of-fit is measured by **stress** – a measure of rank-order disagreement between observed and fitted distances. 
You can follow the same steps as above, just change the ordinatiom method supported in the command. 
::

  # subsample multiple groups phyla
  subgroup = c("Portugal calculus","Espinoso","Prehistoric human","NTC","Soil","Skin","Gut")
  subsample = subset_samples(my_biom_rel, Group2 %in% subgroup)
  # calculate Bray-Curtis (BC) dissimilarities
  distBC = distance(subsample, method = "bray")
  # make PCoA ordination (MDS) with BC dissimilarities
  ordBC = ordinate(subsample, method = "NMDS", distance = distBC)
  # two-dimension plot PCoA (MDS) setting the colors of the points automatically as defined in the metadata "Group2".
  plot_ordination(subsample, ordBC, color = "Group2")
 



Differential taxonomic abundances with DESeq2
*********************************************
A common goal of microbiome studies is to identify differentially abundant taxa (species, genera etc.) between different groups of samples.
One of the tools used to do that is DESeq2, which was originally developed to identify differentially expressed genes in RNAseq data, but it commonly adopted also in microbiome studies. You can read more about DESeq2 in the `original publication`_ and in the `dedicated page`_ of Bioconductor.
You can install `DESeq2 from Bioconductor`_. 

  .. _original publication: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
  .. _dedicated page: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  .. _DESeq2 from Bioconductor: https://bioconductor.org/packages/release/bioc/html/phyloseq.html

Run DESeq2 on the raw abundace data (`here is why`_), filtered for the low-abundance taxa . You can choose to work on a subset of sample, as done above.

  .. _here is why: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

::

  library(DESeq2)
  # subsample multiple groups phyla
  subgroup = c("Portugal calculus","Prehistoric human","NTC","Soil","Skin","Gut")
  subsample = subset_samples(my_biom_flt, Group2 %in% subgroup)

DESeq2 does not tolerate zero-counts, for this reason an offset (+1) is commonly applied to remove all zeroes from your table.
After that, run DESeq2 on your table and contrat two sets of samples to find differential taxa abundances. 
::

  # make a +1 offset to run Deseq2 (which does not tolerate zero counts)
  otu_table(subsample) = otu_table(subsample)+1  
  # make deseq object from phyloseq object
  ds = phyloseq_to_deseq2(subsample, ~ Group2)
  # Run DESeq2
  dds.data = DESeq(ds)
  # With the 'contrast' function you screen two different set of samples (based on your metadata) for differential taxa. 
  res = results(dds.data, contrast=c("Group2","Portugal calculus","NTC"))
  
Then you can explore your results, filter and sort the differential taxa detected based on a False Dicovery Rate (FDR) threshold (e.g. set to 0.01), the log2-fold-change and the base mean. Read more about the `FDR threshold here`_.

  .. _FDR threshold here: https://www.biostars.org/p/209118/
  
::

  # sort based on p-value adjusted:
  resOrdered = res[order(res$padj, na.last=NA), ]
  # set a threshold value for the False Discovery Rate (FDR):
  alpha = 0.01
  # get only significant taxa based on p-value adjusted (the FDR):
  resSig <- subset(resOrdered, padj < alpha)
  # sort significant values based on the log-fold-change:
  resfc = resSig[order(resSig$log2FoldChange),]
  # sort significant values based on abundance:
  resbm = resSig[order(resSig$baseMean),]
  # save the the tables 
  write.table(as.data.frame(resbm), file="/path/table.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=F)
  
.. note::
  A positive log2-fold-change for a comparison of A vs B means that the OTU in A is more abundant than in B (and viceversa).
  For example, a log2-fold-change of −1 means that in A the OTU is of 2^−1 = 0.5 less abundant than in B.


You can plot the log2-fold-changes to visualize for example the number of differential species belonging to different genera (or phyla) in the two groups contrasted, highlitingh the phylum.  
::

  # generate a data frame to handle the data
  resSigMod = cbind(as(resSig, "data.frame"), as(tax_table(subsample)[rownames(resSig), ], "matrix"))
  # Plot log-fold-changes of otus based on Genus
  ggplot(resSigMod, aes(x=genus, y=log2FoldChange, color=phylum)) +
      geom_jitter(size=3, width = 0.2) +
      theme(axis.text.x = element_text(size = 5, angle = -90, hjust = 0, vjust=0.5))
  # Plot log-fold-changes of tus based on Phylum
  ggplot(resSigMod, aes(x=phylum, y=log2FoldChange, color=phylum)) +
      geom_jitter(size=3, width = 0.2) +
      theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5))

You can also generate a heatmap of the transformed abundance data. For visualization or clustering purposes it might be useful to work with transformed versions of the count data.
Various `transformations`_ may be applied on the data: 
1) ntd: Log transformation of the the data
2) vsd: Variance Stabilizing Transformation 
3) rld: Regularized log transformation 

  .. _transformations: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization

When dealing with many samples it is useful to use the vst. To make the various transformatinos run the folowing commands: 
::

  # variance stabilizing transformation
  vsd = varianceStabilizingTransformation(dds.data, blind=FALSE)
  # log transformation
  ntd <- normTransform(dds.data)
  # regularized log transformation 
  rld <- rlog(dds.data, blind=FALSE)

We will plot the heatmap of the transformed data only for a reduced number of taxa (with decreasing counts).
::

  library("pheatmap")
  #adjust max number of taxa to display.
  select <- order(rowMeans(counts(dds.data,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
  # create a dataframe importing using the metadata (e.g. Group2) and incorporating the taxa names as row names.
  df <- as.data.frame(colData(dds.data)[,"Group2"])
  rownames(df) <- colnames(subsample@otu_table) 
  # plot heatmap (on ntd transformation)
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=FALSE, annotation_col=df, fontsize=4)


Principal Component Analysis (PCA)
**********************************
The Principal Component Analysis (PCA) is another multivariate reduction method used for data visualization. You can read more about the general principles of PCA here: `link3`_, `link4`_

  .. _link3: https://builtin.com/data-science/step-step-explanation-principal-component-analysis
  .. _link4: https://archetypalecology.wordpress.com/2018/02/17/principal-component-analysis-in-r/
  
PCA is performed on the abundance data, however an important concept to keep in mind is that **microbiome datasets are compositional** because they have an arbitrary total imposed by the sequencing instrument. This is due to the fact that high-throughput sequencing (HTS) instruments can deliver reads only up to the capacity of the instrument.
Thus, the total read count observed in a HTS run is a fixed-size, random sample of the relative abundance of the molecules in the underlying ecosystem (e.g. oral microbiota). 
Several methods applied include count-based strategies (normalized to a constant area or volume, e.g. the sequencing lane output) such as Bray-Curtis dissimilarity, zero-inflated Gaussian models and negative binomial models, but these do not account for the limitations imposed by the instrument and the so-called principle of true indpendence of species in ecology. 
Read more about compositonal data in the paper by `Gloor et al. (2017)`_.

  .. _Gloor et al. (2017): https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
  
Due to the compositional nature of microbiome datasets, the centered log-ratio (clr) transformation introduced by Aitchison (1986) is often used. The clr-transformed values are scale-invariant, which means that the same ratio is expected to be obtained in a sample with few read counts or an identical sample with many read counts, only the precision of the clr estimate is affected. 
The clr transformation uses the geometric mean of the sample vector as the reference. In particular, sample vectors undergo a transformation based on the logarithm of the ratio between the individual elements and the geometric mean of the vector: 

.. math::

  g(x) = \sqrt[D]{x_{i}, ...,  x_{D}}

Log-ratio transformation is applied to each subject vector *x* with *D* features (OTUs). Basically, the analysis of clr data will reveal how OTUs behave relative to the per-sample average.

.. math::

  clr(x) = [\ln{\frac {x_{i}}{g(x)}},...,\ln{\frac {x_{D}}{g(x)}}]


We will accomplish the PCA of oral microbiota only (while for comparing different microbiota we used the :ref:`PCoA<PCoA>`, or the :ref:`nMDS<nMDS>`). 
For this reason we will subsample the original biom abundance table exclusively for oral groups, for example: 
::

  subgroup = c("Human plaque","Chimpanzee","Prehistoric human","Historic human","Modern human")
  subsample = subset_samples(my_biom, Group2 %in% subgroup)

Then, we filter out low-abundance taxa with a mean < 0.02%. 
::

  minTotRelAbun = 0.0002
  x = taxa_sums(subsample)
  keepTaxa = taxa_names(subsample)[which((x / sum(x)) > minTotRelAbun)]
  subsample_flt = prune_taxa(keepTaxa, subsample)

.. warning::
  The removal of low-abundance taxa is performed only **after** subsampling the abudance table. Doing this before will keep unwanted taxa that are more common in other microbiota (e.g. soil) and irrelevant for describing oral envrionments (or the specific environment you are analysing with the PCA).  

We will do the clr-transformation of the otu table in the biom object with the package `microbiome`. To install the package: 
::

  library(BiocManager)
  BiocManager::install("microbiome")

Then activate the library and transform the data with the following command: 
::
  
  library(microbiome)
  subsample_flt_clr = transform(subsample_flt,"clr")

To refine the analysis you can remove unwanted samples (e.g. those that appeared to be contaminated from Sourcetracker). 
::

  remove = c("VK1","CS45","CS47","CS48","Jomon_5")
  subsample_pruned = prune_samples(!(subsample_flt_clr@sam_data$Sample_short %in% remove), subsample_flt_clr)

Finally, make the PCA with the command `ordinate`, adopting the ordination method `RDA` (redundancy analysis), which without constraints correspond to the PCA in phyloseq. 
::

  ord_clr = ordinate(subsample_pruned, "RDA")

You can plot the variance associated with each PC: 
::
  
  plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
  
And finally plot the PCA in a 2D chart. 
::

  plot_ordination(subsample_flt_clr, ord_clr, color="Group2") + 
    geom_point(size = 2)
 
To display the labels of each symbol: 
::

  plot_ordination(subsample_flt_clr, ord_clr, color="Group2", label="Sample_short") + 
    geom_point(size = 2)

**ALTERNATIVE:** Another package used to perform the PCA is `mixOmics`. To install it: 
::

  ## install BiocManager if not installed 
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  ## install mixOmics 
  BiocManager::install('mixOmics')

Here are the commands to make the PCA with mixOmics. Note the conversion to count-per-million and the +1 offest that is necessary to remove zero-counts (which cannot be handled by the clr-transformation of mixOmics).
::

  # remove contaminated samples
  remove = c("VK1","CS45","CS47","CS48","Jomon_5")
  subsample_pruned = prune_samples(!(subsample_flt@sam_data$Sample_short %in% remove), subsample_flt)
  # isolate the otu_table of the biom as separate object. 
  tab = otu_table(subsample_pruned)
  # make on offset (+1) of your count data, after converting them to counts-per-million (to make the offset irrelevant) 
  tab_cpm_off = (tab*1000000)+1
  tab_cpm_off = t(tab_cpm_off)
  # create e vector from the metadata that you want to use:
  group = subsample_pruned@sam_data$Group2
  # make the PCA with mixOmics. Note the CLR tranformation called in the command. 
  library(mixOmics)
  tune.pca(tab_cpm_off, ncomp = 10, center = TRUE, scale = FALSE, logratio = 'CLR')
  pca.res <- pca(tab_cpm_off, ncomp = 10, center = TRUE, scale = FALSE, logratio = 'CLR')
  # PLot PCA 
  par(mfrow = c(1, 1))
  plot(pca.res$variates$X[,1], pca.res$variates$X[,2], type="n", main="PCA", cex.axis=0.75, cex.lab=0.75, xlab="", ylab="")
  abline(v=0, lty=3, col="grey")
  abline(h=0, lty=3, col="grey")
  # add the points and the symbols (pch) based on the number of groups present in your metadata (those selected when subsampling the oral microbiota)
  points(pca.res$variates$X[,1], pca.res$variates$X[,2], cex=1.2, bg=factor(group), col=factor(group), pch = c(21,22,23,3,4)[as.factor(group)], lwd=0.7)      
  legend("topleft", legend = sort(unique(group)), bty = "n", col = sort(unique(as.factor(group))), pt.cex=1, cex=0.5, pt.bg=sort(unique(as.factor(group))), pch = c(21,22,23,3,4), pt.lwd=0.6)
  text(pca.res$variates$X, labels = subsample_pruned@sam_data$Sample_short, cex=0.4)
