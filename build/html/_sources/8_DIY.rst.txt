##############
DO-IT-YOURSELF
##############


*************************************
Hands-on 1: ancient human mtDNA
*************************************

In this final hands-on session you will analyse shotgun (reduced) sequencing data generated from an ancient human tooth.The genomic library built from the DNA extract was sequenced on an Illumina platform in paired-end mode. Your task is:   

1. Process the raw reads (remove adapters, merge the reads, section 2). 
2. Align the reads to the human mitochondrial DNA (mtDNA) reference sequence, assess the damage of DNA molecules, call the variants (sections 4-5-6).  
3. Run the metagenomic screning of the DNA extract with Kraken using the Minikraken database (section 3).

After reads pre-processing it is up to you whether first aliging the reads or screening the metagenomic content. 

- **Option 1**: You can use your ``vcf`` file to assign an haplogroup to the human samples that you analysed. Some useful tools for haplogroup assignation:  
  
    - Check the variant positions in Phylotree (http://www.phylotree.org/)  
    - Load the ``vcf`` file in Haplogrep (https://haplogrep.uibk.ac.at)
   
- **Option 2**: Run again the metagenomic screening with a Custom Database of Kraken (provided by us), and compare the results with those obtained with Minikraken.




*****************************************
Hands-on 2: ancient human dental calculus
*****************************************

In this hands-on you will use a sequence dataset of ancient dental calculus from a recent study by `Velsko et al.`_ (Microbiome, 2019), which consists of paired-end sequence data. 
After filtering and collapsing the reads, classify the reads of each sample with Kraken2 and estimate the abudances with Bracken. 

  .. _Velsko et al.: https://link.springer.com/article/10.1186/s40168-019-0717-3

After that, you can use R to run analyses and create charts from the abundance tables generated with Bracken.


Abundance table and genome length normalization
***********************************************

Run the following R script to merge the species aundances of all the samples in one table. The script needs as argument the path to the folder containing the bracken results. 
Move to the bracken results folder and type the following command: 
::
  
  Rscript brackenToAbundanceTable.R . 

The script will generate two abundance tables, ``taxa_abundance_bracken_IDs.txt``, which contains the species names as NCBI IDs, and ``taxa_abundance_bracken_names.txt``, which contains the actual species names. 
Since the Minikraken database was built from complete Bacterial, Archaeal and Viral genomes, we must make a normalization for the genome leghts of each species. This normalization is important when we want to analyse relative taxa abundances to characterise full microbiomes. 
To do that will use a Python script that takes three arguments: 

1) The abundance table
2) The table of genome legths
3) The name of the output file

To run the normalization type the command: 
::

  python gL-normalizer-lite.py taxa_abundance_bracken_names.txt prokaryotes_viruses_organelles.table taxa_abundance_bracken_names_normalized.txt


R session
*********

Import files and prepare the abundance tables
=============================================

Download the normalized table ``taxa_abundance_bracken_names_normalized.txt`` in your local machine (e.g. with ``scp``) and open R. Feel free to rename the file at your own convenience (here we will rename it ``mytable_normalized.txt``)

.. highlight:: r

First of all, in R, we must set up the folder which contains the abundance table:   
::

  setwd("/path/to/your/folder")
    
Then we can import the abundance table files, setting the 1st column as row names. Then we remove the 1st, 2nd and 3rd column, which will not be used for the analysis. 
::

  mytable = read.delim("mytable_normalized.txt", header=T, fill=T, row.names=1, sep="\t")
  mytable.final = mytable[,-c(1:3)]

Then we filter the species in the table for their abundance by removing those that are represented below a threshold of 0.02%. 
To do that, we define a function (that we call low.count.removal)
::

  low.count.removal = function(
                          data, 		# OTU count data frame of size n (sample) x p (OTU)
                          percent=0.02	# cutoff chosen
                          ){
      keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
      data.filter = data[,keep.otu]
      return(list(data.filter = data.filter, keep.otu = keep.otu))
  }

We run the function on our table, setting up the threshold at 0.02%. Note that the table is transposed to run the function: 
::

  result.filter = low.count.removal(t(mytable.final), percent=0.02)

We generate a table with the filtered data: 
::

  mytable.final.flt = result.filter$data.filter
  
In the next step, we normalize the reads for the sequencing depth. This means that we will account for the total reads generated for each sample, and normalize the species abundance for that number. 
To do that we will define a **Total Sum Squared** function 
::

  TSS.divide = function(x){
   x/sum(x)
  }

The function is applied to the table, and each row must represent a sample. For this reason we transpose the table.
::

  mytable.final.flt.tss = t(apply(mytable.final.flt, 1, TSS.divide))

We have just generated a table of species abundances of ancient dental calculus samples, normalized for genome lenghts and sequencing depths. 
For comparative analysis, we can now include in our analysis a dataset of normalized species abundances generated with Minikraken (version ``minikraken2_v1_8GB_201904``) representing other microbiomes.  
To do that, we wil repeat all the steps described above. Note that there is no need to define again the functions created above because they are stored in current R session environment. 
::

  setwd("/path/to/table")
  table.lit = read.delim("taxa_abundance_bracken_names_normalized_literature.txt", header=T, fill=T, row.names=1, sep="\t")
  table.lit.final = table.lit[,-c(1:3)]  
  result.filter = low.count.removal(t(table.lit.final), percent=0.02)
  table.lit.final.flt = result.filter$data.filter
  table.lit.final.flt.tss = t(apply(table.lit.final.flt, 1, TSS.divide))

Now we can merge the two tables in one, by merging them for the column containing the species names (this column is selected with ``by=0``)
::

  table.total = merge(t(mytable.final.flt.tss), t(table.lit.final.flt.tss), by=0, all=TRUE)

The following commands are used to finalize the table: 
::

  table.total[is.na(table.total)] <- 0		#removes NA
  row.names(table.total) = table.total[,1]	#copy the species names in 1st columns to row names 
  table.total.final = t(table.total[,-1])	#delete the first column with the species names (now reported as row names)
 
.. note:: 

  You can merge datasets only if they were generated with the same taxonomy database (here the Minikraken 8Gb database). If not, you'll have to run all the samples from the literature with the same database that you used to analyse your samples. 


UPGMA
=====

Once generated the final including both datasets (dental calculus and other microbiomes), we run an UPGMA cluster analysis. We must first install the ``vegan`` and ``ape`` package in R.
::

  install.packages("vegan")		#do it only if the package is not installed yet
  install.packages("ape")		#do it only if the package is not installed yet
  library(vegan)
  library(ape)

Then we use vegan to calculate the **Bray-Curtis** distances, and run the cluster analysis.
::

  bray_dist = vegdist(table.total.final, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
  bray_dist.clust = hclust(bray_dist, method="average", members = NULL)

Finally, we plot the dendrogram: 
::

  plot(as.phylo(bray_dist.clust), type = "unrooted", cex = 0.5, lab4ut="axial", no.margin=T, show.tip.label=T, label.offset=0.02, edge.color = "gray", edge.width = 1, edge.lty = 1)

To visualize better our samples in the following charts, we can define metadata as vectors. We assign group labels on each sample, creating a vector of labels that corresponds to the samples of the dataset that we are analysing. 
For example, for the literature samples, we generate the following vector of metadata describing the kind of sample.
::

  labels_lit = c("Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth",
					"Soil","Soil","Soil","Soil","Soil","Soil","Soil",
					"Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth","Ancient calculus","Ancient tooth",
					"Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque",
					"Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque","Plaque",
					"Plaque","Plaque","Plaque","Plaque","Plaque",
					"Skin","Skin","Skin","Skin",
					"Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut",
					"Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut",
					"Gut","Gut","Gut","Gut",
					"Skin","Skin","Skin","Skin","Skin",
					"Plaque","Plaque")

Create a vector with labels corresponding to the samples that you analysed. Always make sure to follow the order of the samples if you use more than one label. 
For example, the metadata of a table from Velsko's dataset containing ancient samples in the first five rows, and modern samples in the following five rows, will be represented by this vector: 
::

  labels_Velsko = c("Velsko_ancient","Velsko_ancient","Velsko_ancient","Velsko_ancient","Velsko_ancient",
			"Velsko_modern","Velsko_modern","Velsko_modern","Velsko_modern","Velsko_modern")

Merge the labels, again paying attention to the order that you used to merge the tables (first your samples, then the literature dataset)
::
  
  labels = c(labels_Velsko,labels_lit)
  

To have a better look at the correspondence of data we can create a dataframe: 
::

  table.total.final.df = as.data.frame(table.total.final)
  table.total.final.df$group = labels
  
We assign colors to each label. Note that the colors are assigned alphabetically based on the labels that you used.  
::

  coul=c("#E41A1C",		#Ancient calculus		
	"#419681",		#Ancient tooth					
	"#4DAF4A",		#Gut				
	"lightgray",		#Plaque	
	"#984EA3",		#Skin		
	"#FF7F00",		#Soil		
	"goldenrod",		#Velsko-ancient		
	"#994C00")		#Velsko-modern

And finally, we plot the dendrogram by customizing the tips with the color-coded labels:
::

  plot(as.phylo(bray_dist.clust), type = "unrooted", cex = 0.5, lab4ut="axial", no.margin=T, show.tip.label=T, label.offset=0.02, edge.color = "gray", edge.width = 1, edge.lty = 1)
  tiplabels(pch=19, col = coul[factor(labels)], bg = coul[factor(labels)], cex=1, lwd=1)          

And we can add a legend:
::

  legend("topleft", legend = sort(unique(labels)), bty = "n", col = coul, pch = 19, pt.cex=1, cex=0.6, pt.lwd=1)


