#######################################
Antimicrobial Resistance (AMR) analysis 
#######################################

Antimicrobial resistance (AMR) is present in host-associated oral microbiomes, and metagenomic analysis of ancient dental calculus recently demonstrated that it occurs in natural environments and in ancient human and animal samples.
For the analysis we will use the `Comprehensive Antibiotic Resistance Database`_ (CARD). The CARD is a rigorously curated collection of characterized, peer-reviewed resistance determinants and associated antibiotics, organized by the Antibiotic Resistance Ontology (ARO) and AMR gene detection models.

  .. _Comprehensive Antibiotic Resistance Database: https://card.mcmaster.ca

***********************
Build the CARD database
***********************

First, you have to download the CARD database. Make a dedicated folder, download the database into it the db and decompress it. The command below downloads the latest version of the database. 

.. code-block:: bash

  mkdir card_aro_db
  cd card_aro_db
  wget -c https://card.mcmaster.ca/latest/data -O card_db.tar.bz2
  tar -xf data.tar.bz2
  
Of the files extracted, the one to use for building the database is the ``nucleotide_fasta_protein_homolog_model.fasta``. This file contains sequences of antimicrobial resistance genes that do not include mutation as a determinant of resistance - these data are appropriate for BLAST analysis of metagenomic data or searches excluding secondary screening for resistance mutations. 
In contrast, the "protein variant" model includes reference wild type sequences used for mapping SNPs conferring antimicrobial resistance - without secondary mutation screening, analyses using these data will include false positives for antibiotic resistant gene variants or mutants.

Then, it is suggested to replace the empty spaces in the file with underscores, so as to keep sequence IDs as single string (hence, the full info about the hits in the following blastn analysis).

.. code-block:: bash

  sed 's/ /_/g' nucleotide_fasta_protein_homolog_model.fasta > nucleotide_fasta_protein_homolog_model.fasta_mod.fasta

Finally, we build the database with `Blast+`_. The database is created inside a folder named ``card_nucl_db``. 
If you work in G100, you can generate the db in the login node, as the operation is quite fast. 

  .. _Blast+: https://www.ncbi.nlm.nih.gov/books/NBK279690/

.. code-block:: bash

  # In Galileo100 Blast+ is available by activating the dedicated modules. 
  module load profile/bioinf
  module load autoload blast+
  # Generate the CARD database. 
  makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta_mod.fasta -dbtype nucl -out card_nucl_db

The db generated consists of seven files all named with the prefix given in the ``-out`` option of the makeblastdb command.

*********************************************
Build the database of housekeeping (HK) genes
*********************************************

To account for potential taphonomic conditions leading to preservation biases across samples of different chronology, we will use blastn to align the reads against the nucleotide sequences of three microbial housekeeping genes (recA, rpoB, gyrB) deposited in the NCBI.
We first download the sequences with the tools ``esearch`` and ``efetch`` of `EDirect`_ and concatenate them in one fasta file that we will use to build the database.

  .. _EDirect: https://www.ncbi.nlm.nih.gov/books/NBK179288/

.. code-block:: bash

  # Download the sequences of the housekeeping genes
  esearch -db nucleotide -query '"recA"[Protein name]' | efetch -format fasta > recA_nu.fasta
  esearch -db nucleotide -query '"rpoB"[Protein name]' | efetch -format fasta > rpoB_nu.fasta
  esearch -db nucleotide -query '"gyrB"[Protein name]' | efetch -format fasta > gyrB_nu.fasta
  cat recA_nu.fasta rpoB_nu.fasta gyrB_nu.fasta >> recA_rpoB_gyrB_nu.fasta
  # Generate the HK database.   
  makeblastdb -in recA_rpoB_gyrB_nu.fasta -dbtype nucl -out recA_rpoB_gyrB_nu

Just like the CARD db, the HK db generated consists of seven files all named with the prefix given in the ``-out`` option of the makeblastdb command.

*****************************************************
Interrogate the CARD and Housekeeping genes databases
*****************************************************

Before interrogating the databases it is important to isolate from the fastq files that we are analysing the reads classified as Bacteria in Kraken2. To do that we will use ``extract_kraken_reads.py`` from `KrakenTools`_. 
The reads are extracted as fasta files. To run the following commands, make sure that the ``.collapsed.gz`` and the ``.krk`` files have the same name. 

  .. _KrakenTools: https://github.com/jenniferlu717/KrakenTools

.. code-block:: bash

  KRK=/path/to/kraken/output
  OUTPUT=/path/to/output_fasta_files
  for i in *.collapsed.gz; 
  do 
    filename=$(basename "$i")
    fname="${filename%.collapsed.gz}"
    extract_kraken_reads.py -k ${KRK}/${fname}.krk -s $i -o ${OUTPUT}/${fname}.bacteria.fasta -t 2 --include-children -r ${KRK}/${fname}.krk.report
  done


Here are some key options of ``extract_kraken_reads.py``:

================================= ========
BCFtools call options             Function
================================= ========
**-k, --kraken**           		  Kraken output file.
**-s, -s1, -1, -U**               FASTA/FASTQ sequence file (may be gzipped).
**-s2, -2**                       FASTA/FASTQ sequence file (for paired reads, may be gzipped).
**-o**                            output FASTA/Q file with extracted seqs.
**--include-children**            include reads classified at more specific levels than specified taxonomy ID levels (optional).
**--include-parents**             include reads classified at all taxonomy levels between root and the specified taxonomy ID levels (optional).
**-r, --report**                  Kraken report file (required if specifying --include-children or --include-parents).
================================= ========


Finally, we interrogate the two databases with ``blastn`` using the fasta files as query sequences. First the CARD db and then the HK db. 
To run these commands, you can make a list of all the fasta files to process in a file that you can call `samples_fasta.txt`, and loop through the fasta files.
It is a recommended to usa a tag so as to distinguish in the filenames the outputs of the two different databases. 

.. code-block:: bash

  # Interrogate the Card database
  DBNAME=/path/to/card_nucl_db
  OUTPUT=/path/to/output/amr_analysis
  TAG=card
  
  for FASTA in $(cat samples_fasta.txt)
  do
    FILENAME=$(basename "$FASTA")
    SAMPLE=${FILENAME%.fasta}
    blastn -query $FASTA -db $DBNAME -out ${OUTPUT}/${SAMPLE}.${TAG}.blastn.out -outfmt 6
  done

  # Interrogate the HK genes database
  DBNAME=/path/to/recA_rpoB_gyrB_nu
  OUTPUT=/path/to/output/amr_analysis
  TAG=hk
  
  for FASTA in $(cat samples_fasta.txt)
  do
    FILENAME=$(basename "$FASTA")
    SAMPLE=${FILENAME%.fasta}
    blastn -query $FASTA -db $DBNAME -out ${OUTPUT}/${SAMPLE}.${TAG}.blastn.out -outfmt 6
  done

The option ``-outfmt 6`` indicates that we are generating BLASTn outputs in tabular format. The default output contains 12 tab-delimited columns: 

  1. qseqid - query or source (gene) sequence id
  2. sseqid - subject or target (reference genome) sequence id
  3. pident - percentage of identical positions
  4. length - alignment length (sequence overlap)
  5. mismatch - number of mismatches
  6. gapopen - number of gap openings
  7. qstart - start of alignment in query
  8. qend - end of alignment in query
  9. sstart - start of alignment in subject
  10. send - end of alignment in subject
  11. evalue - `expected value`_
  12. bitscore - `bit score`_
 
  .. _expected value:  https://www.metagenomics.wiki/tools/blast/evalue
  .. _bit score:  https://www.metagenomics.wiki/tools/blast/evalue#h.4wxezjs2qtog
  
  
*********************
Normalization of data
*********************

First of all we used a custom python script `aro_blastn_parser_v2.py` (available in a dedicated `repository`_ of Github) to parse the results of the blastn analysis against the HK db in one table. The number of hits normalized for the sequencing depth is also reported in the table. 
 
  .. _repository:  https://github.com/claottoni/toolbox
  
.. code-block:: bash

  amr_blastn_parser_v2.py *.blastn.out > hits_hk_genes.txt
  
Then, we retrieve the `Antibiotic Resistance Ontology`_ (ARO) categories for the hits in the card database by using the ``aro_index.tsv`` fila downloaded together with the CARD database. We generate a file ``.index`` for each sample with the custom script ``aro_blastn_parser.py``.
We can use the code below to loop through all the output files of blastn.

  .. _Antibiotic Resistance Ontology:  https://card.mcmaster.ca/aro/list

.. code-block:: bash

  ARO_INDEX=/path/to/card_nucl_db/aro_index.tsv
  for i in $(find -name "*aro.blastn.out" -type f)
  do 
    filename=$(basename "$i")
    fname="${filename%.out}"
    echo "parsing $i"
    aro_blastn_parser.py $i $ARO_INDEX > $(dirname "$i")/${fname}.out.index
  done
  
Multiple hits for each sequence may be reported in the blastn output and they are sorted by decreasing bitscore (last column). We use the following command to sort and uniq the read-names so as to parse a table where only the first read (in case of multiple hits) with the highest bitscore will be returned.
The final table contains two columns, the read name and the ARO category. 

.. code-block:: bash

  for i in $(find -name "*.out.index" -type f); do echo $i; sort -u -k1,1 $i | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$4}' > ${i}.uniq; done
  

**********************
AMR data analysis in R
**********************
The following analyses are done in R. First, we parse the `.uniq` files of all the samples in one comprehensive abundance table reporting the number hits for each sample. 

.. code-block:: r

  files <- list.files(pattern="*.uniq$", full.names=T,recursive=FALSE)
  for (i in files) {
    if (!exists("tabfinal")){
      tab=read.delim(i, header=F, fill=T, row.names=NULL, sep="\t")
      tabfinal = table(tab$V2)
      print(paste0("parsing ", files[c(1,2)]))
      } else {
    # merge all the others
    tab=read.delim(i, header=F, fill=T, row.names=NULL, sep="\t")
    tabfinal = merge(tabfinal, table(tab$V2), by=1, all=T)}
    #print(paste0("parsing ", i)
  }
  colnames(tabfinal) = c("gene", files)
  write.table(tabfinal, file = 'hits_amr_genes.txt', sep="\t", row.names=F, na="0", col.names=T, quote = FALSE)
  
Now that the hits abundance of all the samples against each database are parsed in two tables, we can import them in order to normalize the number of hits in the card database (AMR genes) with the number of hits in the housekeeping genes database. 

.. code-block:: r
  
  # import the tables.
  amr_genes = read.delim("hits_amr_genes.txt", header=T, fill=T, row.names=NULL, sep="\t")
  hk3_genes = read.delim("hits_hk_genes.txt", header=T, fill=T, row.names=NULL, sep="\t")
  # adjust the layout in the amr_genes table.
  amr_genes[is.na(amr_genes)] <- 0
  amr_genes.final = amr_genes		
  row.names(amr_genes.final) = amr_genes.final[,1]
  amr_genes.final = amr_genes.final[,-1]
  # transpose the amr genes table
  amr_genes.final = t(amr_genes.final)

.. warning::
   
   Make sure that the order of the samples in the two tables is exactly the same, if not the results of the normalization will be inconsistent!

Finally we normalize the hits in the AMR genes db with those in HK genes db hits and create a dataframe. Note that the normalized abundance is reported as counts per million (cpm)

.. code-block:: r

  amr_genes.final.norm = amr_genes.final/hk3_genes.sorted$Hits*1000000
  amr_genes.final.norm = as.data.frame(amr_genes.final.norm)
  write.table(amr_genes.final.norm, "dataset_AMR_hk3norm_cpm.tsv", quote=F, sep="\t", row.names=T, col.names=T)
  
We can add metadata such as chronological groups or other group identifiers as columns in the dataframe. You can create a vector with the data listed in the same order as the samples, or you can create a file with the data listed as column(s) (column1=Sample, column2=Group). 

.. code-block:: r
  
  # option 1: create a vector
  amr_genes.final.norm$group = c("group1","group2","group3","etc")
  # option 2: generate a metadata file (tab-delimited)
  metadata = read.delim("metadata.txt", header=T, row.names=1)
  amr_genes.final.norm$group = metadata$Group


To identify ARO gene families significantly different among the groups we use a Wilcoxon Rank Sum Tests with Benjamini & Hochberg adjusted P-values by selecting the columns in the dataframe corresponding to the gene families (in the example here from column 1 to 33).

.. code-block:: r

  # make wilcoxon test
  wilcox = lapply(amr_genes.final.norm[,c(1:33)], function(x) pairwise.wilcox.test(x, amr_genes.final.norm$group, p.adjust.method = "BH"))
  library(plyr)
  out <- ldply(wilcox, function(x) x$p.value)
  # save table
  write.table(out, "dataset_wilcox_test_hk3_norm_full.tsv", quote=F, sep="\t", row.names=F, col.names=T)
  
We can generate a multipanel plot with the results using ggplot. First we'll prepare the data with ``pivot_longer``, which will be used to select the gene families to plot (in the example below we plot all the gene families from the first, "ABC-F ATP-binding cassette ribosomal protection protein", to the last "Rm3 family beta-lactamase" ). 


.. code-block:: r

  library(tidyverse)
  # Run pivot longer to prepare the data by including all the gene families to display in the plot.
  df.long <- amr_genes.final.norm %>% 
    pivot_longer("ABC-F ATP-binding cassette ribosomal protection protein":"Rm3 family beta-lactamase", names_to = 'variable', values_to = 'value')
  # we changed the order for displaying the charts
  df.long$group <- factor(df.long$group, levels = c("NTC", "Peterborough", "Ireland-Medieval", "UK-18-19th_c.", "Modern"))
  # generate the multipanel plots
  ggplot(data = df.long, aes(x = group, y = value, fill = group)) +
    geom_boxplot(lwd=0.1,
    outlier.size = 0.5,
    outlier.stroke = 0.1) +
    theme(axis.text.x = element_blank(),
    strip.text = element_text(size = 5),
    legend.title = element_text(colour = "black", size = 7, face = "bold"),
    legend.text = element_text(colour = "black", size = 6),
    axis.line = element_line(size=0.1),
    axis.ticks = element_line(size=0.1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=4),  	#color="#993333", face="bold", angle=45
    axis.title.y = element_text( size = 12, face = "bold" )) +
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(facets = ~variable, scales = 'free', 
    #nrow=8, 
    ncol=5
    )


