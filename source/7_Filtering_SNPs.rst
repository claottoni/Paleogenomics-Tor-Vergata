########################################
Filtering, annotating and combining SNPs
########################################

To investigate the genetic variants in the ``vcf`` files we will use the program `snpToolkit`_. 
The ``-h`` option will display the following message:

  .. _snpToolkit: https://github.com/Amine-Namouchi/snpToolkit

::

  $ snptoolkit -h
  usage: snptoolkit [-h] {explore,annotate,combine,viz,expand} ...

      snpToolkit takes vcf files, as well as bam files (optional) as inputs. The vcf files could be generated using samtools/bcftools, gatk HaplotypeCaller or freeBayes.
      Please visit https://github.com/Amine-Namouchi/snpToolkit for more information.

  positional arguments:
    {explore,annotate,combine,viz,expand}
                          commands
      explore             explore your vcf files before annotation
      annotate            Annotate one or multiple vcf files
      combine             combine snpToolkit output files in one alignment in fasta format
      viz                 visualize snptoolkit output files
      expand              expand existent list of polymorphic sites when new SNP output files are availble

  optional arguments:
    -h, --help            show this help message and exit


Five options are possible: ``explore``, ``annotate``, ``combine``, ``viz``, ``expand``. 

***************************
SNPs filtering and annotion
***************************

The ``snpToolkit annotate`` command will display general information about the usage of the program:

.. code-block:: bash

  snptoolkit annotate -h
  usage: snptoolkit annotate [-h] -i IDENTIFIER -g GENBANK [-p PROCESSORS] [-f EXCLUDECLOSESNPS] [-q QUALITY] [-d DEPTH] [-r RATIO] [-e EXCLUDE]

  optional arguments:
    -h, --help           show this help message and exit

  snpToolkit annotate required options:
    -i IDENTIFIER        provide a specific identifier to recognize the file(s) to be analyzed
    -g GENBANK           Pleae provide a genbank file

  snpToolkit annotate additional options:
    -p PROCESSORS        number of vcf files to be annotated in parallel default value [1]
    -f EXCLUDECLOSESNPS  exclude SNPs if the distance between them is lower then the specified window size in bp
    -q QUALITY           quality score to consider as a cutoff for variant calling. default value [20]
    -d DEPTH             minimum depth caverage. default value [3]
    -r RATIO             minimum ratio that correspond to the number of reads that has the mutated allele / total depth in that particular position. default
                        value [0]
    -e EXCLUDE           provide a tab file with genomic regions to exclude in this format: region start stop. region must correspond to the same name(s) of
                        chromsome and plasmids as in the genbank file

Here is a simple example on how to use snpToolkit:
::

  snptoolkit annotate -i vcf -g GCF_000009065.1_ASM906v1_genomic.gbff -d 3 -q 20 -r 0.9

snpToolkit can automatically recogninze ``vcf`` files generated with the following programs: ``samtools mpileup``, ``gatk HaplotyCaller`` and ``freeBayes``. The ``vcf`` files could be gzipped or not. In the command line above, snpToolkit will filter and annotate all SNPs in the ``vcf`` file(s) that fullfil the following criteria: ``quality >= 30``, ``depth of coverage >= 5`` and ``ratio >= 0.9``.

.. note:: For each SNP position, the ratio (r) is calculated as follows:

    r= dm / (dr + dm)

    - dr= Number of reads having the reference allele
    - dm= Number of reads having the mutated allele

The output file(s) of snpToolkit is a tabulated file(s) that you can open with Microsoft Excel and it will look as follow:


.. code-block:: bash

   ##snpToolkit=__version__
   ##commandline= snptoolkit annotate -i vcf -g GCF_000009065.1_ASM906v1_genomic.gbff -d 5 -q 30 -r 0.9 -p 4
   ##VcfFile=sample5.vcf.gz
   ##Total number of SNPs before snpToolkit processing: 406
   ##The options -f and -e were not used
   ##Filtred SNPs. Among the 406 SNPs, the number of those with a quality score >= 30, a depth >= 5 and a ratio >= 0.9 is: 218
   ##After mapping, SNPs were located in:
   ##NC_003131.1: Yersinia pestis CO92 plasmid pCD1, complete sequence 70305 bp
   ##NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   ##The mapped and annotated SNPs are distributed as follow:
   ##Location      Genes   RBS     tRNA    rRNA    ncRNA   Pseudogenes     intergenic      Synonymous      NonSynonumous
   ##SNPs in NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp 155     0       0       1       0       0       57      54      101
   ##SNPs in NC_003131.1: Yersinia pestis CO92 plasmid pCD1, complete sequence 70305 bp    2       0       0       0       0       0       3       1       1
   ##Syn=Synonymous NS=Non-Synonymous
   ##Coordinates   REF     SNP     Depth   Nb of reads REF Nb reads SNPs   Ratio   Quality Annotation      Product Orientation     Coordinates in gene     Ref codon       SNP codon       Ref AA  SNP AA  Coordinates protein     Effect  Location
   82      C       A       36      0       34      1.0     138.0   intergenic      .       +       .       -       -       -       -       -       -       NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   130     G       C       28      0       27      1.0     144.0   intergenic      .       +       .       -       -       -       -       -       -       NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   855     G       A       69      0       62      1.0     228.0   YPO_RS01010|asnC        transcriptional regulator AsnC  -       411     ACC     AC[T]   T       T       137     Syn     NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp


The header of the generated snpToolkit output file includes useful information e.g. raw number of SNPs, Number of filtered SNPs, SNPs distribution, etc... 
The SNPs annotation is organized in tab delimited table. The columns of this table are:

=======================  ========
Column name              Description
=======================  ========
Coordinates              SNP coordinate 
REF                      Reference allele
SNP                      New allele in analyzed sample 
Depth                    Total depth of coverage 
Nb of reads REF          Number of reads with the reference allele
Nb reads SNPs            Number of reads with the new allele
Ratio                    Nb reads SNPs/(Nb of reads REF+Nb reads SNPs)
Quality                  Quality score
Annotation               Distribution within genes or intergenic
Product                  Functional product of the gene
Orientation              Gene orientation
Coordinates in gene      Coordinate of the SNP within the gene
Ref codon                Reference codon, ACC in the example above
SNP codon                New codon, AC[T]
Ref AA                   Amino Acid corresponding to reference codon 
SNP AA                   Amino Acid corresponding to new codon
Coordinates protein      Coordinate of the Amino Acid 
Effect                   Could be Synonymous (Syn) or Non-Synonymous (NS)
Location                 ID of the chromosome and plasmids.
=======================  ========


*********************************************
Compare and combine multiple annotation files
*********************************************

After generating a set of output files, you can run ``snpToolkit combine``:
::

  $ snptoolkit combine  -h
  usage: snptoolkit combine [-h] --location LOCATION [-r RATIO] [--bam BAMFILTER BAMFILTER BAMFILTER] [--snps {ns,s,all,inter}] [-e EXCLUDE]
  
  optional arguments:
    -h, --help            show this help message and exit
  
  snpToolkit combine required options:
    --location LOCATION   provide for example the name of the chromosome or plasmid you want to create fasta alignemnt for
  
  snpToolkit additional options:
    -r RATIO              new versus reference allele ratio to filter SNPs from snpToolkit outputs. default [0]
    --bam BAMFILTER BAMFILTER BAMFILTER
                          provide the depth, ratio and the path to the folder containing the bam files. eg. 3 0.9 path
    --snps {ns,s,all,inter}
                          Specify if you want to concatenate all SNPs or just synonymous (s), non-synonymous (ns) or intergenic (inter) SNPs. default [all]
    -e EXCLUDE            Provide a yaml file with keywords and coordinates to be excluded                          

``snpToolkit combine`` will compare all the SNPs identified in each file and create two additional output files: 

  1) a tabulated files with all polymorphic sites
  2) a ``fasta`` file. 


To combine the snps from different samples in one alignment ``fasta`` file you type the following command:
::

  snptoolkit combine --loc NC_003143.1 -r 0.9 --bam 2 1.0 ../bam

As we will be working with ancient DNA, a small fraction of your genome could be covered. In this case we will use the option ``--bam`` to indicate the path to the folder containing the ``bam`` files. 
The option ``-d`` must be used with the option ``--bam``. By default, all SNPs will be reported. This behaviour can be changed using the option ``--snp``.

.. note :: It is also possible to use the option ``--bam`` with modern data as some genomic regions could be deleted. 

The file reporting the polymorphic sites is organized as follows:

==== =========== === === ============================ ======= ======= ======= =======
ID   Coordinates REF SNP Columns with SNP information sample1 sample2 sample3 sample4
==== =========== === === ============================ ======= ======= ======= =======
snp1 130         A   T                                1       1       1       1
snp2 855         C   G                                0       0       ?       1
snp3 1315        A   C                                1       1       0       0
snp4 12086       G   A                                1       0       ?       0
==== =========== === === ============================ ======= ======= ======= =======
 
The table above reports the distribution of all polymorphic sites in all provided files. 
As we provided the ``bam`` files of the ancient DNA samples, snpToolkit will check if the polymorphic sites (snp2 and snp4) are absent in sample3 
because there is no SNP in that positions or because the region where the snps are located is not covered. In the latter case, snpToolkit will add a question mark ``?`` that reflects a missing data. 
From the table above, it will be possible to generate a ``fasta`` file, like the one below:
::

  >Reference
  ATCGGGTATGCCAATGCGT
  >Sample1
  ACCGGGTATGCCAATGTGT
  >Sample2
  ATTGGGTATGCCAGTGCGT
  >Sample3
  ?TTGAGT?TGTCA?TACGT
  >Sample4
  ATCGGGTATGCCAATGCGT


The ``fasta`` output file will be used to generate a maximum likelihood tree using ``IQ-TREE``


********************************
Phylogenetic tree reconstruction
********************************

There are several tools to build phylogenetic trees. All of these tools, use an alignment file as input file. Now that we have generated an alignment file in ``fasta`` format, we will use ``IQ-TREE`` to build a maximum likelihood tree. 
We use ``IQ-TREE`` for several reasons:

- It performs a composition chi-square test for every sequence in the alignment. A sequence is denoted failed if its character composition significantly deviates from the average composition of the alignment.

- Availability of a wide variety of phylogenetic models. ``IQ-TREE`` uses `ModelFinder`_ to find the best substitution model that will be used directly to build the maximum likelihood phylogenetic tree.

- Multithreading 

To generate the phylogenetic tree type the following command using your ``fasta`` as input:
::

  iqtree -m MFP+ASC -s SNPs_alignment.fasta

==================== ============
IQ-TREE options      Function
==================== ============
**-m**               `Substitution`_ model name to use during the analysis.
**-s**               Alignment file 
==================== ============



More information of IQ-TREE can be found in the program's `tutorial`_

    .. _ModelFinder: https://www.ncbi.nlm.nih.gov/pubmed/28481363
    .. _tutorial: http://www.iqtree.org/doc/Tutorial
    .. _Substitution: http://www.iqtree.org/doc/Substitution-Models

The phylogenetic tree generated can be visualized using ``Figtree``, `download`_ it in your local machine and load the ``treefile`` output from IQ-TREE to visualize the tree.

    .. _download: http://tree.bio.ed.ac.uk/software/figtree