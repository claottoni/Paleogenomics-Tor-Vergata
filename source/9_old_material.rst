##############################
Material from previous courses
##############################


.. _kraken_label:

******
Kraken
******

In this hands-on session we will use `Kraken`_ to screen the metagenomic content of a DNA extract after shotgun sequencing. 
A Kraken database is a directory containing at least 4 files:

  - **database.kdb**: Contains the k-mer to taxon mappings
  - **database.idx**: Contains minimizer offset locations in database.kdb
  - **taxonomy/nodes.dmp**: Taxonomy tree structure + ranks
  - **taxonomy/names.dmp**: Taxonomy names


Minikraken
**********

We will first use a pre-built 8 GB Kraken database, called `Minikraken`_, constructed from complete dusted bacterial, archaeal, and viral genomes in RefSeq (as of October 2017).
You can download the pre-built Minikraken database from the website with ``wget``, and extract the archive content with ``tar``: 

  .. _Kraken: http://ccb.jhu.edu/software/kraken/
  .. _Minikraken: https://ccb.jhu.edu/software/kraken/

::

  wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_8GB_dustmasked.tgz
  tar -xvzf minikraken_20171101_8GB_dustmasked.tgz

Then, we can run the taxonomic assignation of the reads in our sample with the ``kraken`` command
::

  kraken --db minikraken_20171101_8GB_dustmasked --fastq-input filename.gz --gzip-compressed --output filename.kraken



Some of the options available in Kraken:  

=======================  ========
Option                   Function
=======================  ========
**--db** <string>        Path to the folder (database name) containing the database files.  
**--output** <string>    Print output to filename.
**--threads** <integer>  Number of threads (only when multiple cores are used).
**--fasta-input**	     Input is FASTA format.
**--fastq-input**	     Input is FASTQ format.
**--gzip-compressed**    Input is gzip compressed.
=======================  ========


Create report files
******************* 

In Kraken 1, report files are generated with a specific command, after the classification (section 3.1.2: `Create report files`_). Once the taxonomic assignation is done, from the Kraken output file we create a report of the analysis by running the ``kraken-report`` script. Note that the database used must be the same as the one used to generate the output file in the command above. The output file is a tab-delimited file with the following fields, from left to right: 

  1. Percentage of reads covered by the clade rooted at this taxon
  2. Number of reads covered by the clade rooted at this taxon
  3. Number of reads assigned directly to this taxon
  4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
  5. NCBI taxonomy ID
  6. Indented scientific name

Notice that we will have to redirect the output to a file with ``>``.
::

  kraken-report --db Minikraken filename.kraken > filename.kraken.report
 
.. note:: We can use a ``for`` loop to make the taxonomic assignation and create the report file for multiple samples. Notice the assignation of variables ``filename`` and ``fname`` to return output files named after the sample. 
  ::

    for i in *.fastq
    do 
      filename=$(basename "$i")
      fname="${filename%.fastq}"
      kraken --db Minikraken --threads 4 --fastq-input $i --output /${fname}.kraken
      kraken-report --db Minikraken ${fname}.kraken > ${fname}.kraken.report
    done

To visualize the results of the classification in multi-layerd pie charts, use ``Krona``, as described in the section 3.1.2 (see :ref:`krona-label`).


Building a Kraken standard database (on HPC clusters) 
*****************************************************

The pre-built Minikraken database is useful for a quick metagenomic screening of shotgun data. However, by building larger databases (i.e. a larger set of k-mers gathered) we may increase the sensitivity of the analysis. 
One option is to build the Kraken standard database. To create this database we use the command ``kraken-build``, which downlads the ``RefSeq`` complete genomes for the bacterial, archaeal, and viral domains, and builds the database. 
::

  kraken-build --standard --db standardkraken.folder

.. note:: 
  - Usage of the database will require users to keep only the ``database.idx``, ``database.kdb``, ``taxonomy/nodes.dmp`` and ``taxonomy/names.dmp`` files. 
    During the database building process some intermediate file are created that may be removed afterwards with the command: 
    ::
    
      kraken-build --db standardkraken.folder --clean

  - The downloaded RefSeq genomes require 33GB of disk space. The build process will then require approximately 450GB of additional disk space. The final ``database.idx``, ``database.kdb``, and ``taxonomy/`` files require 200 Gb of disk space, and running one sample against such database requires 175 Gb of RAM. 


Building a Kraken custom database (on HPC clusters)
*************************************************** 

Building kraken custum databases is computationally intensive. You will find a ready to use database. 
Kraken also allows creation of customized databases, where we can choose which sequences to include and the final size of the database. For example if you do not have the computational resources to build and run analyses with a full database of bacterial genomes (or you don't need to), you may want to build a custom database with only the genomes needed for your application. 

1. First of all we choose a name for our database and we create a folder with that name using ``mkdir``. Let's call the database ``CustomDB``. This will be the name used in all the dollowing commands after the ``--db`` option. 

2. Download NCBI taxonomy files (the sequence ID to taxon map, the taxonomic names and tree information) with ``kraken-build --download-taxonomy``. 
   The taxonomy files are necessary to associate a taxon to the sequence identifier (the GI number in NCBI) of the ``fasta`` sequences composing our database. For this reason we will build our database only with sequences from the NCBI RefSeq. 
   For more information on the NCBI taxonomy visit click `here`_. This command will create a sub-folder ``taxonomy/`` inside our CustomDB folder:  
   
     .. _here: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi. 
   
   ::
   
     kraken-build --download-taxonomy --threads 4 --db CustomDB
   

3. Install a genomic library. RefSeq genomes in fasta file from five standard groups are made easily available in Kraken with the command **kraken-build --download-library**:

  - *bacteria* : RefSeq complete bacterial genomes
  - *archaea* : RefSeq complete archaeal genomes
  - *plasmid* : RefSeq plasmid sequences
  - *viral* : RefSeq complete viral genomes
  - *human* : GRCh38 human genome

  The following command will download all the RefSeq bacterial genomes (33Gb size) and create a folder ``library/`` with a sub-folder ``bacteria/`` inside your CustomDB folder:  
  ::

    kraken-build --download-library bacteria --threads 4 --db CustomDB

4. We can add any sort of RefSeq ``fasta`` sequences to the library with ``kraken-build --add-to-library``. For example we will add to the library of bacterial genomes the RefSeq sequences of mitochodrial genomes. The sequences will be inside the sub-folder ``added/``.
   ::
   
     kraken-build --add-to-library mitochondrion.1.1.genomic.fna --threads 4 --db CustomDB
     kraken-build --add-to-library mitochondrion.2.1.genomic.fna --threads 4 --db CustomDB

   .. note:: If you have several ``fasta`` files to add you can use a ``for`` loop: 
     ::  

       for i in *.fasta
       do
         kraken-build --add-to-library $i --threads 4 --db CustomDB
       done

         
5. When analyzing a metagenomics sample using a ``Kraken`` database the primary source of false positive hits is represented by low-complexity sequences in the genomes themselves (e.g., a string of 31 or more consecutive A's). 
   For this reason, once gathered all the genomes that we want to use for our custom database, low-complexity regions have to be 'dusted'. 
   The program ``dustmasker`` from `Blast+`_ identifies low-complexity regions and **soft-mask** them (the corresponding sequence is turned to lower-case letters). 
   With a ``for`` loop we run dustmasker on each ``fasta`` file present in the library folder, and we will pipe (``|``) to dustmasker a ``sed`` command to replace the low-complexity regions (lower-case) with Ns. 
   Notice that the output is redirected (``>``) to a temporary file, which is afterwards renamed to replace the original file ``fasta`` file with the command ``mv``.
   
     .. _Blast+: https://www.ncbi.nlm.nih.gov/books/NBK279681/
   
   ::
	
     for i in `find CustomDB/library \( -name '*.fna' -o -name '*.ffn' \)`
     do
       dustmasker -in $i -infmt fasta -outfmt fasta | sed -e '/>/!s/a\|c\|g\|t/N/g' > tempfile
       mv -f tempfile $i
     done

 
  Some of the options available in Dustmasker: 

  ===================== ========
  Option                Function
  ===================== ========
  **-in** <string>      input file name  
  **-infmt** <string>   input format (e.g. fasta)  
  **-outfmt8** <string> output format (fasta)
  ===================== ========
      
 
6. Finally, we build the database with ``kraken-build``. With this command, Kraken uses all the masked genomes contained in the library (bacteria and mtDNA RefSeq) to create a database of 31 bp-long k-mers. 
   We can choose the size of our custom database (hence the number of k-mers included, and the sensitivity) with the  ``--max-db-size`` option (8 Gb here). 
   ::
   
     kraken-build --build --max-db-size 8 --db CustomDB


Taxonomic assignation with Kraken custom database
************************************************* 

Once our custom database is built we can run the command for taxonomic assignation of DNA reads agaisnt the custom database, as in section 1.1 and 1.2. 
::

  kraken --db CustomDB --fastq-input merged.fastq.gz --gzip-compressed --output sample.kraken
  kraken-report --db CustomDB sample.kraken > sample.kraken.report


Or, again, we can loop the commands if we have various samples. 
::

  for i in *.fastq
  do 
    filename=$(basename "$i")
    fname="${filename%.fastq}"
    kraken --db CustomDB --threads 4 --fastq-input $i --output ${fname}.kraken
    kraken-report --db CustomDB ${fname}.kraken > ${fname}.kraken.report
  done


To visualize the results of the classification in multi-layerd pie charts, use ``Krona``, as described in the section 3.1.3: :ref:`krona-label`

