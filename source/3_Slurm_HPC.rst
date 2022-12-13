#########################
Working on HPC with Slurm
#########################

**************************
Connect to the HPC cluster
**************************

An High Performance Computing (HPC) cluster is a collection of many separate servers (computers), called **nodes**. 
HPC clusters are normally organized in:

  - login node, where users login when connecting remotely
  - compute nodes, where the majority of computations is run
  - "fat" compute nodes, which provide more RAM (1Tb or more)
  - GPU nodes, where computations can be run on both on CPU cores and on a Graphical Processing Unit

All cluster nodes have the same components as a laptop or desktop: CPU cores, memory and disk space. The difference between personal computer and a cluster node is in quantity, quality and power of the components.

The HPC cluster used here is called `Galileo100`_ 

  .. _Galileo100: https://www.hpc.cineca.it/hardware/galileo100
  
To login to Galileo the ``ssh`` command (Secure Shell), a cryptographic network protocol for operating network services securely over an unsecured network.
::

  ssh username@login.g100.cineca.it

You will be asked to insert your password, and after that you will be connected to the login node. 

Whenever you want to download a file from the HPC cluster in your local machine, you can do that with the command ``scp`` (Secure Copy), which allows files to be copied to, from, or between different hosts. It uses ssh for data transfer and provides the same authentication and same level of security as ssh. For example, to copy from a remote host (our server) to your computer:
::

  scp username@login.g100.cineca.it:/full_path_to_file /some/local/directory

To copy a folder you need to call the option ``-r``
::

  scp -r username@login.g100.cineca.it:/full_path_to_file /some/local/directory

.. note::

 If you are using a ``pem`` file to connect to the server, you have to use in order to download the files: 
 ::
 
   scp -i filename.pem -r username@remotehost:/full_path_to_file /some/local/directory


*******************
Run jobs with SLURM
*******************

In a HPC clusters you can either work interactively (not included here), or submit jobs via a workload manager/job scheduler. 
`SLURM`_ is a workload manager/job scheduler that is widely used to submit jobs in HPC clusters. 
To execute some commands (or also a script or a program) using SLURM you must prepare a Bash job-script that follows this scheme: 

  .. _SLURM: https://hpc-wiki.info/hpc/SLURM

..
  
  - The **shebang** line ``#!/bin/bash``
  - The SLURM directives identified with ``#SBATCH``
  - list of environment **modules** activated
  - list of environment **variables**
  - execution line (the full list of commands or the name of a script to execute)


Environments modules and variable
*********************************

Many commonly used programs have been pre-installed on the cluster and can be used by everyone through a system called **environment modules**. 
Software packages can be loaded in your HPC work environment through the command ``module``. Software packages are installed in a central, read-only location by the HPC support team and when you load a specific module a modulefile for the corresponding software is built and modifies your environment, so that the software becomes usable to you.
This system makes it easy for everyone to find and load pre-built software into their environment to use with their programs.

Modules are organized in functional categories (compilers, libraries, tools, applications…). A set of multiples modules (so libraries and tools) are collected in pre-organized ``profiles`` (they are created by the HPC support team). 
In Galileo we distinguish **programming profiles**, which contain only compilers, libraries and tools modules used for compilation, debugging, profiling and pre-processing activity:

  - profile/base: this is the default profile that you find once logged in. Basic compilers, libraries and tools (gnu, intel, intelmpi, openmpi--gnu, math libraries, profiling and debugging tools, rcm,…)
  - profile/advanced: Advanced compilers, libraries and tools (pgi, openmpi--intel, etc.)
  - profile/archive: this contains all old versions of applications and tools and the old versions of compilers and libraries not used by other modules

and **production profiles**, which contain only libraries, tools and applications modules used for production activity of specific fields of research. For example:

  - profile/chem: chemistry applications
  - profile/phys: physics applications
  - profile/lifesc: life science applications
  - profile/bioinf: bioinformatic applications
  
To load and unload a profile or a specific module: 
::
  
  module load profile/bioinf
  module unload profile/bioinf

To list the modules loaded in the environment: 
::
  
  module avail
  
To clean the environment from all the loaded modules:
::

  module purge

A series of general useful commands to explore modules available on Galileo: 
::
  
  modmap -categories	# to see all functional categories 
  modmap -c libraries	#Prints all modules available for the selected category (compilers, libraries, tools, applications)
  modmap -all		# to print the full modules map: all profiles, all categories defined for each profile, all modules available for each category
  modmap -profiles	# prints all profiles available
  modmap -p bioinf	# prints all modules available for a profile
  modmap -m samtools	# Prints all versions of the selected module and the profile/s where they are located. Usually, the default module is the latest version.

To show all the variables of a module and its prerequisites (compilers and libraries) and conflict: 
::

  module show samtools
  
Prerequisites of a module can be loaded ahead separately, otherwise you can use the command ``autoload`` to load the module and all the prerequisites: 
::

  module load autoload samtools

.. note::

  | **Environments variables**
  | Loading a module means that a series of useful **shell environment variables** will be set. An environment variable is a variable whose value is set outside of a program. A name (string) or a value (number) is assigned to the variable, and is available as reference for running a program or a script. To see what a module actually is doing you can use ``module show <module name>``. Modules typically load and append things to key environment variables like ``$PATH``, ``$LD_LIBRARY_PATH`` and so on. The variables can be used both in scripts and in the command line. They are usually referenced by putting the special symbols in front of or around the variable name. For instance, to display the user home directory, in most scripting environments, the user has to type:
  
  ::
  
      echo $HOME
      echo $CINECA_SCRATCH
    
  More in general, to create a variable in the environment where you are working: 
  ::
    
    export CHR_1="chr1.fa"	# to assign a value (a filename here) to the variable CHR_1 (the name of the variable is arbitrary)
    echo $CHR_1			# to visualize the value assigned
  
.. warning::

  Multiple profiles can be loaded at the same time. Keep in mind that to use a module placed under other profiles, you have to load the corresponding profile ahead. 
  ::
    
    module load autoload r	# this returns an error (+(0):ERROR:0: Unable to locate a modulefile for 'r') without having previously loaded profile/bioinf, where r is located
    module load profile/bioinf
    module load autoload r	# r is now loaded without errors
    
    
Preparation of Slurm scripts
****************************
As mentioned earlier, the Slurm script starts with a **shebang** header, followed by the `SLURM directives`_. Here is a standard configuration to run FastQC for example: 

  .. _SLURM directives: https://slurm.schedmd.com/sbatch.html

::

  #!/bin/bash
  # Slurm directives
  #SBATCH --job-name=fastqc			# Arbitrary name assigned to job to be executed
  #SBATCH -N 1					# (--nodes) number of computational nodes allocated to this job (min-max number may also be allocated)
  #SBATCH -n 1					# (--ntasks) this allocates the cpus. The controller will allocate one cpu per task. To use for multithreading operations
  #SBATCH --partition=gll_usr_prod		# the partition in which the job is executed
  #SBATCH --mem=10GB				# Specify the real memory required per node. Default units are megabytes. Different units can be specified using the suffix [K|M|G|T]. 
  #SBATCH --time=02:00:00			# minimum time limit on the job allocation. 
  #SBATCH --account=project_account		# name of the project account with computational hours allocated
  #SBATCH --mail-type=ALL			# to receive notifications about the job change of status
  #SBATCH --mail-user=youremail			# email address to receive notifications

Then we continue by activating the modules necessary to run the job: 

::
  
  # Module(s) loading
  module load profile/bioinf
  module load autoload fastqc
  
After that we can define specific variables, for example to define the samples that will be processed. Once the variable is assigned at this initial stage, we can use that all over the commands and the scripts that compose the slurm job (in this case a very easy fastqc command).
::

  # Setting variables
  SAMPLE="SRR5956372"
  
Finally we type the commands: 
::

  # Commands
  mkdir fastqc_output
  fastqc nogroup extract ${SAMPLE}_R1.fastq.gz -o fastqc_output
  fastqc nogroup extract ${SAMPLE}_R2.fastq.gz -o fastqc_output

