# HTS post-sequencing protocol for ‘vsearch’

### Steps for post-sequencing based on Arnold Lab protocol

### Protocol developed by Carolina Sarmiento and modified by Bolívar Aponte Rolón

### September 4, 2020

# Setting up the workspace

• Set PATH or wd where you have the .fastqc.gz files Example:

```{bash}
baponte@baponte-VirtualBox:/media/sf_My_Drive/VBL_users/Grad Students/Bolivar/Dissertation/Leaf_Traits_Panama/Data/Arnold_Lab/Post_Sequencing/MiSeqRun_Bolivar_PoorDemultiplexing_20200228/FastQ_Info/Mock_Communities)
```

• In command line execute: gzip -d \*.fastqc.gz This will decompress all gunzip files (i.e. sample files, usearch download files)

• This decompressing might take a while depending on the amount of files. After this is finished you can proceed to execute FastQC analysis on individual sequences or concatenated sequences.

• For individual sequence reads execute the following: fastqc \*.fastqc

• For concatenate files execute the following first:

```{bash}
cat *R1.fastq > FILENAME.fastq; then fastqc *.fastqc
```

• Once files have been decompressed and FastQC reports have been generated, proceed to execute USearch or Vsearch analyses for MaxEE and CutOffs.

\*Note: Usearch and Vsearch do not have a graphical interface. It only runs in the terminal/command prompt. The executable file has to be in the same folder as the samples to be analyzed are. It will not work otherwise.

# Installing USEARCH

• Download Usearch from: http://www.drive5.com/usearch/download.html

USEARCH is free only for the 32bit version. The 64bit version is not free and not accessible. Download the corresponding version for your computer operating system: Windows, Linux or Mac. • In the file that you downloaded or stored USEARCH proceed to decompress Usearch if you haven't already, execute: gzip -d usearch11.0.667_i86linux32

• So you do not have to call usearch by that long name you can create a symbolic link or change completely in the terminal or rename file in the user interface.

• Execute in Linux:

```{bash}

ln -s \[CURRENT USEARCH VERSION\] usearch
```

• Execute in MacOS:

```{bash}

ls -s \[CURRENT USEARCH VERSION\] usearch
```

• If this does not work you can rename the file all together by executing:

```{bash}
mv usearch11.0.667_i86linux32 usearch
```

Now you should be able to call usearch by just typing usearch. If error occurs due to permissions, run: chmod a+x usearch

This should give you all the authorizations required to execute commands during the session. In Linux, USEARCH will run and execute commands on single files with no problem. A typical command prompt looks like this:

```         
  <usearch -cluster_fast seqs.fasta -id 0.9 -[FILE NAME].fasta>
     |         |           |             |     |               |
     |       Command   Input filename    | Output option with filename
     Binary filename                     Parameter option(s)
```

• See other commands at: https://drive5.com/usearch/manual/cmds_all.htm

To automate the process, Carolina Sarmiento, wrote a simple “for loop” to run commands on multiple files.

```{bash}
       for x in $(ls [FILE PATH]/*)
         do
       vsearch -fastq_eestats2 $x -output $x.txt -length_cutoffs 100,300,10 -ee_cutoffs 0.25,0.5,1.0
         Done
```

Notice the use of 'vsearch' instead of 'usearch'. Theoretically, this loop should run with usearch, but several of us ran into error and trouble executing the command with 'usearch' in Linux/MacOS. The work around this is to install 'vsearch'. For it to work with ‘usearch’ make sure to use ‘usearch’ commands, see manual for reference.

# Installing Vsearch

It is an open-source version of 'usearch'. • “We have implemented a tool called VSEARCH which supports de novo and reference based chimera detection, clustering, full-length and prefix dereplication, rereplication, reverse complementation, masking, all-vs-all pairwise global alignment, exact and global alignment searching, shuffling, subsampling and sorting. It also supports FASTQ file analysis, filtering, conversion and merging of paired-end reads. VSEARCH stands for vectorized search, as the tool takes advantage of parallelism in the form of SIMD vectorization as well as multiple threads to perform accurate alignments at high speed. VSEARCH uses an optimal global aligner (full dynamic programming Needleman-Wunsch), in contrast to USEARCH which by default uses a heuristic seed and extended aligner. This usually results in more accurate alignments and overall improved sensitivity (recall) with VSEARCH, especially for alignments with gaps.” -https://github.com/torognes/vsearch

Download and install • Source distribution To download the source distribution from a release and build the executable and the documentation, use the following commands:

```{bash}
wget https://github.com/torognes/vsearch/archive/v2.15.0.tar.gz 
tar xzf v2.15.0.tar.gz 
cd vsearch-2.15.0 
./autogen.sh 
./configure 
make 
make install # as root or sudo make install
```

Cloning the repository • Cloning the repo instead of downloading the source distribution as a compressed archive, you could clone the repo and build it as shown below. The options to configure as described above are still valid.

```{bash}
git clone https://github.com/torognes/vsearch.git 
cd vsearch .
/autogen.sh 
./configure 
make 
make install # as root or sudo make install

```

• To install, command:

```{bash}
sudo apt install vsearch
```

More information and installing options can be found at: https://github.com/torognes/vsearch

# Filtering and trimming sequences with Vsearch

## Deciding which MaxEE and Cutoff to use?

### Filtering to chosen parameters

You can do these next 2 steps on a per single file basis or you can concatenate all the files into a single master file of all your samples and then proceed to step (B). If you decide to concatenate, here is the code:

```{bash}
	cat *.fastq > NEW_FILE_NAME.fastq
```

```         
A. Looping multiples sequence files of `.fastq` extension to `.fa` extension:
```

```{bash}
for file in $(ls File_path/*)
>do
>vsearch -fastq_filter $file -fastq_maxee 0.5 -fastq_trunclen 270 -fastaout $file.fa -relabel $file
>done
```

```         
B. For single file filtering:
```

```{bash}
>vsearch -fastq_filter File_Name.fastq -fastq_maxee 0.5 -fastq_trunclen 100 -fastaout New_File_Name.fa -relabel New_File
_Name
```

It is useful to include in the new file name the MaxEE and Cutoff parameters used. It can turn into a long name, but you will be certain of the parameters used on the sequence data. For example:

```{bash}
>New_File_Name_maxee05_cutoff.fa
```

# Concatenating all sample files and cleaning file for cluster analysis

If you have not concatenated all the sample files, go back to the previous section and see the corresponding step. If you have already done this then it is good to open the file with a text editor such as Notepad (for Windows) or similar text editor for MacOS. See some alternatives here. What you want to be looking for is the addition of the FILE_PATH and NEW_FILE_NAME to your actual sample identification. This can look like this:
```
> "G:\media\sf_My_Drive\VBL_users\Grad_Students\Bolivar\Dissertation\Leaf_T
raits_Panama\Data\Arnold_Lab\Post_Sequencing\MiSeqRun_Bolivar_PoorDemultipl
exing_20200228\FastQ_Info\Mock_Communities\reverse"
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
> "G:\media\sf_My_Drive\VBL_users\Grad_Students\Bolivar\Dissertation\Leaf_T
raits_Panama\Data\Arnold_Lab\Post_Sequencing\MiSeqRun_Bolivar_PoorDemultipl
exing_20200228\FastQ_Info\Mock_Communities\reverse"
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
```
What we want:
```
> SAMPLE_NAME
exing_20200228\FastQ_Info\Mock_Communities\reverse"
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
AAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAAAAGTCTGTAACTGAAAA
```


You can see the added file path that makes for a long sample name. The second example show the sample name we want. To do this we can open the file with a text editor like VIM for Linux operating systems and clean the sample names in an automated or faster fashion. This can be done with each sample file or the concatenated file of all samples. For searching and replacing terms in the file you can use VIM in Linux operating systems. You can install VIM by executing the following command:

```{bash}
sudo apt install vim
```

Some useful references on how to use ‘vim’: • https://vim.fandom.com/wiki/Search_and_replace • https://www.tutorialspoint.com/unix_commands/vim.htm • https://www.howtoforge.com/vim-basics

To eliminate and replace the complete file path in one command, follow this example:

```         
:%s/media\|sf_My_Drive\|VBL_users\|Grad_Students\|Bolivar\|Dissertation\|Leaf_T
raits_Panama\|Data\|Arnold_Lab\|Post_Sequencing\|MiSeqRun_Bolivar_PoorDemultipl
exing_20200228\|FastQ_Info\|Mock_Communities\|reverse//g
```

This line follows the basic syntax of :%s/search/replace/g where “search’ is the term we want to change with the term ‘replace’. For multiple terms it requires the use of backslash and the vertical bar (\|). To make sure that you are selecting all the desired terms you can activate the highlighted search mode on ‘vim’ by executing: :set hlsearch . This will allow you to see what you are selecting in real time before executing the command and later finding out any errors.

Extra commands On each line, delete all text following the whole word "foo" (to end of line).

```         
:%s/\R1\zs.*//g
```

# Cluster into OTUs

1.  To pull out unique sequences in ‘vsearch’ the command -fastx_uniques does not exist. The equivalent seems to be dereplication that can be used as `--derep_fulllength` or with the `--minseqlength` option. Also, the option `-fastaout` changes for `--output` and here we can remove singletons `minuniquesize 2`:

```{bash}
>vsearch -derep_fulllength FILE_PATH/FILE_NAME.fa --minuniquesize 2 -sizeout -relabel Uniq --output FILE_PATH/NEW_FILE_NAME.fa
```

a.  With this command you will want to access files in a specific folder and move to a new folder to keep them separated and organized. For example, moving from a folder named FILTERED to one named UNIQUE. Again, an additional way to keep track of the progress with the files is to add the step or command to the file name (e.g. NEW_FILE_NAME_unique_nonsingletons.fa).


2.  Use `Unoise3` to denoise and cluster unique sequences at 100% similarity, in `vsearch` commands are different. `Unoise3` is performed instead by `--cluster_unoise`. The caveat is that it does not remove chimeras!
    a.  So first, denoise the reads: Denoising attempts to identify all correct biological sequences in the reads and gives as output a Zotu = a denoised sequence (zero-radius otu):

```{bash}
>vsearch --cluster_unoise FILE_PATH/NEW_FILE_NAME_uniques_nosingletons_0.5_270.fa --centroids FILE_PATH/NEW_FILE_NAME_uniques_nosingletons_unoise_0.5_270.fa
```

3.  Next, we have to remove chimeras - the names of each otu are the 'old identifiers' and it includes the abundance (e.g. Uniq1;size=124513). We can change the names if we want with the option --relabel 'prefix' and Use --sizeout to conserve the abundance annotations:

```{bash}
>vsearch --uchime3_denovo FILE_PATH/FILE_NAME_unique_nosingletons_unoise_0.5_270.fa --relabel Zotu --sizeout --nonchimeras FILE_PATH/NEW_FILE_NAME_unique_nosingletons_unoise_nonchimeras_0.5_270.fa
```

4.  Cluster into 90, 95, 97, 99% similarity: a. Loop:

```{bash}
for id in 99 97 95 90
>do
>vsearch -cluster_smallmem FILE_PATH/NEW_FILE_NAME_unique_nosingletons_unoise_nonchimeras_0.5_270.fa -id 0.$id --sizein -relabel OTU_ --centroids FILE_PATH/NEW_FILE_NAMEzotu$id.fa
>done
```

\*With this command you will want to access files in a specific folder and move to a new folder to keep them separated and organized. For example, moving from a folder named UNIQUE to one named OTU_FILES.

b.  For a single %similarity:

```{bash}
    vsearch --cluster_smallmem FILE_PATH/NEW_FILE_NAME_unique_nosfoingletons_unoise_nonchimeras_0.5_270.fa --id 0.97 --sizeout --relabel OTU_ --centroids FILE_PATH/NEW_FILE_NAMEzotu_97.fa
```

5.  Make the final OTU tables:

a.  For a single %similarity:

```{bash}
vsearch --usearch_global FILE_PATH/FILTERED/FILE_NAME.fa --db FILE_PATH/OTU_FILES/FILE_NAME_zotu97.fa --strand plus --id 0.97 --otutabout FILE_PATH/OTU_TABLES/NEW_FILE_NAME_otuFINAL_97.txt
```

\*Notice that we have selected the original “filtered” file to our selected maxee and cutoff parameters. Then we select the OTU file of our preference (e.g. 97% similarity) and our output is an OTU table based on the selected similarity.

b.  Loop for all:

```{bash}
for id in 99 97 95 90
do
vsearch --usearch_global FILE_PATH/FILTERED/FILE_NAME.fa --db FILE_PATH/OTU_FILES/FILE_NAME_zotu$id.fa --strand plus --id 0.$id --otutabout FILE_PATH/OTU_TABLES/NEW_FILE_NAME_otuFINAL_$id.txt
done
```

Examples of successful loops for clustering

#August 13, 2020 \*Loop used:

```{bash}
for x in 99 97 95 90 \
do \
> vsearch -usearch_global /media/sf_My_Drive/VBL_users/Grad_Students/Bolivar/Dissertation/Leaf_Traits_Panama/Data/Post_Sequencing/MiSeqRun_Bolivar_PoorDemultiplexing_20200228/FastQ_Info/Mock_Communities/Filtered_July17/All_concat_R1_filtered_jul17.fa --db /media/sf_My_Drive/VBL_users/Grad_Students/Bolivar/Dissertation/Leaf_Traits_Panama/Data/Post_Sequencing/MiSeqRun_Bolivar_PoorDemultiplexing_20200228/FastQ_Info/Mock_Communities/OTU_files_July17/All_Mock_R1_jul17_zotu$x.fa --strand plus --id 0.$x --otutabout /media/sf_My_Drive/VBL_users/Grad_Students/Bolivar/Dissertation/Leaf_Traits_Panama/Data/Post_Sequencing/MiSeqRun_Bolivar_PoorDemultiplexing_20200228/FastQ_Info/Mock_Communities/OTU_tables_July17/All_Mock_R1_otuFINAL\_\$x.txt \
> done
```
