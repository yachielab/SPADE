SPADE Installation and User Manual
===================================
SPADE is a software to explore various periodic repeat regions comprehensively from large genomic and protein data resources. The software first automatically extracts multiple sequence entries from an input file (GenBank or FASTA format) and identifies sequence type (DNA or protein) for each entry. Each sequence entry is scanned by a sliding window to count *k*-mers and highly repetitive regions are extracted. The sequence periodicity of each highly repetitive region is then evaluated based on position-period matrix that cumulatively plots distance between neighboring same *k*-mers and their sequence positions. The periodic sequence region is defined and the periodic sequence units are queried for a multiple alignment to identify repetitive motif and its sequence logo. The representative motif sequence is aligned back to the sequence of the periodically repeating region to annotate the repeating units. Finally, the annotations for detected periodic repeats are added to the input information and output in GenBank format with an option of visualizing *k*-mer density, position-periodicity matrix, sequence motif logo and repetitive unit loci with neighboring genes for each periodic repeat (Fig. 1a).

<img src=example/NC_008532.1/nucl_649122_650215/periodic_repeat.png width=300x300>
<img src=example/NC_008532.1/NC_008532.1_SPADE.gb.png width=350x350>

**Figure 1. Example visualizations of periodic repeat regions captured by SPADE. (a) A CRISPR-surrounding region of Streptococcus thermophilus LMD-9. (b) A circular genome map visualization of S. thermophilus LMD-9 genome with all of the periodic biomolecular sequences detected by SPADE.**

Software Dependencies
-------------------------------
SPADE works under Python 2.7.13 or Python 3.6.1 and require BLAST+ (ver 2.6.0 later) and MAFFT (ver 7.221 later) to be installed.


Installation
-------------------------------
1\. Obtain SPADE github packages using the following single command.

```git clone https://github.com/yachielab/SPADE```

To execute SPADE using linux/unix commands, add the SPADE directory to \$PATH and add executable authority to the python scripts in the directory.

````
cd SPADE
chmod u+x *.py
````

2\. Install the necessary Python packages using the following commands. 

````
pip install matplotlib==2.2.3
pip install seaborn==0.8.1
pip install weblogo==3.6.0
pip install biopython
````

If pip is not installed, please download get-pip.py](https://bootstrap.pypa.io/get-pip.py) and install pip using the following command.

```python get-pip.py```

\*If you use default Python in OSX or macOS and encounter the following message when get-pip.py is executed,

````
matplotlib 1.3.1 requires nose, which is not installed.
matplotlib 1.3.1 requires tornado, which is not installed.
````

please use the following commands,

````
pyhton get-pip.py --user
pip install --user -I seaborn
pip install --user -I weblogo
pip install --user -I biopython
````

and set PYTHONPATH to /Users/[USER_NAME]/Library/Python/2.7/lib/python/site-packages

3\. Install MAFFT if it is not installed. MAFFT package is available at the following links.

- [OSX or macOS](https://mafft.cbrc.jp/alignment/software/macosx.html)
- [Linux based distributions](https://mafft.cbrc.jp/alignment/software/linux.html)

Set \$PATH to the MAFFT executable.

4\. Install BLAST+ if it is not installed. Executable BLAST+ package is available at the following links.

- [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


Set \$PATH to the BLAST+ executable.

Example code
-------------------------------
The package contains a GenBank file for the *Streptococcus thermophilus*
LMD-9 genome so users can quickly test the functions of SPADE with the
following operations.

1\. Go to the example folder and execute SPADE to screen periodic repeats in the sequence entries included in the GenBank file.

```SPADE.py -in GCF_000014485.1_ASM1448v1_genomic.gbff```

As this GenBank file contains three sequence entries, SPADE creates three folders NC\_008532.1 for the genomic DNA, NC\_008500.1 and NC\_008501.1 for two plasmid DNAs. See below for the detail of output data format.

2\. If you want to additionally create a circular genome map with cumulative *k*-mer scores and detected periodic repeats as represented in Fig. 1b, go to the folder NC\_008532.1 and simply type the following command.

```genome_circular_plot.py NC_008532.1_SPADE.gb NC_008532.1_kmer_count.txt```

3\. If you want to additionally create a table representing all of the periodic repeats detected by SPADE, type the following command.

```get_spade_annotation.py NC_008532.1_SPADE.gb```

Output data format
-------------------------------
**Structures of output directories**

If your input a sequence file contains multiple sequence entries, SPADE creates multiple folders for different entries under your current directory. For GenBank files, each LOCUS entry is treated as one sequence entry and each folder can contain results for multiple DNA and protein sequences. In each \[entry\] folder, SPADE creates \[entry\].gbk that involves annotations for the detected periodic repeats in addition to the information in the original query file. In the same folder, a descendant folder is created to store all the results for each detected periodic repeat. The directory structure therefore looks as follows.

````
(current directory)/
     +----[entry].gbk
     +----[entry]/
          +----[type]_[start]_[end]/ 
````

where \[type\] is nucl or prot indicating that if the repeat sequence type is nucleotide or protein, respectively, and \[start\] and \[end] denote position of the repeat region in entry sequence.

**Output files in each descendant folder**

For each detected repeat, the following files will be created.

- **repeat.gbk** contains locus information for the detected repeat and its surrounding genes
- **k-mer.tsv** contains a cumulative *k*-mer count distribution across the detected repeat region
- **ppm.tsv** contains a position-period matrix for the detected repeat region
- **ppm4vis.tsv** contains a position-period matrix trimmed for visualization
- **pdist.tsv contains** a distribution of detected periods
- **unit\_seq.fasta** contains multiple sequences that are used to find a repeat motif
- **align.unit\_seq.fasta** contains a multiple alignment result of the detected repeat motif
- **periodic\_repeat.pdf** visualizes the contents in repeat.gbk, k-mer.tsv,
- **ppm4vis.tsv** and pdist.tsv as represented in Fig. 1a
- **weblogo.pdf** represents sequence logo for the detected motif
- **weblogo.txt** is a raw data used produce weblogo.pdf

Software usage
-------------------------------
````help
SYNOPSIS
  SPADE [-h] [--help] [-in input_file] [-f input_file_format] [-t sequence_type] 
    [-Nk kmer_size] [-Nw window_size] [-Ns kmer_score_threshold] [-Ng gap_size] 
    [-Nm region_margin] [-Np period_threshold] [-Nq gap_frequency_threshold] 
    [-Nu motif_letter_consistency] [-Nr non_consensus_length_threshold]
    [-Pk kmer_size] [-Pw window_size] [-Ps kmer_score_threshold] [-Pg gap_size] 
    [-Pm region_margin] [-Pp period_threshold] [-Pq gap_frequency_threshold]
    [-Pu motif_letter_consistency] [-Pr non_consensus_length_threshold]
    [--mafft string] [--blast string] [-n num_threads] [-v string] [-d] [--delete] 
    [-V] [--version]

DESCRIPTION
  SPADE 1.0.0

OPTIONAL ARGUMENTS
 -h, --help
   Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
 -V, --version
   Print software version; ignore all other parameters

 *** Input query options
 -in <File_In>
   Input file name
 -f <String, Permissible values: ‘genbank’ ‘fasta’ ‘auto’>
   Input file type
   Default = ‘auto’
 -t <String, Permissible values: ‘nucl’ ‘prot’ ‘auto’>
   Sequence type, nucleotide (nucl) or protein (prot)
   Default = ‘auto’

 *** General screening options for nucleotide periodic repeats
 -Nk <Integer>
   k-mer size
   Default = 10
 -Nw <Integer>
   Size of sliding window to calculate cumulative k-mer distribution
   Default = 1000
 -Ns <Integer>
   Threshold for peak height of each cumulative k-mer count area
   Default = 20 
 -Ng <Integer>
   Threshold for gap size between significant k-mer count areas
   Default = 200
 -Nm <Integer>
   Size of margin to be evaluated with each detected highly repetitive region
   Default = 1000
 -Np <Real>
   Periodicity score threshold for each detected highly repetitive region
   Default = 0.5
 -Nq <Real>
   Gap frequency threshold for each position of a repeat motif to be removed
   Default = 0.5
 -Nu <Real>
   Threshold for letter consistency score at each position of a repeat motif 
   Default = 0.8
 -Nr <Integer>
   Threshold for length of non-consensus region to be removed from a repeat motif
   Default = 5

 *** General screening options for protein periodic repeats
 -Pk <Integer>
   k-mer size
   Default = 3
 -Pw <Integer>
   Size of sliding window to calculate cumulative k-mer distribution
   Default = 300
 -Ps <Integer>
   Threshold for peak height of each cumulative k-mer count area
   Default = 6
 -Pg <Integer>
   Threshold for gap size between significant k-mer count areas
   Default = 50
 -Pm <Integer>
   Size of margin to be evaluated with each detected highly repetitive region
   Default = 300
 -Pp <Real>
   Periodicity score threshold for each detected highly repetitive region
   Default = 0.3
 -Pq <Real>
   Gap frequency threshold for each position of a repeat motif to be removed
   Default = 0.5
 -Pu <Real>
   Threshold for letter consistency score at each position of a repeat motif 
   Default = 0.8
 -Pr <Integer>
   Threshold for length of non-consensus region to be removed from a repeat motif
   Default = 5

 *** MAFFT and BLAST+ options 
 --mafft <'String'>
   Optional arguments for MAFFT can be defined with single quotations 
   Default = '--auto'
   For MAFFT optional arguments, see 
   https://mafft.cbrc.jp/alignment/software/manual/manual.html
 --blastn <'String'>
   Optional arguments for BLAST+ can be defined with single quotations
   Default = '-strand plus -task blastn-short -penalty -2 –outfmt "6 qseqid qseq sseqid sseq pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore"'
 --blastp <'String'>
   Optional arguments for BLAST+ can be defined with single quotations
   Default = '-task blastp-short –outfmt "6 qseqid qseq sseqid sseq pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore"'
   For BLAST+ optional arguments, see
   https://www.ncbi.nlm.nih.gov/books/NBK279684/

 *** Other options
 -n <Integer>
   Number of CPU threads. If this is set to more than 1, SPADE runs multiple 
   processes for multiple sequence entries in parallel.
   Default = 1
 -v <String, Permissible values: 'Y' 'N'>
   Generate pdf files to visualize results for each detected repeat region
   Default = Y
 -d, --delete
   This option deletes descendant output folders of highly repetitive regions
   that are detected not to contain periodic repeats
````

