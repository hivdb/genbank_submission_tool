# HIVDB GenBank submission helper tools

This program can help you prepare sequences, meta data, and sequence annotations. It's designed for sanger sequencing results or consensus sequence of NGS sequencing.

## How to use the program

1. [Prepare the sequence and meta data files](#prepare-the-sequence-and-meta-data-files)
2. [Run the script](#run-the-script)
3. [Check BankIt files for submission](#check-bankit-files-for-submission)

### Prepare the sequence and meta data files

Two files are required:

1) sequence files in **fasta** format, file name is `unaligned_sequences.fasta`.
2) meta data file in **csv** format, file name is `source_modifier.csv`.


The meta data file should contain a column called `Isolate`, the `Isolate` is used as the unique ID for each sequence, it's also used in the fasta file to identify sequences.

The meta data should contain the columns below:

- `Isolate`
- `Collection_date` (formatted like '2024-01-01')
- `Country`
    - Please check the name is correct, https://www.ncbi.nlm.nih.gov/genbank/collab/country/
- `Isolation source`
- `Genes`

Additional columns can also be added to the meta data file

- `note`
    - this column can contain deidentified patient code, treatment history in raw text

The GenBank document [Preparing a Source Modifiers Table File for All Source Modifiers
](https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html) lists all the headers. If you want to add additional headers other than these headers listed above, **please create a new issue and let us know**.

Note: The program assumes the Organism is *Human immunodeficiency virus 1*, and the Host is *Homo sapiens*, which can be configured using script arguments.


### Run the script

Save `unaligned_sequences.fasta` and  `source_modifier.csv` under `input_files` folder.

```shell
make init
make
```

The program will do several things:

- check the uniqueness of `Isolate`
- align the sequence and get alignment information
- generate two `BankIt files` in `bankit_files` folder
    - `bankit.fasta`
    - `features.txt`

#### Exclude sequence

You can create a file `input_files/aligned_ignore.csv` to exclude some sequences with issues. The header of this file should contain **Isolate** and **gene** columns

### Check BankIt files for submission

Before submitting your sequences and meta data to GenBank, you should double check the generated two files.

- check sequence quality file `input_files/aligned_meta.csv`
- check `bankit.fasta`
    - check the positions of stop codons are correct.
    - check the positions of deletions and gaps are correct.
    - run the [HIVdb Program](https://hivdb.stanford.edu/hivdb/by-sequences/), compare the result of original sequence and the corresponding sequence in this file, the mutations should be the same.
- check `features.txt` file
    - check the start stop positions
    - check sequence with stop codons has misc feature
    - check sequence with gaps has gap information
    - check sequence with untranslatable codons has misc feature
    - check the gene name


## Create BankIt submission


### Prepare other information before submission

- Contact
- Reference
- Sequencing technology

Here is the [BankIt](https://www.ncbi.nlm.nih.gov/WebSub/) tool.

Please see below are some issues you would encounter and how to resolve them:

- In `Nucleotide` tab, you should upload the `bankit.fasta` file
- In `Features` tab, you should upload the `features.txt` file

- Warning: There are one or more significant strings of NNNs (length>10). Please explain what the strings of internal NNNs represent
    - choose *a region of estimated length between the sequenced regions based on an alignment to similar sequences or genome*
- Submission Set/Bach
    - choose *Pop set*
- Warning: Terminal ends of the following sequence(s) are low quanlity (too many ambiguous bases) and have been trimmed.
    - choose *or, click here to undo all trimming, and then click Continue to submit the original untrimmed sequences(s)*
    - note: the BankIt may remove one nucleotide at the end because of ambiguity, which cause the last codon can not be translated
    - after your choice, at the end of the page, it will show *bankit.fasta+(untrimmed+original)*
- Sequence(s) and Definition Line(s), Molecule Type
    - if the sequences were isolated from plasma, choose **genomic RNA**
- Tab `Sequencing Technology`
    - if the sequences were not using NGS methods, please don't choose any of th options `unassembled sequence reads`, `assembed sequences (consisting of two or more sequence reads)`.


## Advanced usage

You can use your prefered alignment tool to prepare the aligned sequence and save in file `input_files/aligned_meta.csv`

Several columns are required

- `Isolate`
- `gene`
    - gene name
- `insertions`
    - insertion mutation list, for example (S255N_K)
- `deletions`
    - deletion position list
- `stops`
    - stop mutation list for example (E170*)
- `gene_AA_length`
    - gene amino acid sequence length
- `gene_NA_length`
    - gene nucleotide sequence length
- `aligne_NA_length`
    - nucleotide sequence length after alignment
- `start_AA_pos`
- `start_NA_pos`
- `stop_AA_pos`
- `stop_NA_pos`
- `translatable`
- `untrans_reason`
    - if the sequence is not translatable, please provide the reason
- `aligned_NA`
    - the aligned_NA should be the same length as `gene_NA_length`, and heading or tailing unsequenced positions should use `.` to indicate tehem.

Then you can use the command line to generate `BankIt files`, please refer to [Makefile](./Makefile) for how to use the scripts.
