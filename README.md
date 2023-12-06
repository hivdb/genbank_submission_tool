# HIVDB Genbank submission helper tools

## How to use the program

### First: prepare files

You should prepare your sequence files in fasta format.
You should prepare file meta data in csv format.

The meta data file should contain a column called `Isolate`, the value of the column should be unique and is used as the **ID** of each sequence in the fasta file.

### Second: save the files

Save the fasta file in `input_files` folder and name it as `unaligned_sequences.fasta`

Save the sequence meta data in `input_files` folder and name it as `source_modifier.csv`. The headers should contain

- `Isolate`
- `Collection_date` (formatted like '2024-01-01')
- `Genes`
- `Country`
    - Please check the name is correct, https://www.ncbi.nlm.nih.gov/genbank/collab/country/
- `Isolation source`

You can provide other informations. If you want additional information to be added into the GenBank submission, please create a new issue and let us know.

### Third: run the script

```shell
make
```

The program will do several things below:

- check the consistency of meta data
- align the sequence and get alignment information
- generate `bankit.fasta` in `bankit_files` folder
- generate `features.txt` in `bankit_files` folder

### Fourth: Check result

- check the aligned sequences in the `bankit.fasta` are correct
    - check the sequences with stop codons
    - check the sequences with deletions
    - check the sequences with gaps
    - compare the HIVDB report with original sequence
- check the `features.txt` file
    - check the non function feature exist for some sequences with issues
    - check the gap annotation
    - check the gene name
    - check the start stop positions

### Final: submit the files

Please submit the `bankit.fasta` and `features.txt` to [BankIt](https://www.ncbi.nlm.nih.gov/WebSub/)
