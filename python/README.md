# pyabPOA: abPOA Python interface
## Introduction
pyabPOA provides a easy-to-use interface to [abPOA](https://github.com/yangao07/abPOA).

## Installation

### Install pyabPOA with conda or pip

pyabPOA can be installed with conda or pip:

```
conda install -c bioconda pyabpoa

pip install pyabPOA
```

### Install pyabPOA from source
Alternatively, you can install pyabPOA from source:
```
git clone https://github.com/yangao07/abPOA.git
cd abPOA/python
make install 
```

### Getting started
After installation, you can run the toy example script to test it:
```
python ./example.py
```

## Usage
```
import pyabPOA as pa
a = pa.msa_aligner()
seqs=[
'CCGAAGA',
'CCGAACTCGA',
'CCCGGAAGA',
'CCGAAGA',
]
res=a.msa(seqs, out_cons=True, out_msa=True, out_pog='pog.png') # perform multiple sequence alignment 
                                                                # generate a figure of alignment graph to pog.png

for seq in res.cons_seq:
    print(seq)  # print consensus sequence

res.print_msa() # print row-column multiple sequence alignment in PIR format
```

## APIs

### Class pyabPOA.msa_aligner
```
pyabPOA.msa_aligner(aln_mode='g', ...)
```
This constructs a multiple sequence alignment handler of pyabPOA, it accepts the following arguments:
* **aln_mode**: alignment mode. 'g': global, 'l': local, 'e': extension. Default: 'g'
* **match**: match score. Default: 2
* **gap_open1**: first gap open penalty. Default: 4
* **gap_ext1**: first gap extension penalty. Default: 2
* **gap_open2**: second gap open penalty. Default: 24
* **gap_ext2**: second gap extension penalty. Default: 1
* **extra_b**: first part of extra band width. Default: 10
* **extra_f**: second part of extra band width. Total extra band width: b+f\*L, L is the sequence length. Default : 0.01
* **is_diploid**: set as 1 if input is diploid data. Default: 0
* **min_freq**: minimum frequency of each consensus to output for diploid data. Default: 0.2

```
pyabPOA.msa_aligner.msa(seqs, out_cons, out_msa, out_pog=None)
```
This method performs mutliple sequence alignment and generates
* consensus sequence if `out_cons` is set as `True`
* row-column multiple sequence alignment in PIR format if `out_msa` is set as `True`
* figure of alignment graph if `out_pog` is set as a file name with suffix as `.png` or `.pdf`

### Class pyabPOA.msa_result
```
pyabPOA.msa_result(seq_n, cons_n, cons_len, ...)
```
This class is the output of `pyabPOA.msa_aligner.msa()`, it contains the information of generated consensus sequence and multiple sequence alignment. An object of this class has the following properties:
* **seq_n**: number of input sequences
* **cons_n**: number of generated consensus sequences
* **cons_len**: an array of consensus sequence length
* **cons_seq**: an array of consensus sequence
* **msa_len**: length of each row in multiple sequence alignment
* **msa_seq**: an array of `seq_n` rows of the multiple sequence alignment

pyabPOA.msa_result() also has a function of `print_msa`. It prints the row-column multiple sequence alignment.