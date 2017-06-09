# invariant-sites
Counts invariant/variant and masks variant sites from multi-FASTA alignment

## Author

Jason Kwong (@kwongjc)

## Dependencies
* Python 3.x
* BioPython
* Pandas

## Usage

```
$ invariant-sites.py -h
usage: 
  invariant-sites.py [--out invariant.fa] FASTA

Counts invariant/variant and masks non-invariant sites from multi-FASTA alignment

positional arguments:
  FASTA       original multi-FASTA alignment file

optional arguments:
  -h, --help  show this help message and exit
  --out FILE  specify output file with non-invariant sites masked
  --mask X    symbol to use for masking (default=X)
  --version   show program's version number and exit
```

***Warning: this tool uses a considerable amount of memory. Ensure you have capacity before running on large alignments.***

## Bugs

Please submit via the [GitHub issues page](https://github.com/kwongj/invariant-sites/issues).  

## Software Licence

[GPLv3](https://github.com/kwongj/invariant-sites/blob/master/LICENSE)
