# invariant-sites
Counts invariant/variant and masks variant sites from multi-FASTA alignment

## Author

Jason Kwong (@kwongjc)

## Dependencies
* Python 2.7.x
* BioPython

## Usage

```
$ invariant-sites.py -h
usage: 
  invariant-sites.py [--out invariant.fa] FASTA

Counts invariant/variant and masks variant sites from multi-FASTA alignment

positional arguments:
  FASTA       original multi-FASTA alignment file

optional arguments:
  -h, --help  show this help message and exit
  --out FILE  specify output file with variant sites masked (default = stdout)
  --mask N    symbol to use for masking; change this if the alignment already contains this symbol
  --version   show program's version number and exit
```

## Bugs

Please submit via the [GitHub issues page](https://github.com/kwongj/invariant-sites/issues).  

## Software Licence

[GPLv3](https://github.com/kwongj/invariant-sites/blob/master/LICENSE)

