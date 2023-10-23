This program reads sumstat files in various formats, and outputs a uniform format.

- install julia, such as with conda install -c conda-forge julia; conda activate julia

- julia julia_install.jl # To install the necessary packages

- julia ProcessSumstats.jl --help # See how to run the program

Columns in output:

- `SNP` - SNP ID in the format Chr:Pos

- `Rsid` - SNP ID as a RS-number

- `A1` - The effect allele

- `A2` - The other allele

- `P`

- `Beta`

- `Z`

- `n`



