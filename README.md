This program reads sumstat files in various formats, and outputs a uniform format.

- install julia: 'conda install -c conda-forge julia; conda activate julia'

- Download the repository: 'git clone https://github.com/ymer/process_sumstats' 

- Install Julia packages: 'julia julia_install.jl'

- See how to run the program: 'julia ProcessSumstats.jl'

- Run with basic settings: 'julia ProcessSumstats.jl [inputfile] [outputfile]'

Columns in output:

- `SNP` - SNP ID in the format Chr:Pos

- `Rsid` - SNP ID as a RS-number

- `A1` - The effect allele

- `A2` - The other allele

- `P`

- `Beta`

- `Z`

- `n`



