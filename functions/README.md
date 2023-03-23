# Functional database

This code builds an SQLite database. 

At the moment we have provided default paths for some of the files, and they should be converted to online paths (e.g. ROADS?)

## Step 1. Building the database

You can build the database with:

```
python ~/GitHubs/CFHackathon2023/functions/create_sqlite.py -d database -o functions.sqlite -a -v
```

I usually run this as a slurm job like so:

```
#!/bin/bash
#SBATCH --job-name=SQLite
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o sqlite-%j.out
#SBATCH -e sqlite-%j.err

python ~/GitHubs/CFHackathon2023/functions/create_sqlite.py -d database -o functions.sqlite -a -v
```


## Step 2. Querying the database

We have a couple of pieces of code to query the database. For example, to look for uniparc IDs we can do this:

### UniParc IDs

Search a file of UniParc IDs and see if we can map this. This does not, yet, combine them into subsystems, but should!

```
python /home/edwa0468/GitHubs/CFHackathon2023/functions/find_uniparc.py -d database/functions.sqlite -f ~/UniRef/missing_ids.txt > uniparc_found.txt
```

### Roles to subsystems

Given a set of FIG functional roles in a file, extract the subsystems information. Note that this can take a column argument, if you have `tsv` data you can extract a specfic column for the role.

For example, this will take the third column (because it is 0 indexed) and append the subsystems

```
python /home/edwa0468/GitHubs/CFHackathon2023/functions/subsystems.py -f patric_functions/CFMag_bin_1.functions -c 2
```


Here is how to take a bunch of genbank files and add subsystem information

```
for i in $(seq 1 104); do echo $i; b=$((i-1)); if [ -e patric/CFMag_bin_${b}.gb ]; then genbank_to -g patric/CFMag_bin_${b}.gb -f patric_functions/CFMag_bin_${b}.functions; python /home/edwa0468/GitHubs/CFHackathon2023/functions/subsystems.py -f patric_functions/CFMag_bin_${b}.functions -c 2 > patric_functions/CFMag_bin_${b}.functions.subsystems; fi; done
```


