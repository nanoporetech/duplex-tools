# Prerequisites for this script

Clone the repository and enter it
```
git clone https://github.com/nanoporetech/read_fillet
cd read_fillet
```

Create an environment and activate it. Install required packages
```
python3.6 -m venv venv --prompt 'read_fillet'
source ./venv/bin/activate
python -m pip install numpy edlib fire pathlib pyfastx natsort tqdm
```

Finally, run the script like this:
```
python split_reads_on_adapter.py path_to_fastq_pass_directory --pattern "*.fastq" 

Optional flags:
--n_bases_to_mask_tail (default 14)
--n_bases_to_mask_head (default 5)
--degenerate_bases (default 11)
--pattern the pattern a fastq has to match to be included (default "*.fastq.gz")
--debug (will output a .fasta with the matches)
--type Native or PCR (PCR requires the existence of certain PCR primers in the sequences)
```

That will, for each fastq, result in a file with an additional "_split" suffix (being output in the same directory) like this:
```
my_basecalls.fastq
my_basecalls_split.fastq
```

The new *_split.fastq will contain two new reads for each read that was split, now with suffix _1 and _2
Reads which were not split will also be added to the new file, so that *_split.fastq can be used as a match for any downstream analysis
```
@<read_id> <remaining_headers>
<sequence>
+
<quality>
@<read_id>_1 <remaining_headers>
<sequence>
+
<quality>
@<read_id>_2 <remaining_headers>
<sequence>
+
<quality>

```


For the assessment script, pomoxis, samtools and seqkit is necessary, so here's how to install conda with those requirements

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p
source ~/miniconda3/bin/activate
conda create --name read_fillet_env
conda activate read_fillet_env
conda install python=3.6
conda install -c bioconda seqkit samtools pomoxis
```
