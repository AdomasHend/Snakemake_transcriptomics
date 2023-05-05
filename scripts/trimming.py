import subprocess
import os
def pair_reads(*args):
    pairs = {}
    for sample in args:
        if '_R1_' in sample:
            prefix1 = sample.split("_R1_")[0]
        else:
            prefix1 = sample.split("_R2_")[0]
        read1 = prefix1 + "_R1_001.fastq.gz"
        read2 = prefix1 + "_R2_001.fastq.gz"
        pairs[read1] = read2
    return pairs

pairs =  pair_reads(*snakemake.input)


for read1, read2 in pairs.items():
    basename1 = os.path.basename(read1)
    basename2 = os.path.basename(read2)
    out_ext = os.path.splitext(basename1)[0]
    out = os.path.splitext(out_ext)[0] + '.fastq'
    out2_ext = os.path.splitext(basename2)[0]
    out2 = os.path.splitext(out2_ext)[0] + '.fastq'
    subprocess.run("bash /home/adomash/Snakemake-transcriptomics/bbmap/bbduk.sh in={} in2={} ref=/home/adomash/Snakemake-transcriptomics/data/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 out=trimmed/{out} out2=trimmed/{out2}".format(read1, read2, out=out, out2=out2),shell=True)

