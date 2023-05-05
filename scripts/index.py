import subprocess

data = snakemake.input.out_bam

for files in data:
    subprocess.run("samtools index -b {file} ".format(file=files),shell=True)
