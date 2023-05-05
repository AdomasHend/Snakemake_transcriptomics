import subprocess
import os
data = snakemake.input.out_bam

for files in data:
    basename1 = os.path.basename(files)
    out_ext = os.path.splitext(basename1)[0]
    out = os.path.splitext(out_ext)[0] + ".counts.txt"
    print(out)
    subprocess.run("featureCounts -p -t exon -g gene_id -a /home/adomash/Snakemake-transcriptomics/data/play_data_ref_annot/chr19_20Mb.gtf -o {output} {file} -s 2 ".format(file=files,output=out),shell=True)
