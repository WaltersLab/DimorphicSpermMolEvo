import os, subprocess

root_path= "/work/jwalters/amongue/PopGen/Manduca/NC_samples/raw_data/"

#I've commented out samples with more than one pair of reads, these scripts will be manually curated
sampledirnames = [
"ACA",
"Q6",
"R10",
"S32",
"S33",
#"S34",
"S35",
#"S36",
#"S37",
"S38",
#"S39",
"S40",
"S42",
"S44"
#"S45"
]


#now write bowtie scripts#

for name in sampledirnames:
    os.chdir(root_path + name)
    file = open (name + "align.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 6\n")
    file.write("#SBATCH --mem=25G\n")
    file.write("#SBATCH -e " + root_path + name + "/" + name + ".er\n")
    file.write("module load Bowtie2\n")
    file.write("module load SAMtools\n")
    file.write("bowtie2 --phred33 \\\n")
    file.write("--very-sensitive-local \\\n")
    file.write("-p 6 \\\n")
    file.write("-x /work/jwalters/amongue/references/references/Manduca_sexta/Manduca_OGS2 \\\n")
    file.write("-1 " + root_path + name + "/*_1.fq.gz \\\n")
    file.write("-2 " + root_path + name + "/*_2.fq.gz  \\\n")
    file.write("| samtools view -bS - > " + root_path + name + "/" + name + ".bam")
    file.close()
    args=['sbatch', name + "align.sh"]
    subprocess.call(args)
