import os, subprocess

root_path= "/work/jwalters/amongue/PopGen/Manduca/NC_samples/raw_data/"

#we're excluding ACA because it's so divergent
sampledirnames = [
"Q6",
"R10",
"S32",
"S33",
"S34",
"S35",
"S36",
"S37",
"S38",
"S39",
"S40",
"S42",
"S44",
"S45"
]


#let's generate some psmc-friendly diploid fasta info
for name in sampledirnames:
    os.chdir(root_path + name)
    file = open (name + "diplo.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 1\n")
    file.write("#SBATCH --mem=35G\n")
    file.write("#SBATCH -e " + root_path + name + "/" + name + ".er\n")
    file.write("module load BCFtools\n")
    file.write("module load SAMtools\n")
    file.write("samtools mpileup -uf /work/jwalters/amongue/references/references/Manduca_sexta/Manduca_OGS2.fa \\\n")
    file.write(root_path + name + "/" + name + ".realigned.bam | \\\n")
    file.write("bcftools call -c | vcfutils.pl vcf2fq > " + root_path + name + "/" + name + ".diploid.fq.gz")
    file.close()
    args=['sbatch', name + "diplo.sh"]
    subprocess.call(args)
