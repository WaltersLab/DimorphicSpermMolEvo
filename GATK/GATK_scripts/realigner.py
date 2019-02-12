import os, subprocess

root_path= "/work/jwalters/amongue/PopGen/Manduca/NC_samples/raw_data/"

sampledirnames = [
"ACA",
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

#now write picard scripts#
for name in sampledirnames:
    os.chdir(root_path + name)
    file = open (name + "realigner.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 4\n")
    file.write("#SBATCH --mem=10G\n")
    file.write("#SBATCH -t 24:00:00\n")
    file.write("module load GATK/3.7-Java-1.8.0_102 \n")
    file.write("java -Xmx10g -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar \\\n")
    file.write("-T IndelRealigner \\\n")
    file.write("-I " + root_path + name + "/" + name + ".rg.rmdup.sorted.bam \\\n")
    file.write("-R /work/jwalters/amongue/references/references/Manduca_sexta/Manduca_OGS2.fa \\\n")
    file.write("-targetIntervals "  + root_path + name + "/" + name + ".rmdup.sorted.bam.list \\\n")
    file.write("-o " + root_path + name + "/" + name + ".realigned.bam") 
    file.close()
    args=['sbatch', name + "realigner.sh"]
    subprocess.call(args)
    
