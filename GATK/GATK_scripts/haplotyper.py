import os, subprocess

root_path= "/work/jwalters/amongue/PopGen/Manduca/NC_samples/raw_data/"

sampledirnames = [
#"ACA",
#"Q6",
#"R10",
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
    file = open (name + "haplotyper.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 12\n")
    file.write("#SBATCH --mem=25G\n")
    file.write("#SBATCH -t 80:00:00\n")
    file.write("module load GATK/3.7-Java-1.8.0_102 \n")
    file.write("java -Xmx25g -jar /nfs/apps/7/arch/generic/GATK/3.7-Java-1.8.0_102/GenomeAnalysisTK.jar \\\n")
    file.write("-nct 12 \\\n")
    file.write("--emitRefConfidence GVCF \\\n")
    file.write("-T HaplotypeCaller \\\n")
    file.write("-I " + root_path + name + "/" + name + ".SHIN.realigned.bam \\\n")
    file.write("-R /work/jwalters/amongue/references/references/Manduca_sexta/Manduca_OGS2.fa \\\n")
    file.write("-o " + root_path + name + "/" + name + ".g.vcf") 
    file.close()
    args=['sbatch', name + "haplotyper.sh"]
    subprocess.call(args)
