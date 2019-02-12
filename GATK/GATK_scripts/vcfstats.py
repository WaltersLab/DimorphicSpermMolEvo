import os, subprocess

root_path= "/work/jwalters/amongue/PopGen/Manduca/NC_samples/raw_data/"

sampledirnames = [
#"ACA",
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
    file = open (name + "vcfstats.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 1\n")
    file.write("#SBATCH --mem=10G\n")
    file.write("#SBATCH -t 80:00:00\n")
    file.write("module load BCFtools \n")
    file.write("bcftools stats " + root_path + name + "/" + name + ".g.vcf > " + root_path + name + "/" + name + "vcf.stats") 
    file.close()
    args=['sbatch', name + "vcfstats.sh"]
    subprocess.call(args)
    
