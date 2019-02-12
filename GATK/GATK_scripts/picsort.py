import os, subprocess

root_path= "/work/jwalters/amongue/PopGen/Manduca/NC_samples/raw_data/"

sampledirnames = [
"ACA",
"Q6",
"R10"
#"S32",
#"S33",
#"S34",
#"S35",
#"S36",
#"S37",
#"S38",
#"S39",
#"S40",
#"S42",
#"S44",
#"S45"
]

#now write picard scripts#
for name in sampledirnames:
    os.chdir(root_path + name)
    file = open (name + "picardsort.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 4\n")
    file.write("#SBATCH --mem=10G\n")
    file.write("#SBATCH -t 24:00:00\n")   
    file.write("java -jar /tools/cluster/6.2/picard-tools/1.87/SortSam.jar INPUT=" + root_path + name+ "/" + name +  ".bam OUTPUT=" + root_path + name + "/" + name + ".pic.sorted.bam SORT_ORDER=coordinate")
    file.close()
    args=['sbatch', name + "picardsort.sh"]
    subprocess.call(args)