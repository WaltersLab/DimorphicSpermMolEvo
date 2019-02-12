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
    file = open (name + "rmdup.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 4\n")
    file.write("#SBATCH --mem=25G\n")
    file.write("#SBATCH -t 24:00:00\n")   
    file.write("java -Xmx25G -jar /tools/cluster/6.2/picard-tools/1.87/MarkDuplicates.jar MAX_RECORDS_IN_RAM=400000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000 INPUT=" + root_path + name+ "/" + name +  ".pic.sorted.bam OUTPUT=" + root_path + name+ "/" + name +  ".rmdup.sorted.bam METRICS_FILE=" + root_path + name+ "/" + name +  "met.txt REMOVE_DUPLICATES=true")
    file.close()
    args=['sbatch', name + "rmdup.sh"]
    subprocess.call(args)  
    
