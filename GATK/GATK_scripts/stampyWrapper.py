import os, subprocess

root_path= "/work/jwalters/amongue/PopGen/Manduca/NC_samples/raw_data/"
ref_path="/work/jwalters/amongue/references/references/Manduca_sexta/Manduca_OGS2"

#here we align our other sphingids to the M.sexta reference

sampledirnames = [
"ACA",
"Q6",
"R10"
]

#now write picard scripts#
for name in sampledirnames:
    os.chdir(root_path + name)
    file = open (name + "stampy.sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH -N 1\n")
    file.write("#SBATCH -n 1\n")
    file.write("#SBATCH -c 1\n")
    file.write("#SBATCH --mem=30G\n")
    file.write("#SBATCH -t 100:00:00\n")
    file.write("module load Stampy/1.0.22\n")
    file.write("module load SAMtools\n")  
    file.write("stampy.py -g" + ref_path +  " -h" + ref_path +  " --substitutionrate=0.1 -M " + root_path + name + "/" + name + ".bam >" + root_path + name + "/" + name + ".stampy.bam")
    file.close()
    args=['sbatch', name + "stampy.sh"]
    subprocess.call(args)
