"""
download HGSVC2_phased_assembly for full-resolution typing assessment

wangshuai July 11, 2022
"""


def get_link():
    f = open(file)
    out = open(command_file, 'w')
    for line in f:
        array = line.split()
        link = prefix_link + array[0]
        command = "wget %s"%(link)
        print (command, file = out)
    
    out.close()
    f.close()



prefix_link = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/assemblies/"
file = "HGSVC2_phased_assembly_download_link.txt"
command_file = "/mnt/d/HLAPro_backup/haplotype/my_HLA/download.sh"
get_link()
