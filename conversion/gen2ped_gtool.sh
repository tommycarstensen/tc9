## convert gen to ped and map
## http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html#Sample_File_Format_
gtool -G --g out_BEAGLE/sep/BeagleOutput.22.gen --s BEAGLE22.sample --ped gtool.BEAGLE.22.ped --map gtool.BEAGLE.22.map --phenotype phenotype_1 --threshold 0.9
