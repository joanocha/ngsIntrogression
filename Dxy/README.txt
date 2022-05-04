# before running snakefile Dxy_new.py you need to ensure the files have the same number of sites and are intersectable
# 1. Please run script to generate intersection file
python intersectMafsByPosition.py file1 file2 file3 ...
# 2. Make sure you sort intersection.txt
cat intersection.txt | sort -V -k1,1 -k2,2 > intersection.sorted.txt
# 3. Run filter_to_intersection.py to ensure you filter your maf files
python filter_to_intersection.py file1 intersection.sorted.bed
# do it for each of your files

# then you have the maf files ready to run snakefile that runs Dxyz on windows for genolikelihoods


