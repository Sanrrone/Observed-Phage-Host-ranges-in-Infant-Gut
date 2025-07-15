cat supp_files/sep_sa.txt | while read h
do
cd $new/sandro/HeP_samples/2_assembly/${h}_virmining
grep "phage" minedviruses.fna | awk -F"_" -v h=$h '{gsub(">","",$0);print h"\t"$0"\t"$1"_"$2"\t"$4"\t"$5}' >> ~/phage_approach/supp_files/prophages_coordinates.tsv

grep "checkv" minedviruses.fna | awk -F"_" '{gsub(">","",$0);print $0"\t"$1"_"$2"\t"$5}' | awk -F"\t" -v h=$h '{split($3,a,"-");split(a[2],b,"/");print h"\t"$1"\t"$2"\t"a[1]"\t"b[1]}' >> ~/phage_approach/supp_files/prophages_coordinates.tsv
cd ..
done
