less -S ../GIAB_131219_D00360_005_BH814YADXX/U0a_CGATGT_L001_001/U0a_CGATGT_L001_001.error0_L1HS_300bp_bin.bed|cut -f1,2,3,4,5 >L1HS_INFO

ls ../*/*/*.error0_L1HS_300bp_bin.bed |while read id 
do 
number=$(less $id|wc -l)
if [ $number != 0 ]
then

    less $id|cut -f6 >tmp
    paste -d '\t' L1HS_INFO tmp >info_tmp
    mv info_tmp L1HS_INFO

fi 
done 

less -S L1HS_INFO|awk '{for(i=6;i<NF;i++)printf("%s ",$i);print $NF}' |awk 'BEGIN{sum = 0}{for(i = 1; i <= NF; i++) {sum += $i} {print sum; sum = 0}}'|paste -d '\t' <(less L1HS_INFO |cut -f 1,2,3,4,5) - |sort -k 4,4V -k 5,5V >L1HS_INFO_sum

less L1HS_INFO_sum|awk '{if ($6>3)print $0"\t"1;else print $0"\t"0}'|sort -k 5,5V|bedtools groupby -g 5 -c 7 -o sum|awk '{print $0"\t"100-100*$2/320}'>GIAB_GM12878_L1HS_INFO_sum_result




