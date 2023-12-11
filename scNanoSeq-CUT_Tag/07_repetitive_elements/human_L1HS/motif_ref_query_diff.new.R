#!/usr/bin/Rscripts
###
 # @Descripttion:scNanoSeq-CUT&Tag LINEs comparison
 # @version:
 # @Author: Yuqing Guo
 # @Date: 2023-7-24
 # @LastEditors: Yuqing Guo
 # @LastEditTime: 2023-8-10
###

argv <- commandArgs(TRUE)
seq<- argv[1]

input_file=paste0(seq,'_ref_query_diff_final')

df=read.table(input_file,header = F,sep = '\t')
df$index=seq(1,dim(df)[1])
if (dim(df[df$V3=='INS',])[1]>0){
df_INS=df[df$V3=='INS',]
rownames(df_INS)=seq(1,dim(df_INS)[1])
df_INS$index2=seq(1,dim(df_INS)[1])
df_INS$diff=c(0,diff(df_INS$index))
df_INS$Group1='NA'
df_INS[df_INS$diff!=1,]$Group1='Split'
df_INS$Group2=df_INS$Group1
df_INS[df_INS$diff!=1,]$Group2=paste0('Split_',c(1:dim(df_INS[df_INS$diff!=1,])[1]))
df_new=df
split_num=dim(df_INS[df_INS$Group1=='Split',])[1]
if (split_num>1){
    for (i in 1:(split_num-1)){
        str_group=paste0('Split_',i)
        end_group=paste0('Split_',i+1)

        str_index=df_INS[df_INS$Group2==str_group,'index']
        end_tmp_index=df_INS[df_INS$Group2==end_group,'index2']-1
        end_index=df_INS[df_INS$index2==end_tmp_index,'index']


        sub_index=str_index-1
        sub_snp=paste(df[(str_index-1):end_index,'V2'],collapse = "")
        df_new[sub_index,'V2']=sub_snp


    }

    i=split_num
    str_group=paste0('Split_',i)

    str_index=df_INS[df_INS$Group2==str_group,'index']
    end_index=tail(df_INS,1)$index
    sub_index=str_index-1
    sub_snp=paste(df[(str_index-1):end_index,'V2'],collapse = "")
    df_new[sub_index,'V2']=sub_snp


}else{
    i=split_num
    str_group=paste0('Split_',i)

    str_index=df_INS[df_INS$Group2==str_group,'index']
    end_index=tail(df_INS,1)$index
    sub_index=str_index-1
    sub_snp=paste(df[(str_index-1):end_index,'V2'],collapse = "")
    df_new[sub_index,'V2']=sub_snp



}
 df_new=df_new[df_new$V3!='INS',]

for (i in 1:dim(df_new)[1]){

    if (df_new[i,]$V2==df_new[i,]$V1){
        df_new[i,]$V3='.'
    }else if (df_new[i,]$V2=='-'){
        df_new[i,]$V3='GAP'
    }else{
        df_new[i,]$V3=df_new[i,]$V2
    }
}
   
df_new=df_new[,1:3]
colnames(df_new)=c('Ref',seq,paste0(seq,'_Diff'))
output=paste0(seq,'_ref_query_diff_final_modify')
write.table(df_new,output,sep = '\t',quote = F,row.names = F)
    
}else{
    df_new=df[,1:3]
    colnames(df_new)=c('Ref',seq,paste0(seq,'_Diff'))
    output=paste0(seq,'_ref_query_diff_final_modify')
    write.table(df_new,output,sep = '\t',quote = F,row.names = F)
    
}
