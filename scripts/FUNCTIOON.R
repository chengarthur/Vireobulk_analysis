stat_output_matrix<-function(df_ma,string_type,cutoff)
{ stat_ma<-matrix(nrow = 512,ncol =4,data = 0 )
  colnames(stat_ma)<-c("nSNP","sig_counts","total","percentage")
  stat_ma[,1]<-c(1:512)
if(string_type!="COR")
  {for (i in 1:512)
    {stat_ma[i,2]=sum(ma[ma[,"nSNP"]==i,string_type]<=cutoff) 
    stat_ma[i,3]=sum(df_ma[,"nSNP"]==i)    
    stat_ma[i,4]=stat_ma[i,2]/stat_ma[i,3]
    }
  }
if(string_type=="COR")
  {for (i in 1:512)
    {stat_ma[i,2]=sum(ma[ma[,"nSNP"]==i,string_type]>=cutoff) 
    stat_ma[i,3]=sum(df_ma[,"nSNP"]==i)    
    stat_ma[i,4]=stat_ma[i,2]/stat_ma[i,3]
  
    }
  }
return(stat_ma)
}

binsize_output<-function(stat_ma1)
{ bin_ma<-matrix(nrow = 9,ncol =4,data = 0 )
  colnames(bin_ma)<-c("bin","sig_counts","total","percentage")
  for(i in 1:9)
    
    {
    bin_ma[i,1]=i
    tem_sum_sig=0
     tem_sum_total=0
     if(2^(i+1)<=length(colnames(stat_ma1)))
        {
        for (j in 2^(i-1): 2^(i+1))
          {tem_sum_sig=tem_sum_sig+stat_ma1[j,2]
          tem_sum_total=tem_sum_total+stat_ma1[j,3]    
          }
          bin_ma[i,2]=tem_sum_sig
          bin_ma[i,3]=tem_sum_total
        }
     if(2^(i+1)>length(colnames(stat_ma1)))
       {
       for (k in 2^(i-1): 2^(i))
       {tem_sum_sig=tem_sum_sig+stat_ma1[k,2]
       tem_sum_total=tem_sum_total+stat_ma1[k,3]    
       }
       bin_ma[i,2]=tem_sum_sig
       bin_ma[i,3]=tem_sum_total
       }
  }
  bin_ma[,4]=bin_ma[,2]/bin_ma[,3]
  return(bin_ma)
}

cutoffcor<-function(df_demutix,df_umi,meta,string_type,cutoff){
  if (identical(rownames(df_demutix),rownames(df_umi))==FALSE){
    print("please input identical rowname matrix")
    break
  }
  else{
    if(string_type!="JSD" && string_type!="Pvalue"){
      tem_demuti<-as.data.frame(df_demutix[meta[,string_type]>=cutoff,])
      tem_umi<-as.data.frame(df_umi[meta[,string_type]>=cutoff,])
      tem_umi$gene<-rownames(tem_umi)
      tem_demuti$gene<-rownames(tem_demuti)
      
      tem_umi<-gather(tem_umi, key = "variable", value = "value",-gene)
      tem_demuti<-gather(tem_demuti, key = "variable", value = "value",-gene)
      mel_tem<-cbind(tem_demuti,tem_umi)
      mel_tem<-mel_tem[,c(1,2,3,6)]
      colnames(mel_tem)<-c("gene","donor","predicted","umi")
      return(mel_tem)
    }
    if(string_type=="JSD"||string_type=="Pvalue"){
      tem_demuti<-as.data.frame(df_demutix[meta[,string_type]<=cutoff,])
      tem_umi<-as.data.frame(df_umi[meta[,string_type]<=cutoff,])
      tem_umi$gene<-rownames(tem_umi)
      tem_demuti$gene<-rownames(tem_demuti)
      
      tem_umi<-gather(tem_umi, key = "variable", value = "value",-gene)
      tem_demuti<-gather(tem_demuti, key = "variable", value = "value",-gene)
      mel_tem<-cbind(tem_demuti,tem_umi)
      mel_tem<-mel_tem[,c(1,2,3,6)]
      colnames(mel_tem)<-c("gene","donor","predicted","umi")
      return(mel_tem)
    
  }
}
}


meta_generate<-function(demuti_ex,n_donors,umi_count){
  inter_gene<-unique(intersect(rownames(umi_count),rownames(demuti_mat)))
  umi_tem<-as.data.frame(umi_count[inter_gene,])
  demuti_tem<-demuti_ex[inter_gene,1:n_donors]
  demuti_tem2<-demuti_ex[inter_gene,]
  
  
  ma<-matrix(data = 0, nrow = length(inter_gene), ncol = 4)
  rownames(ma)<-rownames(demuti_tem)
  for (i in 1:length(inter_gene)){
    ma[i,1]<-JSD(as.matrix(rbind(demuti_tem[i,],umi_tem[i,])))
    ma[i,2]<-cor(x = as.numeric(demuti_tem[i,]),y = as.numeric(umi_tem[i,]))
  }
  
  colnames(ma)<-c("JSD","COR","Pvalue","nSNP")
  ma<-as.data.frame(ma)
  identical(rownames(ma),rownames(demuti_tem))
  ma$Pvalue<-demuti_tem2$p
  ma$nSNP<-demuti_tem2$numSNP
  return(ma)
  
}



  