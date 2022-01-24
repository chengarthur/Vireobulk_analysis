esayheatmapoutput<-function(out_form,genelist,gloden_ratio,pvalue_cutoff,ratio_cutoff,number_cutoff){
  nSamples<-length(gloden_ratio)
  templematrix<-out_form[genelist,1:nSamples]
  temp2<-out_form[genelist,-(1:nSamples)]
  temp2<-na.omit(temp2)
  templematrix<-na.omit(templematrix)
  for (i in 1: nSamples)
  { if (gloden_ratio[i]>=ratio_cutoff){
    
    templematrix[,i]<-templematrix[,i]/gloden_ratio[i]
  }
   
   else {
     
     templematrix[,i]<-templematrix[,i]/ratio_cutoff
   }
    
    
    }
  
  tem_matrix2<-cbind(templematrix,temp2)
  tem_matrix2<-tem_matrix2[tem_matrix2$p<=pvalue_cutoff,]
  tem_matrix2<-tem_matrix2[tem_matrix2$numSNP>=number_cutoff,]
  return(tem_matrix2)
  }
esayheatmapoutput_no_modified<-function(out_form,genelist,gloden_ratio){
  nSamples<-length(gloden_ratio)
  templematrix<-out_form[genelist,1:nSamples]
  temp2<-out_form[genelist,-(1:nSamples)]
  temp2<-na.omit(temp2)
  templematrix<-na.omit(templematrix)
  tem_matrix2<-cbind(templematrix,temp2)
  return(tem_matrix2)
  
}
return_time_order_one_gene<-function(list1,gene){
  namelist<-names(list1)
  ndonor<-length(colnames(as.data.frame(list1[1])))-3
  tem_matrix<-matrix(data = NA,nrow = length(namelist),ncol = ndonor)
  tem_matrix<-as.data.frame(tem_matrix)
  rownames(tem_matrix)<-namelist
  colnames(tem_matrix)<-colnames(as.data.frame(list1[1]))[1:ndonor]
  for (i in 1:length(namelist)){
    tem_matrix[i,]<-data.frame(list1[namelist[i]])[gene,1:ndonor]
    
    
  }
  return(tem_matrix)
}

formated_data<-function(tem_matrix){
  
  data<-tem_matrix
  data<-cbind(data,"days"=rownames(data))
  data$days<-mapvalues(data$days,c( "ips_p2","ips_p3","ips_day7","ips_day18","ips_day28"),c(-7,0,7,18,28))
  data$days<-as.numeric(data$days)
  long_data<-melt(data,id.vars = "days")

    return(long_data)
}