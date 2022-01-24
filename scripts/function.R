ConvertoPdFormat<-function(csvdf){
  df<-data.frame("variants"=paste(csvdf$Chr,csvdf$Start,csvdf$Ref,csvdf$Alt,sep = "_"),"function"=csvdf$Func.refGene,"genes"=csvdf$Gene.refGene,"avsnp"=csvdf$avsnp150)
  return(df)
}
StaticsPD<-function(df){
  
  flist<-unique(df$function.)
  df_s<-matrix(data = NA,nrow = 2,ncol = length(flist))
  for ( i in 1:length(flist))
  {
    df_s[1,i]=sum(df$function.==flist[i])
    
  }
  for ( i in 1:length(flist))
  {
    df_s[2,i]=df_s[1,i]/sum(df_s[1,])
    
  }
  
  colnames(df_s)<-flist
  return(df_s)
}