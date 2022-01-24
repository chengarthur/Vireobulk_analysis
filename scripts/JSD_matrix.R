JSD_matrix<-function(data.frame1){
  m1<-as.matrix(data.frame1)
  m0<-matrix(data = NA,nrow =length(colnames(m1)),ncol = length(colnames(m1)) )
  colnames(m0)<-colnames(m1)
  rownames(m0)<-colnames(m1)
  for (i in 1:length(colnames(m1))){
    for (j in 1:length(colnames(m1))){
      
      m0[i,j]<-as.numeric(JSD(rbind(m1[,i], m1[,j])))
    }
    
    
  }
  return(m0)
}