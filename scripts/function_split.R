cell_type_partiont<-function(seurat_ob,split_meta,ratio_meta){
  s<-split_meta
  r<-ratio_meta
  meta_raw<-seurat_ob@meta.data
  list_split<-unique(meta_raw[,s])
  list_ratio<-unique(meta_raw[,r])
  matrix_result<-matrix(data = 0,ncol = length(list_split),nrow = length(list_ratio))
  matrix_result<-as.data.frame(matrix_result)
  colnames(matrix_result)<-list_split
  rownames(matrix_result)<-list_ratio
  for (i in list_split){
    for (j in list_ratio){
      matrix_result[j,i]=sum(meta_raw[s]==i & meta_raw[r]==j)/sum(meta_raw[s]==i)
          }
  }
  
  
  return(matrix_result)
}

