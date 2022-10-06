####  object2
### abundance ratio
## list that contains time point  bulk-demutiplex output
# ips_total
#  $ips_day7
#   $ips_day14
#
#

### genelist
## vector that contains gene modules 
#c("CD14","CD34","NANOG")
### expression
## tpm or counts of RNA expression matrix
### composition
## time point（batch）* samples matrix
### replicate
## timepoint(sample) * gene matrix

setClass("gene_set_m",slots=list(name="character",abundance_ratio="list",genelist="vector",expression="matrix",composition="matrix",replicate="matrix",timepoint="vector",samples="vector",results='list',results_expression="list"),)

new_gene_set_m<-function(names,abundance_ratio1,genelist1,expression1,composition1,replicate1,list1){
  genelist2<-unique(genelist1)
  genelist2<-genelist2[genelist2 %in%rownames(expression1)]
  module<-new(Class = "gene_set_m",name=names,abundance_ratio=abundance_ratio1,expression=expression1,genelist=genelist2,composition=composition1,results=vector(mode="list", length=length(genelist2)),results_expression=vector(mode="list", length=length(genelist2)))
  module@samples<-colnames(module@composition)
  module@timepoint<-rownames(module@composition)
  names(module@results)<-genelist2
  names(module@results_expression)<-genelist2
  colnames(module@expression)<-module@timepoint
  return(module)
  ###
  ### hint for unmatched rownaes/colnames
  
  
  
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
gene_module_ratio_add<-function(object2){
list_gene<-object2@genelist
tem_object<-object2
for (i in list_gene ){
  tem_object@results[[i]]<-return_time_order_one_gene(tem_object@abundance_ratio,i)
  
}
  return(tem_object)
  
}
gene_module_ratio_NAremove<-function(object2,pcut_cutoff,retortN)
{
  nr= dim(object2@results[[1]])[1]
  nc= dim(object2@results[[1]])[2]
  tp=object2@timepoint
  tem_object2=object2
  #compare composion with  
  list_gene<-object2@genelist
  for  (x in list_gene ) 
  {
    for(y in tp)
    {
      if(is.na(tem_object2@results[[x]][y,1]))
      {
        tem_object2@results[[x]][y,]=tem_object2@composition[y,] 
        next;
      } 
      else 
      {
        if (tem_object2@abundance_ratio[[y]][x,"p"]>=pcut_cutoff)
        {
          
          tem_object2@results[[x]][y,]=tem_object2@composition[y,]
          next;
        } 
        else if (tem_object2@abundance_ratio[[y]][x,"p"]<pcut_cutoff) 
        {  
          for(i in 1:nc)
          { 
            v1=tem_object2@results[[x]][y,i]/tem_object2@composition[y,i]
            v2=tem_object2@composition[y,i]/tem_object2@results[[x]][y,i]
            if (v1>retortN|v2>retortN)
            {
              tem_object2@results[[x]][y,i]<-tem_object2@composition[y,i]
              next;
            }
            
          }
        }
      } 
    }
  } 
  
  
  return(tem_object2)   
}

gene_abundance_result<-function(object2)
{
  list_gene<-object2@genelist
  tem_object<-object2
  for (i in list_gene )
  {
    tem_object@results_expression[[i]]<-tem_object@expression[i,]*object2_mod2@results[[i]]
  }
  return(tem_object)
  
}





