##
demutiplex_ratio<-cell_type_partiont(query_kidney_raw,"orig.ident","predicted.id")
ref@meta.data$batch_ip<-paste(ref$iPSC.line,ref$Protocol,sep = "_")
ref_ratio<-cell_type_partiont(ref,"batch_ip","Putative_Cell_Types")
demu_ratio<-fill_TEcc(ref_ratio,demutiplex_ratio)

demu_ratio_merge<-rbind(
  "Interstitial Cells-1"=demu_ratio["Interstitial Cells-1",]+demu_ratio["Interstitial Cells-2",]+demu_ratio["Interstitial Cells-3",]+demu_ratio["Interstitial Cells-4",]+demu_ratio["Interstitial Cells-5",],
  "Proliferating-1"=demu_ratio["Proliferating-1",],
  "Proximal Convoluted Tubule-like"=demu_ratio["Proximal Convoluted Tubule-like",],
  "Distal Nephron"=demu_ratio["Distal Nephron",],
  "Collecting Duct-like"=demu_ratio["Collecting Duct-like",],
  "Podocyte-like"=demu_ratio["Podocyte-like",],
  "Nephron Progenitor Cells-Distal"=demu_ratio["Nephron Progenitor Cells-Distal",],
  "Muscle-like"=demu_ratio["Muscle-like-1",],
  "Neuronal precursor"=demu_ratio["Neuronal precursor",],
  "Melanoma"=demu_ratio["Melanoma",],
  "Neuron-like"=demu_ratio["Neuron-like",]
)
demuti_jsd_raw_matrix<-JSD_matrix(demutiplex_ratio)
demuti_jsd_merge_matrix<-JSD_matrix(demu_ratio_merge)

###
total_ratio<-fill_TEcc(demu_ratio,experiment_organoids)
total_ratio_merge<-rbind(
  "Interstitial Cells-1"=total_ratio["Interstitial Cells-1",]+total_ratio["Interstitial Cells-2",]+total_ratio["Interstitial Cells-3",]+total_ratio["Interstitial Cells-4",]+total_ratio["Interstitial Cells-5",]+total_ratio["Interstitial Cells-6",]+total_ratio["Interstitial Cells-7",]+total_ratio["Interstitial Cells-8",],
  "Proliferating-1"=total_ratio["Proliferating-1",],
  "Proximal Convoluted Tubule-like"=total_ratio["Proximal Convoluted Tubule-like",],
  "Distal Nephron"=total_ratio["Distal Nephron",],
  "Collecting Duct-like"=total_ratio["Collecting Duct-like",],
  "Podocyte-like"=total_ratio["Podocyte-like",],
  "Nephron Progenitor Cells-Distal"=total_ratio["Nephron Progenitor Cells-Distal",],
  "Endothelial Cells"=total_ratio["Endothelial Cells",],
  "Muscle-like"=total_ratio["Muscle-like-1",]+total_ratio["Muscle-like-2",],
  "Neuronal precursor"=total_ratio["Neuronal precursor",],
  "Melanoma"=total_ratio["Melanoma",],
  "Neuron-like"=total_ratio["Neuron-like",]
)
##
total_jsd_raw_matrix<-JSD_matrix(total_ratio)
total_jsd_merge_matrix<-JSD_matrix(total_ratio_merge)