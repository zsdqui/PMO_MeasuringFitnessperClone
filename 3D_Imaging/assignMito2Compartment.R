# other_center_coord="0_prediction_c0.model.p_Cells_Centers.csv";
# mito_center_coord="0_target_Cells_Centers.csv"
assignMito2Compartment<-function(mito_center_coord, other_center_coord){
  setwd("~/Projects/PMO/MeasuringFitnessPerClone/data/3Dbrightfield/allencell/G03_CellposeOutput/Cells_center_coordinates")
  ZSTACK_DISTANCE=0.29
  mito_coord=read.csv(file=mito_center_coord,check.names = F,stringsAsFactors = F)
  other_coord=read.csv(file=other_center_coord,check.names = F,stringsAsFactors = F)
  mito_coord$z=mito_coord$z*ZSTACK_DISTANCE
  other_coord$z=other_coord$z*ZSTACK_DISTANCE
  
  ##@TODO: should not nbe necessary!
  mito_center_coord=strsplit(mito_center_coord,"_Cells_")[[1]][1]
  other_center_coord=strsplit(other_center_coord,"_Cells_")[[1]][1]
  
  d=flexclust::dist2(mito_coord, other_coord)
  ##Assign mitochondria to closest nucleus
  i_nuc=apply(d,1,which.min)
  ## Aggregate mitochondria all_cells_coordinates for each nucleus 
  ## and name output file according to nucleus ID
  for(i in unique(i_nuc)){
    i_mito=which(i==i_nuc)
    dm=c()
    for(j in i_mito){
      f_mito=list.files(paste0("../All_Cells_coordinates/",mito_center_coord),pattern =paste0("cell_",j,"_"),full.names = T )
      ## @TODO: should not be necessary
      if(isempty(f_mito)){
        next
      }
      coord=read.csv(file=f_mito,check.names = F,stringsAsFactors = F)
      coord$mito_id=j
      dm=rbind(dm,coord)
    }
    #### WRITE LINKED COMPARTMENTS:
    ## Copy other compartment file (e.g. nucleus)
    f_nuc=list.files(paste0("../All_Cells_coordinates/",other_center_coord),pattern =paste0("cell_",i,"_"),full.names = T )
    ## @TODO: should not be necessary
    if(isempty(f_nuc) || is.null(dm)){
      next
    }
    file.copy(f_nuc,paste0("../../G05_multiOrganelles_Linked/",other_center_coord,"_cell_",i,"_coordinates.csv"))
    ## Also write mitochondria file
    write.csv(dm,paste0("../../G05_multiOrganelles_Linked/",mito_center_coord,"_mito_cell_",i,"_coordinates.csv"), row.names = F, quote = F)
  }
}
