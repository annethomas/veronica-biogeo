library(mvMORPH)

TTR_dir="~/biogeobears/TTR"
traitset=commandArgs(trailingOnly=TRUE)[1] ##TTR4, TTR9
model=commandArgs(trailingOnly=TRUE)[2] ## BM, OU
dirname=commandArgs(trailingOnly=TRUE)[3] ## e.g. sortadate50_merged_TTRtips
bsm_dirname=commandArgs(trailingOnly=TRUE)[4] ## e.g. BSMs_run1
BSM_range=commandArgs(trailingOnly=TRUE)[5] ## BSM number(s) comma separated in quotes

outdir=file.path(TTR_dir,"mvmorph_output",dirname)
dir.create(outdir,recursive=TRUE)

print(BSM_range)
BSM_range=as.numeric(unlist(strsplit(BSM_range,",")))
print(BSM_range)
for(i in BSM_range[1]:BSM_range[2]){
  input=readRDS(file.path(TTR_dir,"mvmorph_input/mvmorph_simmap_input",dirname,bsm_dirname,paste0("mvmorph_siminput_",traitset,"_BSM",i,".Rdata")))
  model_list=list()

  print(model)
  print(Sys.time())
  if(model=="BM"){
    mvfit=mvBM(input$simmap,input$traits_matrix,model="BMM",echo=TRUE)
  } else if(model == "OU"){
    mvfit=mvOU(input$simmap,input$traits_matrix,model="OUM",echo=TRUE)
  } else{
    stop(paste("don't recognize model",model))
  }
  print(i)
  print(Sys.time())
  #saveRDS(mvfit,file.path(TTR_dir,paste0("mvmorph_mvfit_",traitset,model,"M.Rdata")))
  model_list$mvfit=mvfit
      
  mvanc=estim(input$simmap,input$traits_matrix,mvfit,asr = TRUE)
  row.names(mvanc$estimates) <- input$simmap$node.label
  colnames(mvanc$estimates) <- colnames(input$traits_matrix)
  model_list$mvanc=mvanc
  
  saveRDS(model_list,file.path(outdir,paste0("mvmorph_output_",traitset,model,"M_BSM",i,".Rdata")))

}




