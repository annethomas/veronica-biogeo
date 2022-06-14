library(mvMORPH)

TTR_dir="~/biogeobears/TTR"
traitset=commandArgs(trailingOnly=TRUE)[1] ##TTR4, TTR9
model=commandArgs(trailingOnly=TRUE)[2] ##BM, OU, EB

input=readRDS(file.path(TTR_dir,"mvmorph_input",paste0("mvmorph_input_",traitset,".Rdata")))
model_list=list()

model_list=list()
print(model)
print(Sys.time())
if(model=="BM"){
  mvfit=mvBM(input$traits_phy$tree,input$traits_matrix,model="BM1",echo=TRUE)
} else if(model == "OU"){
  mvfit=mvOU(input$traits_phy$tree,input$traits_matrix,model="OU1",echo=TRUE)
} else if(model == "EB"){
  mvfit=mvEB(input$traits_phy$tree,input$traits_matrix,echo=TRUE,control = list(maxit = 100000))
} else{
  stop(paste("don't recognize model",model))
}
print(Sys.time())

model_list$mvfit=mvfit
    
mvanc=estim(input$traits_phy$tree,input$traits_matrix,mvfit,asr = TRUE)
row.names(mvanc$estimates) <- input$traits_phy$tree$node.label
colnames(mvanc$estimates) <- colnames(input$traits_matrix)
model_list$mvanc=mvanc

saveRDS(model_list,file.path(TTR_dir,"mvmorph_output",paste0("mvmorph_output_",traitset,model,".Rdata")))

