library(RPANDA)

## paths
source('~/Cambridge/Research/scripts/biogeography_github/rpanda/plot_fit_env_custom_rpanda.R')
data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
bgb_dir=file.path(data_dir,"biogeobears")
TTR_dir=file.path(data_dir,"TTR")
paleo_dir=file.path(data_dir,"NZ_geo/NZ_paleoelevation")

phy=read.tree(file.path(bgb_dir,"sortadate50_tree.tre"))

### elevation data
load(file.path(paleo_dir,"/NZalps_elev.Rdata"))
plot(elevation_data,xlab="Time (Ma)", ylab="Elevation (masl)",cex.lab=1.5,cex.axis=1.5)
lines(smooth_spline_fit,lwd=2)

## reverse axis
plot(elevation_data[,1],elevation_data[,2],xlim=rev(range(elevation_data[,1])),xlab="Time (Ma)", ylab="Elevation (masl)",cex.lab=1.5,cex.axis=1.5)
lines(smooth_spline_fit,lwd=2,xlim=rev(range(elevation_data[,1])))

## linear interpolation
ss=smooth_spline_fit
x=seq(0,18,length.out=1000)
linear_elev=approx(ss$x, ss$y, xout = x)$y

linear_elev=data.frame(time=x,linear_elev=linear_elev)
plot(linear_elev)
saveRDS(linear_elev,file.path(paleo_dir,"linear_interpolation_elevation.Rdata"))

rpanda_models=list()
####################
## fit bd to time ## 
####################

# exponential speciation
f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
# linear speciation
f.lamb.lin<-function(t,y){y[1]+y[2]*t}
lamb_par_init<-c(0.05,0.01)

# zero extinction (Yule)
f.mu.z<-function(t,y){0}
mu_par_init.z<- c()
# constant extinction
f.mu.c<-function(t,y){y[1]}
mu_par_init.c<-c(0.005)

tot_time<-max(node.age(phy)$ages)

## exponential speciation, zero extinction
res.yule.exp<-fit_bd(phy,tot_time,f.lamb.exp,f.mu.z,lamb_par_init,
                   mu_par_init.z,f=105/124,expo.lamb=TRUE,fix.mu=TRUE) #don't think mu_par_init is actually necessary
plot_fit_bd(res.yule.exp,tot_time)
res.yule.exp
rpanda_models$res.yule.exp=res.yule.exp

## linear speciation, zero extinction
res.yule.lin<-fit_bd(phy,tot_time,f.lamb.lin,f.mu.z,lamb_par_init,
                   mu_par_init.z,f=105/124,fix.mu=TRUE)
plot_fit_bd(res.bd.lin,tot_time)
res.yule.lin
rpanda_models$res.yule.lin=res.yule.lin

## exponential speciation, constant extinction
res.bd.exp<-fit_bd(phy,tot_time,f.lamb.exp,f.mu.c,lamb_par_init,
                   mu_par_init.c,f=105/124,expo.lamb=TRUE,cst.mu=TRUE)
plot_fit_bd(res.bd.exp,tot_time)
res.bd.exp
rpanda_models$res.bd.exp=res.bd.exp
  
## linear speciation, constant extinction
res.bd.lin<-fit_bd(phy,tot_time,f.lamb.lin,f.mu.c,lamb_par_init,
                   mu_par_init.c,f=105/124,cst.mu=TRUE)
plot_fit_bd(res.bd.lin,tot_time)
res.bd.lin
rpanda_models$res.bd.lin=res.bd.lin

###########################
## fit to elevation data ##
###########################

# exponential speciation
f.lamb.exp<-function(t,x,y){y[1]*exp(y[2]*x)}
# linear speciation
f.lamb.lin<-function(t,x,y){y[1]+y[2]*x}
lamb_par_init<-c(0.04,0.001)

# constant extinction
f.mu.c<-function(t,x,y){y[1]}
mu_par_init<-c(0.005)
# zero extinction
f.mu.z<-function(t,x,y){0}
mu_par_init<- c()

tot_time<-max(node.age(phy)$ages)

## exponential speciation, zero extinction
res.env.exp<-fit_env(phy,linear_elev,tot_time,f.lamb.exp,f.mu.z,
                     lamb_par_init,f=105/124,fix.mu=TRUE,dt=1e-3)
res.env.exp
rpanda_models$res.env.exp=res.env.exp

## linear speciation, zero extinction
res.env.lin<-fit_env(phy,linear_elev,tot_time,f.lamb.lin,f.mu.z,
                     lamb_par_init,f=105/124,fix.mu=TRUE,dt=1e-3)
res.env.lin
rpanda_models$res.env.lin=res.env.lin

## exponential speciation, constant extinction
res.env.exp.mu<-fit_env(phy,linear_elev,tot_time,f.lamb.exp,f.mu.c,
                        lamb_par_init,mu_par_init,f=105/124,dt=1e-3)
res.env.exp.mu
rpanda_models$res.env.exp.mu=res.env.exp.mu

## linear speciation, constant extinction
res.env.lin.mu<-fit_env(phy,linear_elev,tot_time,f.lamb.lin,f.mu.c,
                        lamb_par_init,mu_par_init,f=105/124,dt=1e-3)
res.env.lin.mu
rpanda_models$res.env.lin.mu=res.env.lin.mu

#################
## save models ##
#################
saveRDS(rpanda_models,file=file.path(paleo_dir,"rpanda_models.Rdata"))
rpanda_res=data.frame(model=character(),logL=numeric(),aicc=numeric(),
                      lambda1=numeric(),lambda2=numeric(),
                      mu=numeric())
for(m in 1:length(rpanda_models)){
  print(m)
  model=rpanda_models[[m]]
  model.name=names(rpanda_models)[m]
  print(model.name)
  logL=model$LH
  aicc=model$aicc
  lambda1=model$lamb_par[1]
  lambda2=model$lamb_par[2]
  mu=model$mu_par[1]
  if(is.null(mu)) mu=NA
  newrow=data.frame(model.name=model.name,logL=logL,aicc=aicc,lambda1=lambda1,lambda2=lambda2,mu=mu)
  rpanda_res=rbind(rpanda_res,newrow)
}
write.csv(rpanda_res,file=file.path(paleo_dir,"rpanda_model_params.csv"),row.names = FALSE)
rpanda_res[rpanda_res$aicc==min(rpanda_res$aicc),]
saveRDS(res.env.lin,file=file.path(paleo_dir,"rpanda_elev/lin_elev_model.Rdata"))

#####################
## plot best model ##
#####################
res.env.lin = readRDS(file.path(paleo_dir,"rpanda_elev/lin_elev_model.Rdata"))
### elevation data
load(file.path(paleo_dir,"/NZalps_elev.Rdata"))

pdf(file.path(paleo_dir,"rpanda_elev/elevation_spec_rate_suppfig.pdf"),width=7,height=4)
layout(matrix(1:2,ncol=2))

plot(elevation_data[,1],elevation_data[,2],xlim=rev(range(elevation_data[,1])),xlab="Time (Ma)", ylab="Elevation (masl)",cex.lab=1.2,cex.axis=1.2)
lines(smooth_spline_fit,lwd=2,xlim=rev(range(elevation_data[,1])))
mtext("(a)",side=3,line=.75,adj=-.25,cex=1.2)

linear_elev=readRDS(file.path(paleo_dir,"linear_interpolation_elevation.Rdata"))

phy=read.tree(file.path(bgb_dir,"sortadate50_tree.tre"))
tot_time<-max(node.age(phy)$ages)

# plot_fit_env(res.env.lin, linear_elev,tot_time)
# plot_fit_env_custom(res.env.lin, linear_elev,tot_time,env.xlab = "Elevation (masl)")
plot_fit_env_custom(res.env.lin, linear_elev,tot_time,spec.rate.only=TRUE)
mtext("(b)",side=3,line=.75,adj=-.25,cex=1.2)
dev.off()

res.env.lin
# model :
#   [1] "environmental birth death"
# 
# LH :
#   [1] -177.8768
# 
# aicc :
#   [1] 359.8713
# 
# lamb_par :
#   [1] -0.9483144992  0.0005132614
# 

