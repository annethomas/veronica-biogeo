plot_fit_env_custom<-function (fit.env, env_data, tot_time,env.xlab,cex=1.2,spec.rate.only=FALSE) 
{
  if (!inherits(fit.env, "fit.env")) 
    stop("object is not of class \"fit.env\"")
  t <- seq(0, tot_time, length.out = 100)
  #dev.new()
  # plot(-t, fit.env$f.lamb(t), type = "l", xlab = "Time (Ma)", 
  #      ylab = "Speciation rate", main = "Fitted speciation rate",cex.lab=cex,cex.axis=cex)
  
  ###############################
  ### the target plot for ms ####
  ###############################
  
  plot(t,xlim=rev(range(t)), fit.env$f.lamb(t), type = "l", xlab = "Time (Ma)", 
       ylab = "Speciation rate (#/Ma/lineage)",cex.lab=cex,cex.axis=cex)
  df <- smooth.spline(env_data[, 1], env_data[, 2])$df
  spline_result <- pspline::sm.spline(env_data[, 1], env_data[, 2], 
                             df = df)
  env_func <- function(t) {
    predict(spline_result, t)
  }
  #dev.new()
  # plot(env_func(t), fit.env$f.lamb(t), xlab = env.xlab, 
  #      ylab = "Speciation rate", main = "Fitted speciation rate",cex.lab=cex,cex.axis=cex)

  #################
  ## other plots ##
  #################
  if(!spec.rate.only){
    plot(env_func(t), fit.env$f.lamb(t), xlab = env.xlab, 
         ylab = "Speciation rate",cex.lab=cex,cex.axis=cex)
    
    if ("f.mu" %in% attributes(fit.env)$names) {
    #dev.new()
    plot(-t, fit.env$f.mu(t), type = "l", xlab = "Time (Ma)", 
         ylab = "extinction rate", main = "Fitted extinction rate",cex.lab=cex,cex.axis=cex)
    #dev.new()
    plot(env_func(t), fit.env$f.mu(t), xlab = env.xlab, 
         ylab = "extinction rate", main = "Fitted extinction rate",cex.lab=cex,cex.axis=cex)
    r <- function(t) {
      fit.env$f.lamb(t) - fit.env$f.mu(t)
    }
    #dev.new()
    plot(-t, r(t), type = "l", xlab = "Time (Ma)", 
         ylab = "net diversification rate", main = "Fitted net diversification rate",
         cex.lab=cex,cex.axis=cex)
    #dev.new()
    plot(env_func(t), r(t), xlab = env.xlab, 
         ylab = "net diversification rate", main = "Fitted net diversification rate",
         cex.lab=cex,cex.axis=cex)
    }
    else {
    #dev.new()
    plot(-t, fit.env$f.lamb(t), type = "l", xlab = "Time (Ma)", 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
    #dev.new()
    plot(env_func(t), fit.env$f.lamb(t), xlab = env.xlab, 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
   }
  }
  
}