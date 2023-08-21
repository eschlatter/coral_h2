############################################################
fn_plot_one_priorandpost <- function(prior_params,post_samples,name='variance component',xlim=1){
############################################################
  # prior_params: prior distribution parameters, like the ones passed to mcmcglmm
  # post_samples: samples from a posterior distribution formatted as a one-column mcmc object
  # name: character string of variance component name
  # xlim: max value on x-axis
  
  xvals = seq(0,xlim,length.out=length(post_samples))
  
  #generate prior density (y values to match the existing xvals vector)
  if(length(prior_params)==4){ #parameter-expanded priors
    prior_density = df(xvals/prior_params$alpha.V, df1 = 1, df2 = prior_params$nu, ncp = (prior_params$alpha.mu^2)/prior_params$alpha.V)
  } else if(length(prior_params)==2){ #regular inverse-gamma priors
    prior_density = dinvgamma(xvals,shape = prior_params$nu/2,scale = prior_params$nu*prior_params$V/2)
  } else print('Warning: cannot handle this type of prior')
  
  #get prior mode to use for plotting    
  prior_mode = xvals[which.max(prior_density)]
  if(prior_mode==xvals[length(xvals)]) print('Warning: xlim is too small')
  
  #grab posterior samples and convert to density
  prior_density_noinf <- prior_density[prior_density!=Inf]
  n_samps = 100*max((max(max(post_samples,na.rm=TRUE),max(prior_density_noinf,na.rm=TRUE)))/xlim,1)
  post_density = density(post_samples,n=n_samps) #make sure the density samples enough to match prior
  post_mode = post_density$x[which.max(post_density$y)]
  
  #do some scaling if necessary
  multfactor=1
  if((max(post_density$y) > max(prior_density_noinf,na.rm=TRUE)*10)){ #if the posterior density is on a much larger scale than the prior
    multfactor = 0.5*(max(post_density$y)/max(prior_density_noinf,na.rm=TRUE))
  }
  
  #make a plot
  colors=c('Prior'='blue','Posterior'='red')
  pd = data.frame(x=post_density$x,y=post_density$y)
  
  g <- ggplot()+
    theme_light()+
    scale_color_manual(values=colors)+
    geom_line(data=pd,aes(x=x,y=y/multfactor,color='Posterior'))+ #posterior, scaled by multfactor
    geom_line(data=data.frame(x=xvals,y=prior_density),aes(x=x,y=y,color='Prior'))+ #prior
    xlim(0,xlim)+
    scale_y_continuous(sec.axis = sec_axis(~.*multfactor))+
    theme(axis.text.y.right=element_text(color='red',margin=margin(l=-1,r=-1,unit='pt')),
          axis.text.y.left=element_text(color='blue',margin=margin(l=-1,r=-1,unit='pt')),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.x=element_text(size=6),
          #legend.position = 'none'
    )+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle(as.character(name))+
    theme(aspect.ratio = 1)
  
  return(g)
  
} #fn_plot_one_priorandpost


############################################################
fn_plot_all_priorandpost <- function(model,prior,xlim=1){
############################################################
  # model: mcmcglmm model object
  # prior: the object given to mcmcglmm
  # xlim: max value on x-axis

  param_list = colnames(model$VCV)
  figs = list()
  
  for(i in 1:length(param_list)){
    #get prior parameters for the current variance component
    if(i==length(param_list)) {prior_params = prior$R #last one is the residual
    } else prior_params = prior$G[[i]] #others are elements of G
    
    g <- fn_plot_one_priorandpost(prior_params,model$VCV[,i],as.character(param_list[i]),xlim)
    print(g)
    #figs[[i]] <- g
    
  } #variance components
  
  #patchwork::wrap_plots(figs, nrow = 1)
    
} #fn_plot_all_priorandpost


###########################################################
fn_herit_plot <- function(model,modeltype){
###########################################################  
  #model: mcmc object
  #modeltype: animal model or half-sib design; determines how heritability is calculated
  
  traces <- as.data.frame(model$VCV)
  
  if(modeltype=='animal'){
    herit = traces$animal/rowSums(traces)
    t = 'h^2: VA/(VA+VD+VS+VE)'
  } else if(modeltype=='halfsib'){
    herit = 2*traces$sire/rowSums(traces)
    t = 'h^2: 2*VD/(VD+VS+VD:S+VE)'
  } else print('Please specify model type: animal or halfsib')
  
  g <- ggplot(as.data.frame(herit), aes(x=as.data.frame(herit)[,1])) + 
        theme_light()+
        geom_density()+
        xlim(0,1)+
        theme(axis.text.y=element_text(color='black',size=5,margin=margin(l=-1,r=-1,unit='pt')),
              axis.ticks = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x=element_text(size=6),
              legend.position = 'none')+
        xlab(NULL)+
        ylab(NULL)+
        ggtitle(as.character(t))
  
  return(g)
} #fn_herit_plot