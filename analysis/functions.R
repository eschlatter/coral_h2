priorandpost_plot <- function(model,prior,xlim=1){
  # model: mcmcglmm model object
  # prior: the object given to mcmcglmm
  # xlim: max value on x-axis
  
  xvals = seq(0,xlim,length.out=nrow(model$VCV))  
  param_list = colnames(model$VCV)
  figs = list()
  
  for(i in 1:length(param_list)){
    #get prior parameters for the current variance component
    if(i==length(param_list)) {prior_params = prior$R #last one is the residual
    } else prior_params = prior$G[[i]] #others are elements of G
    
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
    post_samples = model$VCV[,i]
    post_density = density(post_samples)
    #post_density$y = post_density$y/post_density$n
    post_mode = post_density$x[which.max(post_density$y)]
    
    #do some scaling if necessary
    multfactor=1
    if((max(prior_density,na.rm=TRUE)*10) < max(post_density$y)){ #if the prior density is on a much smaller scale than the posterior
      multfactor = 0.4*(max(post_density$y)/max(prior_density,na.rm=TRUE))
      prior_density = prior_density * multfactor
    }
    
    colors=c('Prior'='blue','Posterior'='red')
    
  g <- ggplot(as.data.frame(model$VCV), aes(x=as.data.frame(model$VCV)[,i])) + 
        theme_light()+
        geom_density(aes(color='Posterior'))+
        xlim(0,xlim)+
        geom_line(aes(x = xvals, 
                      y = prior_density,
                      color='Prior'))+
        labs(color='Legend')+
        scale_color_manual(values=colors)+
        scale_y_continuous(sec.axis = sec_axis(~./multfactor)) +
        theme(axis.text.y.right=element_text(color='blue',margin=margin(l=-1,r=-1,unit='pt')),
              axis.text.y.left=element_text(color='red',margin=margin(l=-1,r=-1,unit='pt')),
              #axis.text.y = element_text(size = 5),
              axis.ticks = element_blank(),
              axis.title.y = element_blank(),
              #axis.text.x=element_text(size=6),
              legend.position = 'none')+
        xlab(NULL)+
        ylab(NULL)+
        ggtitle(as.character(param_list[i]))+
        theme(aspect.ratio = 1)
   
    figs[[i]] <- g
  } #variance components i
  
  patchwork::wrap_plots(figs, nrow = 1)
  return(figs)
}