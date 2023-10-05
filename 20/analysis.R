setwd('~/Dropbox/intergenerational/simulations/')

### General ggplot2 scatterplot function ###
scatter_plot = function(x,x_SE,y,y_SE,pch,pch_name,color,color_name,xlab,ylab,outfile,hlines=TRUE){
  require(ggplot2)
  require(gridExtra)
  require(reshape2)
  # Plot data frame
  plot_df = data.frame(x=x,y=y,point_type=pch,colour=color)
  # Plot
  splot = ggplot(plot_df,aes(x = x, y = y, color=colour, pch=point_type)) +
    geom_point()+
    theme_minimal() + theme(axis.line = element_line(color="black"),
                            axis.ticks = element_line(color="black"),
                            panel.border = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust=1))+
    geom_abline(slope=1,intercept=0)+xlab(xlab)+ylab(ylab)+labs(color = color_name, pch=pch_name)
  if (!x_SE==0){
    x_lower = x+qnorm(0.025)*x_SE
    x_upper = x-qnorm(0.025)*x_SE
    splot = splot+geom_errorbarh(aes(xmin = x_lower,xmax=x_upper))
  }
  if (!y_SE==0){
    y_lower = y+qnorm(0.025)*y_SE
    y_upper = y-qnorm(0.025)*y_SE
    splot = splot+geom_errorbar(aes(ymin = y_lower,ymax=y_upper))
  }
  if (hlines==TRUE){splot = splot+geom_hline(yintercept=c(0,1),linetype=2)}
  ggsave(outfile,plot=splot,width=7,height=5)
  return(splot)
}

########## Functions to compute correlations from pedigree ###########
compute_cousin_cor=function(pedigree,last_gen){
  ngen = as.integer(strsplit(pedigree[dim(pedigree)[1],1],'_')[[1]][1])
  parent_gen = pedigree[sapply(pedigree[,1],function(x) strsplit(x,'_')[[1]][1]==as.character(ngen-1)),]
  parent_fams = unique(parent_gen[,1])
  parent_sibs = t(sapply(parent_fams,function(x) parent_gen[parent_gen[,1]==x,'IID']))
  c1_phenotypes = t(sapply(parent_sibs[,1],function(x) last_gen[last_gen$FATHER_ID==x,'PHENO']))
  c2_phenotypes = t(sapply(parent_sibs[,2],function(x) last_gen[last_gen$MOTHER_ID==x,'PHENO']))
  return(sum(cor(c1_phenotypes,c2_phenotypes))/4)
}

compute_gpar_cor=function(pedigree,last_gen){
  gpar_ids = cbind(pedigree[match(last_gen$FATHER_ID,pedigree$IID),c('FATHER_ID','MOTHER_ID')],
                   pedigree[match(last_gen$MOTHER_ID,pedigree$IID),c('FATHER_ID','MOTHER_ID')])
  gpar_phenotypes = cbind(pedigree[match(gpar_ids[,1],pedigree$IID),'PHENO'],
                          pedigree[match(gpar_ids[,2],pedigree$IID),'PHENO'],
                          pedigree[match(gpar_ids[,3],pedigree$IID),'PHENO'],
                          pedigree[match(gpar_ids[,4],pedigree$IID),'PHENO'])
  return(sum(cor(last_gen$PHENO,gpar_phenotypes))/4)
}

compute_corrs=function(pedigree){
  ngen = strsplit(pedigree[dim(pedigree)[1],1],'_')[[1]][1]
  last_gen = pedigree[sapply(pedigree[,1],function(x) strsplit(x,'_')[[1]][1]==ngen),]
  sib_cor = cor(last_gen[seq(1,dim(last_gen)[1],2),'PHENO'],
                     last_gen[seq(2,dim(last_gen)[1],2),'PHENO'])
  #### FIX MOTHER PHENO ###
  last_gen$MOTHER_PHENO = pedigree[match(last_gen$MOTHER_ID,pedigree[,2]),'PHENO']
  last_gen$FATHER_PHENO = pedigree[match(last_gen$FATHER_ID,pedigree[,2]),'PHENO']
  po_cor = (cor(last_gen$PHENO,last_gen$FATHER_PHENO)+cor(last_gen$PHENO,last_gen$MOTHER_PHENO))/2
  cousin_cor = compute_cousin_cor(pedigree,last_gen)
  gpar_cor = compute_gpar_cor(pedigree,last_gen)
  outcors = c(sib_cor,cousin_cor,po_cor,gpar_cor)
  return(outcors)
}

##### Predicted correlation functions #####
predicted_cor=function(k,m,v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  r_out=((1+r_delta)/2)^(2*k+m+1)*v_g+
    ((1+r_eta)/2)^(2*k+m)*v_eg+
    0.5*((1+r_delta)/2)^k*((1+r_eta)/2)^k*(((1+r_delta)/2)^m+((1+r_eta)/2)^m)*c_ge
  return(r_out/v_y)
}

predicted_sib_cor=function(v_g,v_eg,c_ge,v_y,r_delta){
  return(predicted_cor(0,0,v_g,v_eg,c_ge,v_y,r_delta,0))
}

predicted_cousin_cor=function(v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  return(predicted_cor(1,0,v_g,v_eg,c_ge,v_y,r_delta,r_eta))
}

predicted_cousin_cor=function(v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  r_cousin=((1+r_delta)/2)^3*v_g+((1+r_eta)/2)^2*v_eg+((1+r_delta)*(1+r_eta)/4)*c_ge
  return(r_cousin/v_y)
}

predicted_ancestor_cor=function(k,v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  r_out=((1+r_delta)/2)^(k)*v_g+
    ((1+r_eta)/2)^(k)*v_eg+
    0.5*(((1+r_delta)/2)^k+((1+r_eta)/2)^(k-1))*c_ge
  return(r_out/v_y)
}

predicted_c_ge_eq=function(v_g,v_eg,r_delta,r_eta,r_delta_eta_c,r_delta_eta_tau){
  c_ge_eq = sqrt(2*v_g*v_eg/((1-r_delta)*(1-r_eta)))*(r_delta_eta_c+r_delta_eta_tau)
  return(c_ge_eq)
}

predicted_r_po=function(v_g,v_eg,c_ge,v_y,r_delta,r_eta,r_delta_eta_c,r_y){
  r=predicted_ancestor_cor(1,v_g,v_eg,c_ge,v_y,r_delta,r_eta)
  T_eq = (v_g+v_eg+c_ge)/v_y
  r=r+r_y*(1-T_eq)*((v_g+c_ge/2+v_eg)/(2*v_y)+r_delta_eta_c*sqrt((v_g*v_eg)/(2*v_y^2*(1+r_eta))))
  return(r)
}

###################################### Read results #########################################
statistic = c('v_g','v_g_eq','v_eg','v_eg_eq','c_ge','c_ge_eq','v_y_eq',
               'r_sib','r_cousin','r_po','r_gpar')

results = data.frame(r_y = rep(c(0, 0.25, 0.5, 0.75),4),
                     v_indir=c(rep(0,4),rep(0.25,12)),
                     r_dir_indir=c(rep(0,8),rep(0.5,4),rep(1,4)))

pgi_results_true = cbind(results,data.frame(delta=rep(NA,dim(results)[1]),delta_SE=rep(NA,dim(results)[1]),alpha=rep(NA,dim(results)[1]),alpha_SE=rep(NA,dim(results)[1]),rk=rep(NA,dim(results)[1]),rk_SE=rep(NA,dim(results)[1]),k=rep(NA,dim(results)[1]),k_SE=rep(NA,dim(results)[1]),r=rep(NA,dim(results)[1]),r_SE=rep(NA,dim(results)[1]),h2_eq=rep(NA,dim(results)[1]),h2_eq_SE=rep(NA,dim(results)[1]),rho=rep(NA,dim(results)[1]),rho_SE=rep(NA,dim(results)[1]),alpha_delta=rep(NA,dim(results)[1]),alpha_delta_SE=rep(NA,dim(results)[1]),v_eta_delta=rep(NA,dim(results)[1]),v_eta_delta_SE=rep(NA,dim(results)[1])))
pgi_results_v1 = pgi_results_true
pgi_results_v10 = pgi_results_true
pgi_results_v100 = pgi_results_true

pop_pgi_results_true = pgi_results_true
pop_pgi_results_v1 = pgi_results_true
pop_pgi_results_v10 = pgi_results_true
pop_pgi_results_v100 = pgi_results_true

results=cbind(results,data.frame(v_g=rep(NA,dim(results)[1]),
                                 v_g_eq=rep(NA,dim(results)[1]),
                                 v_eg=rep(NA,dim(results)[1]),
                                 v_eg_eq=rep(NA,dim(results)[1]),
                                 c_ge=rep(NA,dim(results)[1]),
                                 c_ge_eq=rep(NA,dim(results)[1]),
                                 v_y_eq=rep(NA,dim(results)[1]),
                                 r_delta=rep(NA,dim(results)[1]),
                                 r_eta=rep(NA,dim(results)[1]),
                                 r_delta_eta_c=rep(NA,dim(results)[1]),
                                 r_delta_eta_tau=rep(NA,dim(results)[1]),
                                 r_sib=rep(NA,dim(results)[1]),
                                 r_cousin=rep(NA,dim(results)[1]),
                                 r_po=rep(NA,dim(results)[1]),
                                 r_gpar=rep(NA,dim(results)[1])))

sim = 'r_y_0.5'

for (i in 1:dim(results)[1]){
  if (results$v_indir[i]==0){
    simname = paste('r_y',results$r_y[i],sep='_')
  } else {
    simname = paste('v_indir',results$v_indir[i],'r_dir_indir',results$r_dir_indir[i],'r_y',results$r_y[i],sep='_')
  }
  ped = read.table(paste(simname,'ped',sep='.'),sep='',header=T,stringsAsFactors = F)
  # Variance component results
  vcs = read.table(paste(simname,'_VCs.txt',sep=''),header=T,stringsAsFactors = F)
  ngen = dim(vcs)[1]
  if (results$v_indir[i]==0){
    results[i,3+1:11] = c(vcs[1,'v_g'],vcs[ngen,'v_g'],0,0,0,0,vcs[ngen,'v_y'],
                         vcs[ngen,'r_delta'],0,0,0)
