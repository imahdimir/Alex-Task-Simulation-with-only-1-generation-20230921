setwd('~/Dropbox/1-GitHub/Alex-Task-Simulation-with-only-1-generation-20230921/ngen_1')

### General ggplot2 scatterplot function ###

scatter_plot = function(x,x_SE,y,y_SE,pch,pch_name,color,color_name,xlab,ylab,outfile,hlines=TRUE){
  require(ggplot2)
  require(gridExtra)
  require(reshape2)
  # Plot data frame
  plot_df = data.frame(x=x,y=y,point_type=pch,colour=color)
  # CIs
  x_lower = x+qnorm(0.025)*x_SE
  x_upper = x-qnorm(0.025)*x_SE
  y_lower = y+qnorm(0.025)*y_SE
  y_upper = y-qnorm(0.025)*y_SE
  # Plot
  splot = ggplot(plot_df,aes(x = x, y = y, color=colour, pch=point_type)) +
    geom_point()+geom_errorbarh(aes(xmin = x_lower,xmax=x_upper))+
    geom_errorbar(aes(ymin = y_lower,ymax=y_upper))+
    theme_minimal() + theme(axis.line = element_line(color="black"),
                            axis.ticks = element_line(color="black"),
                            panel.border = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust=1))+
    geom_abline(slope=1,intercept=0)+xlab(xlab)+ylab(ylab)+labs(color = color_name, pch=pch_name)
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

###################################### Read results #########################################
statistic = c('v_g','v_g_eq','v_eg','v_eg_eq','c_ge','c_ge_eq','v_y_eq',
               'r_sib','r_cousin','r_po','r_gpar')

results = data.frame(r_y = rep(c(0, 0.25, 0.5, 0.75),4),
                     v_indir=c(rep(0,4),rep(0.25,12)),
                     r_dir_indir=c(rep(0,8),rep(0.5,4),rep(1,4)))

pgi_results_true = cbind(results,data.frame(delta=rep(NA,dim(results)[1]),delta_SE=rep(NA,dim(results)[1]),alpha=rep(NA,dim(results)[1]),alpha_SE=rep(NA,dim(results)[1]),rk=rep(NA,dim(results)[1]),rk_SE=rep(NA,dim(results)[1]),k=rep(NA,dim(results)[1]),k_SE=rep(NA,dim(results)[1]),r=rep(NA,dim(results)[1]),r_SE=rep(NA,dim(results)[1]),h2_eq=rep(NA,dim(results)[1]),h2_eq_SE=rep(NA,dim(results)[1]),rho=rep(NA,dim(results)[1]),rho_SE=rep(NA,dim(results)[1]),alpha_delta=rep(NA,dim(results)[1]),alpha_delta_SE=rep(NA,dim(results)[1]),v_eta_delta=rep(NA,dim(results)[1]),v_eta_delta_SE=rep(NA,dim(results)[1])))
pgi_results_v1 = pgi_results_true

pop_pgi_results_true = pgi_results_true
pop_pgi_results_v1 = pgi_results_true

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
  print(i)
  if (results$v_indir[i]==0){
    simname = paste('r_y',results$r_y[i],sep='_')
  } else {
    simname = paste('v_indir',results$v_indir[i],'r_dir_indir',results$r_dir_indir[i],'r_y',results$r_y[i],sep='_')
  }
  ped = read.table(paste(simname,'ped',sep='.'),sep='',header=T,stringsAsFactors = F)
  # Variance component results
  vcs = read.table(paste(simname,'_VCS.txt',sep=''),header=T,stringsAsFactors = F)
  ngen = dim(vcs)[1]
  if (results$v_indir[i]==0){
    results[i,3+1:11] = c(vcs[1,'v_g'],vcs[ngen,'v_g'],0,0,0,0,vcs[ngen,'v_y'],
                         vcs[ngen,'r_delta'],0,0,0)
  } else {
    results[i,3+1:11] = c(vcs[1,'v_g'],vcs[ngen,'v_g'],
                          vcs[1,'v_eg'],vcs[ngen,'v_eg'],
                          vcs[1,'c_ge'],vcs[ngen,'c_ge'],
                          vcs[ngen,'v_y'],
                          vcs[ngen,'r_delta'],vcs[ngen,'r_eta'],
                          vcs[ngen,'r_delta_eta_c'],vcs[ngen,'r_delta_eta_tau'])
  }
  
  results[i,3+12:15] = compute_corrs(ped)
  
  # PGI results
  pgi_results_true[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'pgi_v0.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
  pgi_results_v1[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'pgi_v1.am_adj_pars.txt',sep='_'),header=T,row.names=1)))

  if (results$v_indir[i]>0){
    pop_pgi_results_true[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'population_pgi_v0.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
    pop_pgi_results_v1[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'population_pgi_v1.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
  }
  
  print(simname)
}

save.image('sim_results.RData')

write.csv(results,'results_table.csv',quote=F,row.names=F)

write.csv(pgi_results_true,'true_direct_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pgi_results_v1,'v1_direct_effect_pgi_results_table.csv',quote=F,row.names=F)

write.csv(pop_pgi_results_true,'true_pop_direct_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pop_pgi_results_v1,'v1_pop_effect_pgi_results_table.csv',quote=F,row.names=F)


#################################### Plot Results ###################################

load('sim_results.RData')


# plots the inferred r_delta against the true r_delta with standard errors (plot a in fig 4 in manuscript)
v1_r = scatter_plot(pgi_results_true$r,pgi_results_true$r_SE,
                    pgi_results_v1$r,pgi_results_v1$r_SE,
                    pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indir',
                    color=pgi_results_true$r_y,color_name='phenotype correlation',
                    xlab='True r (r_delta)',ylab='Estimated r (noise=1) (r^_{delta})',
                    outfile='plots/r_v1.png')

# plot b in fig 4 in manuscript
v1_h2_eq = scatter_plot(results$v_g_eq[results$r_y>0]/results$v_y_eq[results$r_y>0],0,
                        pgi_results_v1$h2_eq[results$r_y>0],pgi_results_v1$h2_eq_SE[results$r_y>0],
                        pch=results$v_indir[results$r_y>0]==0,pch_name='v_indir==0',
                        color=results$r_y[results$r_y>0],color_name='phenotype correlation',
                        xlab='True h2_eq',ylab='Estimated h2_eq (noise=1)',
                        outfile='plots/h2_eq_v1.png')

# plot c in fig 4 in manuscript
v1_alpha_delta = scatter_plot(pgi_results_true$alpha_delta,pgi_results_true$alpha_delta_SE,
                              pgi_results_v1$alpha_delta,pgi_results_v1$alpha_delta_SE,
                              pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                              color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                              xlab='True alpha_delta',ylab='Estimated alpha_delta (noise=1)',
                              outfile='plots/alpha_delta_v1.png')

# plot d in fig 4 in manuscript
v1_eta_delta = scatter_plot(pgi_results_true$v_eta_delta,pgi_results_true$v_eta_delta_SE,
                            pgi_results_v1$v_eta_delta,pgi_results_v1$v_eta_delta_SE,
                            pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                            color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                            xlab='True v_eta_delta',ylab='Estimated v_eta_delta (noise=1)',
                            outfile='plots/v_eta_delta_v1.png')


