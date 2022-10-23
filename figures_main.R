# This file contains the code to generate the figures from the paper given
# the data from the computation script in the data_folder.

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

# DISCLAIMER
# The data in the folder "data/" was generated with a different naming of the 
#   methods; method 1 is the same, method 2 is an alternative method not analysed
#   in the paper and method 3 is method 2 of the paper. This is only relevant
#   for the naming in the data files, the code here parses all the data s.t. it has the
#   naming corresponding to the paper (c.f. function "dict" in this file).
data_folder = file.path("data")


#---- useful plot functions ----
dict = function(x){return(switch(x,
   "mean,1,BIC" = "P.1, refined, BIC",
   "mean,1,TEST" = "P.1, refined, IRC",
   "mean,1simp,BIC" = "P.1, simple, BIC",
   "mean,1simp,TEST" = "P.1, simple, IRC",
   "mean,2,BIC" = "alt., refined, BIC",
   "mean,2,TEST" = "alt., refined, IRC",
   "mean,2simp,BIC" = "alt., simple, BIC",
   "mean,2simp,TEST" = "alt., simple, IRC",
   "mean,3,BIC" = "P.2, refined, BIC",
   "mean,3,TEST" = "P.2, refined, IRC",
   "mean,3simp,BIC" = "P.2, simple, BIC",
   "mean,3simp,TEST" = "P.2, simple, IRC",
   
   "GIES,GIES,TEST" = "GIES"
))}



# ---- figure 1a (high-dim skeleton recovery) ----
hd_skel1 = readRDS(file.path(data_folder, "1a_01.Rds"))
hd_skel1 = hd_skel1$df[hd_skel1$df$method %in% c("mean","median","ltest","baseObs"), ]
f1 = ggplot(hd_skel1, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.") + xlab("number of samples")
hd_skel2 = readRDS(file.path(data_folder, "1a_02.Rds"))
hd_skel2 = hd_skel2$df[hd_skel2$df$method %in% c("mean","median","ltest","baseObs"), ]
f2 = ggplot(hd_skel2, aes(ndatasets, SHD, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.") + xlab("number of interv. datasets")
hd_skel3 = readRDS(file.path(data_folder, "1a_03.Rds"))
hd_skel3 = hd_skel3$df[hd_skel3$df$method %in% c("mean","median","ltest","baseObs"), ]
f3 = ggplot(hd_skel3, aes(interventionSize, SHD, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.") + xlab("size of interventions")
fig1a = f1 + f2 + f3 + plot_layout(nrow = 1, guides="collect")
fig1a



# ---- figure 1b (low- and high-dim polytrees and DAGs with GIES) ----
legend_str = "Method"

# low-dim polytree
dfplotPT = readRDS(file.path(data_folder, "1b_ld_pt.Rds"))$df
rand_shd_mean_ldpt = mean(dfplotPT[dfplotPT$method == dfplotPT$method[1], "shd_rand"])
dfplotPT$method = sapply(dfplotPT$method, dict, USE.NAMES = FALSE)
f_ld_PT = ggplot(dfplotPT, aes(totalSamples, SHD, fill=factor(method))) +
  geom_boxplot() + labs(fill=legend_str, title = "Low-dim. Polytrees")+ xlab(NULL)+
  geom_hline(yintercept = rand_shd_mean_ldpt, linetype="dotted") +
  theme(plot.title = element_text(face="bold"))

# low-dim DAGs
dfplotDAG = readRDS(file.path(data_folder, "1b_ld_dags.Rds"))$df
rand_shd_mean_lddags = mean(dfplotDAG[dfplotDAG$method == dfplotDAG$method[1], "shd_rand"])
dfplotDAG$method = sapply(dfplotDAG$method, dict, USE.NAMES = FALSE)
f_ld_dags = ggplot(dfplotDAG, aes(totalSamples, SHD, fill=factor(method))) +
  geom_boxplot() + labs(fill=legend_str, title = "Low-dim. DAGs") + xlab(NULL)+ylab(NULL)+
  geom_hline(yintercept = rand_shd_mean_lddags, linetype="dotted") +  
  theme(plot.title = element_text(face="bold"))

# high-dim polytree
l03 = readRDS(file.path(data_folder, "1b_hd_pt.Rds"))
l03_GIES = readRDS(file.path(data_folder, "1b_hd_pt_GIES.Rds"))
df03_all = rbind(l03$df, l03_GIES$df)
df03_all = df03_all[df03_all$interventionSize == 10, ]
rand_shd_mean2 = mean(df03_all[df03_all$method == df03_all$method[1], "shd_rand"])
df03_all$method = sapply(df03_all$method, dict, USE.NAMES = FALSE)
f_hd_PT = ggplot(df03_all, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  geom_hline(yintercept = rand_shd_mean2, linetype="dotted") + xlab("number of samples")+
  labs(title="High-dim. Polytrees",fill=legend_str) + 
  theme(plot.title = element_text(face="bold")) +
  coord_trans(y="log10")

# high-dim dags
l7 = readRDS(file.path(data_folder,"1b_hd_dags.Rds"))
l7_GIES = readRDS(file.path(data_folder,"1b_hd_dags_GIES.Rds"))
nbh5 = rbind(l7$df, l7_GIES$df)
nbh5plot = nbh5[nbh5$totalSamples != "3000", ]
rand_shd_mean = mean(nbh5plot[nbh5plot$method == nbh5plot$method[1], "shd_rand"])
nbh5plot$method = sapply(nbh5plot$method, dict, USE.NAMES = FALSE)
f_hd_DAG = ggplot(nbh5plot, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  geom_hline(yintercept = rand_shd_mean, linetype="dotted") +  
  labs(title="High-dim. DAGs",fill=legend_str) + xlab("number of samples") + ylab(NULL) +coord_trans(y="log10") +
  theme(plot.title = element_text(face="bold"))
  
# assemble figure
fig1b = f_ld_PT + f_ld_dags + f_hd_PT + f_hd_DAG + 
  plot_layout(nrow=2,ncol = 2, guides="collect")
fig1b



# ---- table 1 (runtimes in high-dim) ----
create_runtime_tbl = function(df,fct){
  methods = unique(df$method)
  ss = unique(df$totalSamples)
  tmp = array(0, length(methods) * length(ss))
  i = 1
  for(s in ss) {
    for(m in methods){
      tmp[i] = fct(df[df$method == m & df$totalSamples == s, "time_s"])
      i = i + 1
    }
  }
  mat = matrix(tmp, length(methods), length(ss))
  rownames(mat) = methods
  colnames(mat) = ss
  return(mat)
}
runtime_pt_mean = create_runtime_tbl(df03_all, mean)
runtime_pt_max = create_runtime_tbl(df03_all, max)
runtime_dags_mean = create_runtime_tbl(nbh5plot, mean)
runtime_dags_max = create_runtime_tbl(nbh5plot, max)

# display results
runtime_pt_mean
runtime_pt_max
runtime_dags_mean
runtime_dags_max

