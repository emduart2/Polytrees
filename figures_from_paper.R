library(RColorBrewer)
library(patchwork)
library(scales)
library(latex2exp)
librarhy(ggplot2)

data_folder = "data/"  # with final '/'
figure_width = 1000

#---- useful plot functions ----
# function to add optimal and random shd to the methods 
create_dagOpt_dagRand = function(dfAll){
  df_optSHD = dfAll[dfAll$method == dfAll$method[1],] # select all entries from one method to get all simulated graphs exactly once
  df_optSHD$method = "shd_opt"
  df_optSHD$SHD = df_optSHD$SHD_optimal
  
  df_randSHD = dfAll[dfAll$method == dfAll$method[1],] # select all entries from one method to get all simulated graphs exactly once
  df_randSHD$method = "shd_rand"
  df_randSHD$SHD = df_randSHD$shd_rand
  
  return(rbind(df_optSHD, df_randSHD))
}

create_randDag = function(dfAll){
  df_randSHD = dfAll[dfAll$method == dfAll$method[1],] # select all entries from one method to get all simulated graphs exactly once
  df_randSHD$method = "shd_rand"
  df_randSHD$SHD = df_randSHD$shd_rand
  return(df_randSHD)
}

# function to compare the unique values of two dfs
compare_dfs = function(df1,df2){
  for(c in names(df1)){
    print(c)
    print(unique(df1[,c]))
    print(unique(df2[,c]))
    print("-----------------------------")
  }
}


# colorcode of figures
ind = c(3,4,7,8)
reds = brewer.pal(n = 9, name = 'Reds')[c(3,4,7,8)]
blues = brewer.pal(n = 9, name = 'Greens')[c(3,4,7,8)]
greens = brewer.pal(n = 9, name = 'Blues')[c(3,4,7,8)]
pattern = c(4,3,2,1)
colors1 = c(reds[pattern], blues[pattern], greens[pattern])
dict = function(x){return(switch(x,
   "gtruth,1,BIC" = "P.1, refined, BIC",
   "gtruth,1,TEST" = "P.1, refined, IRC",
   "gtruth,1simp,BIC" = "P.1, simple, BIC",
   "gtruth,1simp,TEST" = "P.1, simple, IRC",
   "gtruth,2,BIC" = "P.2, refined, BIC",
   "gtruth,2,TEST" = "P.2, refined, IRC",
   "gtruth,2simp,BIC" = "P.2, simple, BIC",
   "gtruth,2simp,TEST" = "P.2, simple, IRC",
   "gtruth,3,BIC" = "P.2, refined, BIC",
   "gtruth,3,TEST" = "P.2, refined, IRC",
   "gtruth,3simp,BIC" = "P.2, simple, BIC",
   "gtruth,3simp,TEST" = "P.1, simple, IRC",
   
   "mean,1,BIC" = "P.1, refined, BIC",
   "mean,1,TEST" = "P.1, refined, IRC",
   "mean,1simp,BIC" = "P.1, simple, BIC",
   "mean,1simp,TEST" = "P.1, simple, IRC",
   "mean,2,BIC" = "P.2, refined, BIC",
   "mean,2,TEST" = "P.2, refined, IRC",
   "mean,2simp,BIC" = "P.2, simple, BIC",
   "mean,2simp,TEST" = "P.2, simple, IRC",
   "mean,3,BIC" = "P.2, refined, BIC",
   "mean,3,TEST" = "P.2, refined, IRC",
   "mean,3simp,BIC" = "P.2, simple, BIC",
   "mean,3simp,TEST" = "P.2, simple, IRC",
   
   "GIES,GIES,TEST" = "GIES"
))}
legend_str = "Method"




# ---- figure 1 (high-dim skeleton recovery) ----
hd_skel1 = readRDS(paste(data_folder, "02-03.Rds",sep=""))
hd_skel1 = hd_skel1$df[hd_skel1$df$method %in% c("mean","median","ltest","baseObs"), ]
f1 = ggplot(hd_skel1, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.") + xlab("number of samples")
hd_skel2 = readRDS(paste(data_folder, "02-04.Rds",sep=""))
hd_skel2 = hd_skel2$df[hd_skel2$df$method %in% c("mean","median","ltest","baseObs"), ]
# hd_skel2$ndatasets = factor(as.numeric(as.character(hd_skel2$ndatasets))-1)
f2 = ggplot(hd_skel2, aes(ndatasets, SHD, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.") + xlab("number of interv. datasets")
hd_skel3 = readRDS(paste(data_folder, "02-05.Rds",sep=""))
hd_skel3 = hd_skel3$df[hd_skel3$df$method %in% c("mean","median","ltest","baseObs"), ]
f3 = ggplot(hd_skel3, aes(interventionSize, SHD, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.") + xlab("size of interventions")
fig1 = f1 + f2 + f3 + plot_layout(nrow = 1, guides="collect")
fig1



# ----- figure 2 (high-dim orientation recovery) ----
hd_ort_nds = rbind(
  readRDS(paste(data_folder, "02-01.Rds",sep=""))$df)
hd_ort_nss = rbind(
  readRDS(paste(data_folder, "02-02.Rds",sep=""))$df)
hd_ort_nds$method = sapply(hd_ort_nds$method, dict, USE.NAMES = FALSE)
hd_ort_nss$method = sapply(hd_ort_nss$method, dict, USE.NAMES = FALSE)

hd_ort_nds = hd_ort_nds[(hd_ort_nds$proc %in% c("1","1simp","3","3simp")),]
hd_ort_nss = hd_ort_nss[(hd_ort_nss$proc %in% c("1","1simp","3","3simp")),]

y = "time_s"
ggplot(hd_ort_nds, aes(ndatasets, hd_ort_nds[, y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+ 
  xlab("number of interv. datasets")+ 
  # coord_trans(y="log10") +
  ggplot(hd_ort_nss, aes(totalSamples, hd_ort_nss[,y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+ 
  xlab("number of samples") + 
  # coord_trans(y="log10") +
  plot_layout(ncol = 1, guides="collect")






# ---- figure 3 ----
legend_str = "Method"
# ts_str = TeX(r"(n)")
ts_str = "number of samples"

# low-dim polytree
dfplotPT = readRDS(paste(data_folder, "lowdim_polytrees.Rds",sep=""))$df
rand_shd_mean_ldpt = mean(dfplotPT[dfplotPT$method == dfplotPT$method[1], "shd_rand"])
dfplotPT$method = sapply(dfplotPT$method, dict, USE.NAMES = FALSE)
f_ld_PT = ggplot(dfplotPT, aes(totalSamples, SHD, fill=factor(method))) +
  geom_boxplot() + labs(fill=legend_str, title = "Low-dim. Polytrees")+ xlab(NULL)+
  geom_hline(yintercept = rand_shd_mean_ldpt, linetype="dotted") +
  theme(plot.title = element_text(face="bold"))
f_ld_PT

# low-dim DAGs
dfplotDAG = readRDS(paste(data_folder, "lowdim_dags.Rds",sep=""))$df
opt_shd_mean_lddags = mean(dfplotDAG[dfplotDAG$method == dfplotDAG$method[1], "SHD_optimal"])
rand_shd_mean_lddags = mean(dfplotDAG[dfplotDAG$method == dfplotDAG$method[1], "shd_rand"])
dfplotDAG$method = sapply(dfplotDAG$method, dict, USE.NAMES = FALSE)
f_ld_dags = ggplot(dfplotDAG, aes(totalSamples, SHD, fill=factor(method))) +
  geom_boxplot() + labs(fill=legend_str, title = "Low-dim. DAGs") + xlab(NULL)+ylab(NULL)+
  geom_hline(yintercept = rand_shd_mean_lddags, linetype="dotted") +  
  # geom_hline(yintercept = opt_shd_mean_lddags, linetype="dashed") +
  theme(plot.title = element_text(face="bold"))
  


# high-dim polytree
# l01 = readRDS(paste(load_dir, "04-01.Rds",sep=""))
# l01_gies = readRDS(paste(load_dir, "04-01_GIES.Rds",sep=""))
l03 = readRDS(paste(load_dir, "04-03.Rds",sep=""))
l03_GIES = readRDS(paste(load_dir, "05-03_GIES.Rds",sep=""))
ndsBIC = readRDS(paste(load_dir, "05-01.Rds",sep=""))
nspBIC = readRDS(paste(load_dir, "05-03.Rds",sep=""))
df03_all = rbind(l03$df, l03_GIES$df, nspBIC$df)
df03_all = df03_all[df03_all$interventionSize == 10, ]
df03_all = df03_all[df03_all$method %in% c("GIES,GIES,TEST","mean,1simp,TEST","mean,3simp,TEST"),]
rand_shd_mean2 = mean(df03_all[df03_all$method == df03_all$method[1], "shd_rand"])
df03_all$method = sapply(df03_all$method, dict, USE.NAMES = FALSE)
f_hd_PT = ggplot(df03_all, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  geom_hline(yintercept = rand_shd_mean2, linetype="dotted") + xlab(ts_str)+
  labs(title="High-dim. Polytrees",fill=legend_str) + 
  theme(plot.title = element_text(face="bold")) +
  coord_trans(y="log10")



# high-dim dags
l7 = readRDS(paste(load_dir,"05-07.Rds",sep=''))
l7_GIES = readRDS(paste(load_dir,"05-07_GIES.Rds",sep=''))
nbh5 = rbind(l7$df, l7_GIES$df)
# nbh5 = rbind(nbh5, create_dagOpt_dagRand(nbh5))
nbh5plot = nbh5[nbh5$totalSamples != "3000", ]
opt_shd_mean = mean(nbh5plot[nbh5plot$method == nbh5plot$method[1], "SHD_optimal"])
rand_shd_mean = mean(nbh5plot[nbh5plot$method == nbh5plot$method[1], "shd_rand"])
nbh5plot$method = sapply(nbh5plot$method, dict, USE.NAMES = FALSE)
f_hd_DAG = ggplot(nbh5plot, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  geom_hline(yintercept = rand_shd_mean, linetype="dotted") +  
  # geom_hline(yintercept = opt_shd_mean, linetype="dashed") + 
  labs(title="High-dim. DAGs",fill=legend_str) + xlab(ts_str) + ylab(NULL) +coord_trans(y="log10") +
  theme(plot.title = element_text(face="bold"))
  


# assemble figure
fig3 = f_ld_PT + f_ld_dags + f_hd_PT + f_hd_DAG + 
  plot_layout(nrow=2,ncol = 2, guides="collect")
fig3
# size: 900 by 450




# ---- table 1 (runtime in high-dim DAG setting) ----
df = nbh5plot
methods = unique(df$method)
ss = unique(df$totalSamples)
tmp = array(0, length(methods) * length(ss))
i = 1
for(s in ss) {
  for(m in methods){
    tmp[i] = mean(df[df$method == m & df$totalSamples == s, "time_s"])
    i = i + 1
  }
}
mat = matrix(tmp, length(methods), length(ss))
rownames(mat) = methods
colnames(mat) = ss
mat




# ---- orientation high-dimensional ??? -----
# ll1 = readRDS("cluster_computation/03_results_with500/03-05.Rds")
ll1 = readRDS("cluster_computation/04_results/04-05.Rds")$df
df_ssBIC = ll1[ll1$kindOfIntervention == "perfect",]

ll2 = readRDS("cluster_computation/04_results/04-06.Rds")$df

tmp = readRDS("local_computations/01-02.Rds")$df # lowdimensional
tmp$time = tmp$time_s
tmp$time_s = NULL



f1 = ggplot(rbind(dfBIC, dfTESTGIES), aes(totalSamples, SHD, fill=method)) +
  geom_boxplot() + labs(title = "low-dim Polytrees")

ndsBIC = readRDS("cluster_computation/03_results_with500/03-04.Rds")
nssBIC = readRDS("cluster_computation/03_results_with500/03-05.Rds")




nds vary: 1,2,3,3simp * TEST
















# ---- SUPPLEMENTARY MATERIAL PLOTS -----

# low-dimensional polytrees (all stuff)
load_dir = 'cluster_computation/Results_monday/'
# ll_nds = rbind(
#   readRDS("local_computations/01-01.Rds")$df,
#   readRDS("local_computations/01-03.Rds")$df,
#   readRDS("local_computations/01-04.Rds")$df)
ll_ss = rbind(
  readRDS("local_computations/01-02.Rds")$df,
  readRDS("local_computations/01-05.Rds")$df,
  readRDS("local_computations/01-06.Rds")$df)
ll_plot = ll_ss[ll_ss$method %in% c("GIES,GIES,TEST","mean,1simp,TEST","mean,3simp,TEST"),]
f1 = ggplot(ll_plot, aes(totalSamples, SHD, fill=method)) + 
  geom_boxplot() + labs(title = "low-dim Polytrees")


# ---- low-dim orientations (all) -----
l1 = readRDS("cluster_computation/01_results/l101.Rds")
l2 = readRDS("cluster_computation/01_results/l102.Rds")

l1$df$method = sapply(l1$df$method, dict, USE.NAMES = FALSE)
l2$df$method = sapply(l2$df$method, dict, USE.NAMES = FALSE)
ggplot(l1$df, aes(ndatasets, l1$df[, y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+
  ggplot(l2$df, aes(totalSamples, l2$df[,y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+
  plot_layout(ncol = 1, guides="collect") 


# ---- high-dim dags runtime ----
dfplot = nbh5[nbh5$method != "shd_opt",]
dfplot = dfplot[dfplot$method != "shd_rand",]
dfplot = dfplot[dfplot$totalSamples != "3000",]
dfplot$method = sapply(dfplot$method, dict, USE.NAMES = FALSE)
f4 = ggplot(dfplot, aes(totalSamples, time_s, fill=factor(method))) + geom_boxplot() +
  labs(title="High-dim. DAGs",fill=legend_title) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))
f4

# merge figures (first row polytrees, second row DAGs)
# fig3_1 = f1 + f2 + f3 + f4 + plot_layout(nrow = 2, ncol=2, guides="collect") 
(fig3_onlyHD = f2 + f3 + f4 + plot_layout(nrow = 1, guides="collect") ) 
# size: 1250 450



# high-dimensional polytrees from data from cluster computations