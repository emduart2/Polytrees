library(RColorBrewer)
library(patchwork)

# useful plot functions
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


# ---- figure 2 -----
l1 = readRDS("cluster_computation/01_results/l101.Rds")
l2 = readRDS("cluster_computation/01_results/l102.Rds")

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
   "gtruth,3,BIC" = "P.3, refined, BIC",
   "gtruth,3,TEST" = "P.3, refined, IRC",
   "gtruth,3simp,BIC" = "P.3, simple, BIC",
   "gtruth,3simp,TEST" = "P.3, simple, IRC",
   
   "mean,1,BIC" = "P.1, refined, BIC",
   "mean,1,TEST" = "P.1, refined, IRC",
   "mean,1simp,BIC" = "P.1, simple, BIC",
   "mean,1simp,TEST" = "P.1, simple, IRC",
   "mean,2,BIC" = "P.2, refined, BIC",
   "mean,2,TEST" = "P.2, refined, IRC",
   "mean,2simp,BIC" = "P.2, simple, BIC",
   "mean,2simp,TEST" = "P.2, simple, IRC",
   "mean,3,BIC" = "P.3, refined, BIC",
   "mean,3,TEST" = "P.3, refined, IRC",
   "mean,3simp,BIC" = "P.3, simple, BIC",
   "mean,3simp,TEST" = "P.3, simple, IRC",
   
   "GIES,GIES,TEST" = "GIES"
))}
l1$df$method = sapply(l1$df$method, dict, USE.NAMES = FALSE)
l2$df$method = sapply(l2$df$method, dict, USE.NAMES = FALSE)
legend_str = "Method"

ggplot(l1$df, aes(ndatasets, l1$df[, y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+
  ggplot(l2$df, aes(totalSamples, l2$df[,y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+
  plot_layout(ncol = 1, guides="collect") 


# ---- figure 3 ----
legend_title = "Method"
load_dir = 'cluster_computation/Results_monday/'
# f1: low-dimensional?
# f1 = ggplot(rbind(dfBIC, dfTESTGIES), aes(totalSamples, SHD, fill=method)) + 
#   geom_boxplot() + labs(title = "low-dim Polytrees")

ll_nds = rbind(
  readRDS("local_computations/01-01.Rds")$df,
  readRDS("local_computations/01-03.Rds")$df,
  readRDS("local_computations/01-04.Rds")$df)
ll_ss = rbind(
  readRDS("local_computations/01-02.Rds")$df,
  readRDS("local_computations/01-05.Rds")$df,
  readRDS("local_computations/01-06.Rds")$df)
f1 = ggplot(ll_ss, aes(totalSamples, SHD, fill=method)) + 
  geom_boxplot() + labs(title = "low-dim Polytrees")
f1

# f2: high-dimpolytree
# l01 = readRDS(paste(load_dir, "04-01.Rds",sep=""))
# l01_gies = readRDS(paste(load_dir, "04-01_GIES.Rds",sep=""))
l03 = readRDS(paste(load_dir, "04-03.Rds",sep=""))
l03_GIES = readRDS(paste(load_dir, "05-03_GIES.Rds",sep=""))
ndsBIC = readRDS(paste(load_dir, "05-01.Rds",sep=""))
nspBIC = readRDS(paste(load_dir, "05-03.Rds",sep=""))
df03_all = rbind(l03$df, l03_GIES$df, nspBIC$df)
df03_all = df03_all[df03_all$interventionSize == 10, ]
df03_all = df03_all[df03_all$method %in% c("GIES,GIES,TEST","mean,1simp,TEST","mean,3simp,TEST"),]
# df03_all = rbind(df03_all, create_randDag(df03_all))
rand_shd_mean2 = mean(df03_all[df03_all$method == df03_all$method[1], "shd_rand"])
df03_all$method = sapply(df03_all$method, dict, USE.NAMES = FALSE)
f2 = ggplot(df03_all, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  geom_hline(yintercept = rand_shd_mean2, linetype="dotted") + 
  labs(title="High-dim. Polytrees",fill=legend_title) + coord_trans(y="log10")
f2



#   scale_y_continuous(trans = log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x)))
# f2
# 
# + coord_trans(y="log2")


# f3: high-dim dags
l7 = readRDS(paste(load_dir,"05-07.Rds",sep=''))
l7_GIES = readRDS(paste(load_dir,"05-07_GIES.Rds",sep=''))
nbh5 = rbind(l7$df, l7_GIES$df)
# nbh5 = rbind(nbh5, create_dagOpt_dagRand(nbh5))
nbh5plot = nbh5[nbh5$totalSamples != "3000", ]
opt_shd_mean = mean(nbh5plot[nbh5plot$method == nbh5plot$method[1], "SHD_optimal"])
rand_shd_mean = mean(nbh5plot[nbh5plot$method == nbh5plot$method[1], "shd_rand"])
nbh5plot$method = sapply(nbh5plot$method, dict, USE.NAMES = FALSE)
f3 = ggplot(nbh5plot, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  geom_hline(yintercept = rand_shd_mean, linetype="dotted") +  
  geom_hline(yintercept = opt_shd_mean, linetype="dashed") + 
  labs(title="High-dim. DAGs",fill=legend_title)


# f4 runtime of one of the two
library(scales)
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





# ---- fig2 high-dimensional -----

ll1 = readRDS("cluster_computation/03_results_with500/03-05.Rds")
ll1 = readRDS("cluster_computation/04_results/04-05.Rds")
dfBIC = ll1$df[ll1$df$method %in% c("mean,1simp,BIC","mean,3simp,BIC"),]
dfTESTGIES = readRDS("local_computations/01-02.Rds")$df
f1 = ggplot(rbind(dfBIC, dfTESTGIES), aes(totalSamples, SHD, fill=method)) +
  geom_boxplot() + labs(title = "low-dim Polytrees")


ndsBIC = readRDS("cluster_computation/03_results_with500/03-04.Rds")
nssBIC = readRDS("cluster_computation/03_results_with500/03-05.Rds")

