load_dir = 'cluster_computation/Results_monday/'
library(ggpubr)
library(patchwork)

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



#---- polytrees -----
l01 = readRDS(paste(load_dir, "04-01.Rds",sep=""))
l01_gies = readRDS(paste(load_dir, "04-01_GIES.Rds",sep=""))
l03 = readRDS(paste(load_dir, "04-03.Rds",sep=""))
l03_GIES = readRDS(paste(load_dir, "05-03_GIES.Rds",sep=""))
ndsBIC = readRDS(paste(load_dir, "05-01.Rds",sep=""))
nspBIC = readRDS(paste(load_dir, "05-03.Rds",sep=""))

df01_all = rbind(l01$df, l01_gies$df, ndsBIC$df)
df01_all = df01_all[df01_all$interventionSize == 10,]
df01_all = rbind(df01_all, create_randDag(df01_all))
df03_all = rbind(l03$df, l03_GIES$df, nspBIC$df)
df03_all = df03_all[df03_all$interventionSize == 10, ]
df03_all = rbind(df03_all, create_randDag(df03_all))

(fig2_hd = ggplot(df01_all, aes(ndatasets, SHD, fill=factor(method))) + geom_boxplot() +
    ggplot(df03_all, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))
(fig2_hd = ggplot(df01_all, aes(ndatasets, SID, fill=factor(method))) + geom_boxplot() +
    ggplot(df03_all, aes(totalSamples, SID, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))
(fig2_hd = ggplot(df01_all, aes(ndatasets, time_s, fill=factor(method))) + geom_boxplot() +
    ggplot(df03_all, aes(totalSamples, time_s, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))


# plot zoom in
# time_s
df01plot = df01_all[df01_all$proc != "GIES",]
df01plot = df01plot[df01plot$method != "shd_rand",]
df03plot = df03_all[df03_all$proc != "GIES",]
df03plot = df03plot[df03plot$method != "shd_rand",]
(fig2_hd = ggplot(df01plot, aes(ndatasets, SHD, fill=factor(method))) + geom_boxplot() +
    ggplot(df03plot, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))
(fig2_hd = ggplot(df01plot, aes(ndatasets, SID, fill=factor(method))) + geom_boxplot() +
    ggplot(df03plot, aes(totalSamples, SID, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))
(fig2_hd = ggplot(df01plot, aes(ndatasets, time_s, fill=factor(method))) + geom_boxplot() +
    ggplot(df03plot, aes(totalSamples, time_s, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))




# ---- check out dags -----
l6 = readRDS(paste(load_dir,"05-06.Rds",sep=''))
l6_GIES = readRDS(paste(load_dir,"05-06_GIES.Rds",sep=''))
nbh1 = rbind(l6$df, l6_GIES$df)
nbh1 = rbind(nbh1, create_dagOpt_dagRand(nbh1))
l7 = readRDS(paste(load_dir,"05-07.Rds",sep=''))
l7_GIES = readRDS(paste(load_dir,"05-07_GIES.Rds",sep=''))
nbh5 = rbind(l7$df, l7_GIES$df)
nbh5 = rbind(nbh5, create_dagOpt_dagRand(nbh5))

dfAll = rbind(nbh1, nbh5)
y = "time_s"
ggplot(dfAll, aes(totalSamples, dfAll[,y], fill=method)) + geom_boxplot() + ylab(y)+
  facet_grid(cols = ggplot2::vars(dag_nbh), labeller = labeller(.default = label_both))


# only get simple plots of 1simp, 3simp
df_plot = dfAll[dfAll$proc != "GIES",]
df_plot = df_plot[df_plot$method != "shd_opt",]
df_plot = df_plot[df_plot$method != "shd_rand",]
cs = hue_pal()(5)[c(2,3)]
y = "time_s"
ggplot(df_plot, aes(totalSamples, df_plot[,y], fill=method)) + geom_boxplot() + ylab(y)+
  facet_grid(cols = ggplot2::vars(dag_nbh), labeller = labeller(.default = label_both)) + 
  scale_fill_manual(values=cs)

# only get nbh=1
df_plot = dfAll[dfAll$dag_nbh == 1, ]
y = "SID"
ggplot(df_plot, aes(totalSamples, df_plot[,y], fill=method)) + geom_boxplot() + ylab(y)+
  facet_grid(cols = ggplot2::vars(dag_nbh), labeller = labeller(.default = label_both))








