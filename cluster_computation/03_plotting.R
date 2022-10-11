# load_dir = "cluster_computation/03_results_with500/"
load_dir = "cluster_computation/Results_monday/"
experiment_id = "04"
library(ggpubr)

#---- fig 2 but high-dim ----
# - pick one of interventionsize (whatever look sbetter)
# - analyse the trends of varying totalSamples & ndatasets

# varying ndatases with GIES
# - check runtime of GIES 
# - if it terminates, we want to say that we are faster

l01 = readRDS(paste(load_dir, experiment_id,"-01.Rds",sep=""))
l01_gies = readRDS(paste(load_dir, experiment_id,"-01_GIES.Rds",sep=""))
l03 = readRDS(paste(load_dir, experiment_id,"-03.Rds",sep=""))
l03_GIES = readRDS(paste(load_dir, "05-03_GIES.Rds",sep=""))
ndsBIC = readRDS(paste(load_dir, "05-01.Rds",sep=""))
nspBIC = readRDS(paste(load_dir, "05-03.Rds",sep=""))

df01_all = rbind(l01$df, l01_gies$df, ndsBIC$df)
df01_all = df01_all[df01_all$interventionSize == 10,]
df03_all = rbind(l03$df, dfGIES, nspBIC$df)
df03_all = df03_all[df03_all$interventionSize == 10, ]

(fig2_hd = ggplot(df01_all, aes(ndatasets, SHD, fill=factor(method))) + geom_boxplot() +
  ggplot(df03_all, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))
(fig2_hd = ggplot(df01_all, aes(ndatasets, SID, fill=factor(method))) + geom_boxplot() +
    ggplot(df03_all, aes(totalSamples, SID, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))
(fig2_hd = ggplot(df01_all, aes(ndatasets, time_s, fill=factor(method))) + geom_boxplot() +
    ggplot(df03_all, aes(totalSamples, time_s, fill=factor(method))) + geom_boxplot() +
    plot_layout(nrow=1, guides="collect") + plot_annotation(title=l01$str))



#---- fig 2 but with all intervention types ----
# - we want to see if there is any difference with different intervention types
l04 = readRDS(paste(load_dir, experiment_id,"-04.Rds",sep=""))
l05 = readRDS(paste(load_dir, experiment_id,"-05.Rds",sep=""))

y = "SHD"
dfAll = rbind(l04$df, l05$df)
y_limit = c(min(dfAll[,y]), max(dfAll[,y]))
df04per = l04$df[l04$df$kindOfIntervention == "perfect", ]
df05per = l05$df[l05$df$kindOfIntervention == "perfect", ]
f_perfect = ggplot(df04per, aes(ndatasets, df04per[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  ggplot(df05per, aes(totalSamples, df05per[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  plot_layout(nrow=1, guides="collect") + plot_annotation(title=paste("perfect, ",l04$str))
df04imp = l04$df[l04$df$kindOfIntervention == "imperfect", ]
df05imp = l05$df[l05$df$kindOfIntervention == "imperfect", ]
f_imperfect = ggplot(df04imp, aes(ndatasets, df04imp[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  ggplot(df05imp, aes(totalSamples, df05imp[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  plot_layout(nrow=1, guides="collect") + plot_annotation(title=paste("imperfect, ",l04$str))
df04inh = l04$df[l04$df$kindOfIntervention == "inhibitory", ]
df05inh = l05$df[l05$df$kindOfIntervention == "inhibitory", ]
f_inhibitory = ggplot(df04inh, aes(ndatasets, df04inh[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  ggplot(df05inh, aes(totalSamples, df05inh[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  plot_layout(nrow=1, guides="collect") + plot_annotation(title=paste("inhibitory, ",l04$str))
df04rand = l04$df[l04$df$kindOfIntervention == "random", ]
df05rand = l05$df[l05$df$kindOfIntervention == "random", ]
f_random = ggplot(df04rand, aes(ndatasets, df04rand[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  ggplot(df05rand, aes(totalSamples, df05rand[,y], fill=factor(method))) + geom_boxplot() +
  ylab(y) + ylim(y_limit) +
  plot_layout(nrow=1, guides="collect") + plot_annotation(title=paste("random, ",l04$str))
# f_perfect
# f_imperfect
# f_inhibitory
# f_random
f_perfect
ggarrange(f_perfect, f_imperfect, f_inhibitory, f_random, nrow=2, ncol=2)


# #---- potential fig 3: high-dim DAG setting for varying nsamples ----
# - pick one interventionsize (whatever looks better)
# - compare runtime of algorithms, ours should be faster than GIES
l06 = readRDS(paste(load_dir, experiment_id,"-06.Rds",sep=""))
l06_GIES = readRDS(paste(load_dir, experiment_id,"-06_GIES.Rds",sep=""))
l07 = readRDS(paste(load_dir, experiment_id,"-07.Rds",sep=""))
l07_GIES = readRDS(paste(load_dir, experiment_id,"-07_GIES.Rds",sep=""))
l08 = readRDS(paste(load_dir, experiment_id,"-08.Rds",sep=""))
l08_GIES = readRDS(paste(load_dir, experiment_id,"-08_GIES.Rds",sep=""))

intervSize = 10 # 2 or 10

nbh1_all = rbind(l06_GIES$df, l06$df)
nbh5_all = rbind(l07_GIES$df, l07$df)
nbh10_all = rbind(l08_GIES$df, l08$df)
dfAll = rbind(nbh1_all, nbh5_all, nbh10_all)
dfAll = dfAll[dfAll$interventionSize == intervSize,]
dfAll$dag_nbh = factor(dfAll$dag_nbh)
df_optSHD = dfAll[dfAll$method== dfAll$method[1],] # select all entries from one method to get all simulated graphs exactly once
df_optSHD$method = "shd_opt"
df_optSHD$SHD = df_optSHD$SHD_optimal
(f_shd = ggplot(rbind(dfAll, df_optSHD), aes(dag_nbh, SHD, fill=factor(method))) + geom_boxplot() +
  labs(title = paste(l06$str,", intervSize: ",intervSize)))
(f_sid = ggplot(dfAll, aes(dag_nbh, SID, fill=factor(method))) + geom_boxplot() +
  labs(title = paste(l06$str,", intervSize: ",intervSize)))
(f_time = ggplot(dfAll, aes(dag_nbh, time, fill=factor(method))) + geom_boxplot() +
  labs(title = paste(l06$str,", intervSize: ",intervSize)))

f_shd +labs(title="")+ f_sid+labs(title="") + f_time+labs(title="") + plot_layout(nrow=1, guides="collect") +
  plot_annotation(title = paste(l06$str,", intervSize: ",intervSize))


