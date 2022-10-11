library(ggplot2)
library(patchwork)

load_dir = "cluster_computation/Results_monday/"

# plot fig2 in high-dim with GIES and BIC
l01 = readRDS(paste(load_dir, "04-01.Rds",sep=""))
l01_gies = readRDS(paste(load_dir, "04-01_GIES.Rds",sep=""))
l03 = readRDS(paste(load_dir, "04-03.Rds",sep=""))
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





# compare GIES vs our algos on DAGs
nbh1 = readRDS(paste(load_dir, "05-06.Rds",sep=""))
nbh5 = readRDS(paste(load_dir, "05-07.Rds",sep=""))
nbh10 = readRDS(paste(load_dir, "05-08.Rds",sep=""))
dfAll = rbind(nbh1$df, nbh5$df, nbh10$df)

ggplot(dfAll, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() + 
  facet_grid(cols = ggplot2::vars(dag_nbh), labeller = labeller(.default=label_both)) +
  labs(main=nbh1$str)


