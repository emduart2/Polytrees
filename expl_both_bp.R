library(patchwork)


# get a feeling for the algorithm
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(200),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("mean","1"),c("mean","3"))
l301 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str

f300 = ggplot(l200$df, aes(totalSamples, shd_all,fill=factor(method))) +
  geom_boxplot() + labs(title=l$str)
f301 = ggplot(l200$df, aes(totalSamples, time_all, fill = factor(method)))+
  geom_boxplot()
f302 = ggplot(l200$df, aes(totalSamples, time_ort, fill = factor(method)))+
  geom_boxplot()
f303 = (f300+theme(legend.position="NONE")) + (f301+theme(legend.position="NONE")) +
  f302 + plot_layout(nrow=1)


# compare runtimes
df_params <- expand.grid(
  tsize = c(60),
  totalSamples = c(350),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l302 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str
f302 <- ggplot(l302$df, aes(totalSamples, shd_all,fill=factor(method))) +
  geom_boxplot() + labs(title=l302$str) + 
  ggplot(l302$df, aes(totalSamples, time_all, fill = factor(method)))+
  geom_boxplot() + plot_layout(nrow=1,guides='collect')
f302
