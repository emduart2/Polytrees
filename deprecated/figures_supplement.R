# computations for old figure 2
# ---- figure 2
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_solo),
  interventionSize = c(intervSize),
  ndatasets = c(nds_range),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("gtruth","1simp"),c("gtruth","1"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("TEST","BIC"))
saveRDS(l, file = paste(save_dir,"02-01.Rds",sep=""))


# varying nsamples
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_range),
  interventionSize = c(intervSize),
  ndatasets = c(nds_solo),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("gtruth","1simp"),c("gtruth","1"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("TEST","BIC"))
saveRDS(l, file = paste(save_dir,"02-02.Rds",sep=""))



# ----- figure 1b (high-dim orientation recovery) ----

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
hd_ort_nds = rbind(
  readRDS(paste(data_folder, "02-01.Rds",sep=""))$df)
hd_ort_nss = rbind(
  readRDS(paste(data_folder, "02-02.Rds",sep=""))$df)
hd_ort_nds$method = sapply(hd_ort_nds$method, dict, USE.NAMES = FALSE)
hd_ort_nss$method = sapply(hd_ort_nss$method, dict, USE.NAMES = FALSE)

hd_ort_nds = hd_ort_nds[(hd_ort_nds$proc %in% c("1","1simp","3","3simp")),]
hd_ort_nss = hd_ort_nss[(hd_ort_nss$proc %in% c("1","1simp","3","3simp")),]

y = "SHD"
ggplot(hd_ort_nds, aes(ndatasets, hd_ort_nds[, y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+ 
  xlab("number of interv. datasets")+ 
  # coord_trans(y="log10") +
  ggplot(hd_ort_nss, aes(totalSamples, hd_ort_nss[,y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+ 
  xlab("number of samples") + 
  # coord_trans(y="log10") +
  plot_layout(ncol = 1, guides="collect")





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