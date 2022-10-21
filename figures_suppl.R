library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

save_dir = "/home/lenny/Documents/Mathematische Statistik Lehrstuhl/fig_suppl"

#---- helpful plot functions ----
# dictionary for converting axis labels
dict_axisLabs = function(x){
  val = switch(x,
               "ndatasets" = "d",
               "totalSamples" = "n",
               "interventionSize" = "k",
               "time_s" = "time (s)",
               "dag_nbh" = "e"
  )
  if(is.null(val)){
    return(x)
  }
  return(val)
}



# function to make method names pretty
pretty_names = function(df){
  dict_method = function(x){return(switch(x,
     "gtruth,1,BIC" = "P.1, refined, BIC",
     "gtruth,1,TEST" = "P.1, refined, IRC",
     "gtruth,1simp,BIC" = "P.1, simple, BIC",
     "gtruth,1simp,TEST" = "P.1, simple, IRC",
     "gtruth,3,BIC" = "P.2, refined, BIC",
     "gtruth,3,TEST" = "P.2, refined, IRC",
     "gtruth,3simp,BIC" = "P.2, simple, BIC",
     "gtruth,3simp,TEST" = "P.2, simple, IRC",
     "gtruth,2,BIC" = "alt., refined, BIC",
     "gtruth,2,TEST" = "alt., refined, IRC",
     "gtruth,2simp,BIC" = "alt., simple, BIC",
     "gtruth,2simp,TEST" = "alt., simple, IRC",
     
     "gtruth,1,BIC" = "P.1, refined, BIC",
     "gtruth,1,IRC" = "P.1, refined, IRC",
     "gtruth,1simp,BIC" = "P.1, simple, BIC",
     "gtruth,1simp,IRC" = "P.1, simple, IRC",
     "gtruth,3,BIC" = "P.2, refined, BIC",
     "gtruth,3,IRC" = "P.2, refined, IRC",
     "gtruth,3simp,BIC" = "P.2, simple, BIC",
     "gtruth,3simp,IRC" = "P.2, simple, IRC",
     "gtruth,2,BIC" = "alt., refined, BIC",
     "gtruth,2,IRC" = "alt., refined, IRC",
     "gtruth,2simp,BIC" = "alt., simple, BIC",
     "gtruth,2simp,IRC" = "alt., simple, IRC",
     
     "mean,1,BIC" = "P.1, refined, BIC",
     "mean,1,TEST" = "P.1, refined, IRC",
     "mean,1simp,BIC" = "P.1, simple, BIC",
     "mean,1simp,TEST" = "P.1, simple, IRC",
     "mean,3,BIC" = "P.2, refined, BIC",
     "mean,3,TEST" = "P.2, refined, IRC",
     "mean,3simp,BIC" = "P.2, simple, BIC",
     "mean,3simp,TEST" = "P.2, simple, IRC",
     "mean,2,BIC" = "alt., refined, BIC",
     "mean,2,TEST" = "alt., refined, IRC",
     "mean,2simp,BIC" = "alt., simple, BIC",
     "mean,2simp,TEST" = "alt., simple, IRC",
     
     "mean,1,BIC" = "P.1, refined, BIC",
     "mean,1,IRC" = "P.1, refined, IRC",
     "mean,1simp,BIC" = "P.1, simple, BIC",
     "mean,1simp,IRC" = "P.1, simple, IRC",
     "mean,3,BIC" = "P.2, refined, BIC",
     "mean,3,IRC" = "P.2, refined, IRC",
     "mean,3simp,BIC" = "P.2, simple, BIC",
     "mean,3simp,IRC" = "P.2, simple, IRC",
     "mean,2,BIC" = "alt., refined, BIC",
     "mean,2,IRC" = "alt., refined, IRC",
     "mean,2simp,BIC" = "alt., simple, BIC",
     "mean,2simp,IRC" = "alt., simple, IRC",
     
     "GIES,GIES,TEST" = "GIES",
     "GIES,GIES,BIC" = "GIES",
     "GIES,GIES,IRC" = "GIES",
     
     "shd_opt" = "SHD opt.",
     "shd_rand" = "SHD rand.",
     
     "baseObs" = "baseObs",
     "ltest" = "ltest",
     "mean" = "mean",
     "median" = "median",
     "baseAll" = "pooled"
  ))}
  df$method = sapply(df$method, dict_method, USE.NAMES = FALSE)
  return(df)
}

# function to plot orientation plots with different interventions
plot_orientations = function(l, x, y){
  
  # remove procedure 2 and extract mean of random
  df = l$df
  df = df[ !(df$proc %in% c("2","2simp")), ]
  rand_SHD = mean(df[df$proc == "random", "SHD"])
  df = df[ df$proc != "random", ]
  df = pretty_names(df)
  
  # setup colors
  reds = brewer.pal(n = 9, name = 'Reds')[c(3,4,7,8)]
  greens = brewer.pal(n = 9, name = 'Greens')[c(3,4,7,8)]
  blues = brewer.pal(n = 9, name = 'Blues')[c(3,4,7,8)]
  pattern = c(4,3,2,1)
  colors1 = c(greens[pattern], blues[pattern]) # colormap for orientations
  
  # plot
  fig = ggplot(df, aes(df[,x], df[,y],fill=factor(method))) +
    geom_boxplot() + xlab(x) + ylab(y) + labs(fill="Method",title=l$str) +
    scale_fill_manual(values=colors1) +
    geom_hline(yintercept = rand_SHD, linetype="dotted") +
    facet_grid(cols = ggplot2::vars(kindOfIntervention),
               labeller = labeller(.default = label_both))
  return(fig)
}

# function for plotting skeleton in a nice way
plot_skeleton = function(df, x, y){
  
  
  
  
  fig = ggplot(df, aes(df[,x], df[,y],fill=factor(method))) +
    geom_boxplot() + xlab(dict_axisLabs(x)) + ylab(dict_axisLabs(y)) + 
    labs(fill="Aggr. fct.") +
    facet_grid(cols = ggplot2::vars(kindOfIntervention),
               labeller = labeller(.default = label_both))
  return(fig)
}


# function that plots single orientation plot in pretty way
pplot_ort = function(df, x, y){
  
  # remove procedure 2 and extract mean of random
  df = df[ !(df$proc %in% c("2","2simp")), ]
  rand_SHD = mean(df[df$proc == "random", "SHD"])
  df = df[ df$proc != "random", ]
  df = pretty_names(df)
  
  # setup colors
  ind = c(2,3,6,7)
  reds = brewer.pal(n = 9, name = 'Reds')[ind]
  greens = brewer.pal(n = 9, name = 'Greens')[ind]
  blues = brewer.pal(n = 9, name = 'Blues')[ind]
  pattern = c(4,3,2,1)
  colors1 = c(greens[pattern], blues[pattern]) # colormap for orientations
  
  # return
  fig = ggplot(df, aes(df[,x], df[,y],fill=factor(method))) +
    geom_boxplot() + xlab(dict_axisLabs(x)) + ylab(dict_axisLabs(y)) + labs(fill="Method") +
    scale_fill_manual(values=colors1)
  if(y != "time_s"){
    fig = fig + geom_hline(yintercept = rand_SHD, linetype="dotted")
  } 
  return(fig)
}

# function to extract parameter value range
unq = function(df){
  cols = c("tsize","totalSamples","interventionSize","ndatasets","kindOfIntervention","use_dags","dag_nbh","method")
  for(col in cols){
    print(paste(col,": ", paste(unique(df[,col]), collapse=",")))
  }
}

# color scheme for skeleton recovery
colors_skel = c(
  "baseObs" = hue_pal()(4)[1],
  "ltest" = hue_pal()(4)[2],
  "mean" = hue_pal()(4)[3],
  "median" = hue_pal()(4)[4],
  "pooled" = "orange"
)
greens = brewer.pal(n = 9, name = 'Greens')[c(8,7,4,3)]
blues = brewer.pal(n = 9, name = 'Blues')[c(8,7,4,3)]
colors_ort = c(
  "P.1, refined, BIC" = greens[1],
  "P.1, refined, IRC" = greens[2],
  "P.1, simple, BIC" = greens [3],
  "P.1, simple, IRC" = greens[4],
  "P.2, refined, BIC" = blues[1],
  "P.2, refined, IRC" = blues[2],
  "P.2, simple, BIC" = blues[3],
  "P.2, simple, IRC" = blues[4],
  "GIES" = hue_pal()(3)[1],
  "SHD rand." = "orange",
  "SHD opt." = "yellow"
)





#---- ort: main analysis -----
l_nds = readRDS("data_new/Orientation/03-01.Rds") #ndatasets
l_smpls = readRDS("data_new/Orientation/03-02.Rds"); # totalsamples
l_intvS = readRDS("data_new/Orientation/03-03.Rds");  # intervSize
fig3 = 
  pplot_ort(ld1$df[ld1$df$kindOfIntervention == "perfect",], "totalSamples","SHD") + 
    labs(title="(a) p=20, d=11, k=2") +  # low-dim base setting
  pplot_ort(l_smpls$df[l_smpls$df$kindOfIntervention == "perfect" & l_smpls$df$totalSamples != 2000,], "totalSamples","SHD") + 
    labs(title="(b) p=500, d=21, k=10") +  # high-dim base setting
  
  pplot_ort(ld2$df[ld2$df$kindOfIntervention == "perfect",], "interventionSize","SHD") + 
  labs(title="(c) p=20, n=500, d=11") +
  pplot_ort(l_intvS$df[l_intvS$df$kindOfIntervention == "perfect" & (!(l_intvS$df$interventionSize %in% c(2,4))),], "interventionSize","SHD") + 
  labs(title="(d) p=500, n=1000, d=21") +
  
  pplot_ort(ld3$df[ld3$df$kindOfIntervention == "perfect",], "ndatasets","SHD") + 
  labs(title="(e) p=20, n=500, k=2") +
  pplot_ort(l_nds$df[l_nds$df$kindOfIntervention == "perfect",], "ndatasets","SHD") +
  labs(title="(f) p=500, n=1000, k=10") +
  
  pplot_ort(ld2$df[ld2$df$kindOfIntervention == "perfect",], "interventionSize","time_s") + 
  labs(title="(g) p=20, n=500, d=11") +
  pplot_ort(l_intvS$df[l_intvS$df$kindOfIntervention == "perfect" & (!(l_intvS$df$interventionSize %in% c(2,4))),], "interventionSize","time_s") + 
  labs(title="(h) p=500, n=1000, d=21") +
  
  plot_layout(ncol=2, guides = "collect") & theme(legend.position = "bottom", plot.title = element_text(face="bold"))
fig3




# ---- dags: different neighborhoods ----
# load_dir = 'old_cluster_comp/'
# l6 = readRDS(paste(load_dir,"05-06.Rds",sep=''))
# l6_GIES = readRDS(paste(load_dir,"05-06_GIES.Rds",sep=''))
# nbh1 = rbind(l6$df, l6_GIES$df)
# nbh1 = rbind(nbh1)
# l7 = readRDS(paste(load_dir,"05-07.Rds",sep=''))
# l7_GIES = readRDS(paste(load_dir,"05-07_GIES.Rds",sep=''))
# nbh5 = rbind(l7$df, l7_GIES$df)
# nbh5 = rbind(nbh5)
# l8 = readRDS(paste(load_dir,"05-08.Rds",sep=''))
# nbh10 = rbind(l8$df)
# 
# 
# dfAll = rbind(nbh1, nbh5, nbh10)
# y = "SHD"
# means = c("1" = 1000, "5" =2000, "10"=3000)
# ggplot(dfAll, aes(totalSamples, dfAll[,y], fill=method)) + geom_boxplot() + ylab(y)+
#   geom_hline(aes(yintercept = value), data=means) +  
#   facet_grid(cols = ggplot2::vars(dag_nbh), labeller = labeller(.default = label_both))
# 
# ylims = c(min(dfAll$SHD), max(dfAll$SHD))
# nbh1_rand = mean(nbh1[nbh1$method == "shd_rand","SHD"])
# nbh1_plt = nbh1[ !(nbh1$method %in% c("shd_opt","shd_rand")),]
# # nbh1_plt = pretty_names(nbh1_plt)
# fnbh1 = ggplot(nbh1_plt, aes(totalSamples, nbh1_plt[,y], fill=factor(method))) + geom_boxplot()+
#   xlab(dict_axisLabs("totalSamples")) + ylab(dict_axisLabs(y)) + ylim(ylims)+
#   geom_hline(yintercept = nbh1_rand, linetype = "dotted") + 
#   labs(title = "p=500, d=21, k=10, h=1")
# nbh5_rand = mean(nbh5[nbh5$method == "shd_rand","SHD"])
# nbh5_plt = nbh5[ !(nbh5$method %in% c("shd_opt","shd_rand")),]
# # nbh5_plt = pretty_names(nbh5_plt)
# fnbh5 = ggplot(nbh5_plt, aes(totalSamples, nbh5_plt[,y], fill=factor(method))) + geom_boxplot()+
#   xlab(dict_axisLabs("totalSamples")) + ylab(dict_axisLabs(y)) +ylim(ylims)+
#   geom_hline(yintercept = nbh5_rand, linetype = "dotted") + 
#   labs(title = "p=500, d=21, k=10, h=5")
# nbh10_rand = mean(nbh10[nbh10$method == "shd_rand","SHD"])
# nbh10_plt = nbh10[ !(nbh10$method %in% c("shd_opt","shd_rand")),]
# # nbh10_plt = pretty_names(nbh10_plt)
# # nbh10_plt$method = factor(nbh10_plt$method, levels = c("GIES","P.1,simple,IRC","P.2,simple,IRC"))
# fnbh10 = ggplot(nbh10_plt, aes(totalSamples, nbh10_plt[,y], fill=factor(method))) + geom_boxplot()+
#   xlab(dict_axisLabs("totalSamples")) + ylab(dict_axisLabs(y)) +ylim(ylims)+
#   geom_hline(yintercept = nbh10_rand, linetype = "dotted") + 
#   labs(title = "p=500, d=21, k=10, h=10")
# fnbhs = fnbh1 + fnbh5 + fnbh10 + plot_layout(nrow=1, guides="collect") &
#   theme(plot.title = element_text(face="bold"))
# fnbhs
# 
# 
# 
# # all mean,3simp,TEST
# l = readRDS("data_new/Results/05-vlarge_011.Rds"); unq(l$df) # PT: 1000 nodes, totalSamples
# l = readRDS("data_new/Results/05-vlarge_012.Rds"); unq(l$df) # DAG: 1000 nodes, nbh=1,5,10,20
#   # vlarge for GIES did not finish
#   #  ATM: we are running vlarge with GIES
# 
# l = readRDS("data_new/Results/05-large_011.Rds"); unq(l$df) # DAG, 1000 nodes, 20 nbh
# # l = readRDS("data_new/Results/05-large_022.Rds"); unq(l$df) # DAG: 2000 nodes, 20 nbh
# 
# 
# tmp = l$df
# tmp$method = "SHD_rand"
# tmp$SHD = tmp$shd_rand
# tmp2 = l$df
# tmp2$method = "SHD_opt"
# tmp2$SHD = tmp2$SHD_optimal
# dfAll = rbind(l$df, tmp, tmp2)
# dfAll$method = factor(dfAll$method, levels=c("SHD_rand","mean,3simp,TEST","SHD_opt"))
# ggplot(dfAll, aes(factor(dag_nbh), SHD, fill=factor(method))) + geom_boxplot()
# 


# actual plotting
randDag = function(dfAll, opt=FALSE){
  df_randSHD = dfAll[dfAll$method == dfAll$method[1],] # select all entries from one method to get all simulated graphs exactly once
  df_randSHD$method = "shd_rand"
  df_randSHD$SHD = df_randSHD$shd_rand
  if(opt){
    df_opt = df_randSHD
    df_opt$method = "shd_opt"
    df_opt$SHD = df_opt$SHD_opt
    df_randSHD = rbind(df_opt,df_randSHD)
  }
  return(df_randSHD)
}
l1 = readRDS("data_new/gies_comp/1-large_GIES_1.Rds")$df;  # pt, 
l2 = readRDS("data_new/gies_comp/1-large_GIES_3.Rds")$df; # dags, nbh5, GIES
l3 = readRDS("data_new/gies_comp/2-large_GIES_nb_1.Rds")$df;  # dags, nbh1, 
l4 = readRDS("data_new/gies_comp/2-large_GIES_nb_10.Rds")$df;  # dags, nbh10
l5 = readRDS("data_new/gies_comp/ours.Rds"); 
ours = l5[[1]]$df
ours = ours[ours$dag_nbh %in% c(1,5,10,20),]

orderMethods = c("SHD rand.","GIES","P.2, simple, IRC")
dfNbhs = rbind(l2,randDag(l2),l3,randDag(l3),l4,randDag(l4),ours, 
               randDag(ours[ours$dag_nbh == 20,]))
dfNbhs = pretty_names(dfNbhs)
dfNbhs$method = factor(dfNbhs$method, levels=orderMethods)
ggplot(dfNbhs, aes(factor(dag_nbh), SHD, fill=(method))) + geom_boxplot() +
  scale_fill_manual(values=colors_ort, breaks=orderMethods) + 
  labs(title="DAGs with p=500, n=1000, d=21, k=10", fill="Method") +
  theme(plot.title = element_text(face="bold")) +
  xlab(dict_axisLabs("dag_nbh")) + coord_trans(y="log10")







# ---- skel: runtime  ----
data_folder = "data/"  # with final '/'
y = "time_s"
aggrFcts = c("ltest","mean","median")
dfAll = rbind(hd_skel1, hd_skel2, hd_skel3)
ylimits = c(min(dfAll[,y],na.rm=TRUE), max(dfAll[,y],na.rm=TRUE))
hd_skel1 = readRDS(paste(data_folder, "1a_01_runtime.Rds",sep=""))
hd_skel1 = hd_skel1$df[hd_skel1$df$method %in% aggrFcts, ]
f1 = ggplot(hd_skel1, aes(totalSamples, hd_skel1[,y], fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.", title = "p=500, d=20, k=1") + 
  xlab(dict_axisLabs("totalSamples")) + ylab(dict_axisLabs(y)) + ylim(ylimits) +
  scale_fill_manual(values = colors_skel, breaks = aggrFcts)
hd_skel2 = readRDS(paste(data_folder, "1a_02_runtime.Rds",sep=""))
hd_skel2 = hd_skel2$df[hd_skel2$df$method %in% aggrFcts, ]
f2 = ggplot(hd_skel2, aes(ndatasets, time_s, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.", title = "p=500, n=500, k=1") + 
  xlab(dict_axisLabs("ndatasets")) + ylab(dict_axisLabs(y))+ ylim(ylimits) +
  scale_fill_manual(values = colors_skel, breaks = aggrFcts)
hd_skel3 = readRDS(paste(data_folder, "1a_03_runtime.Rds",sep=""))
hd_skel3 = hd_skel3$df[hd_skel3$df$method %in% aggrFcts, ]
f3 = ggplot(hd_skel3, aes(interventionSize, time_s, fill=factor(method))) + geom_boxplot() + 
  labs(fill="Aggr. fct.", title = "p=500, n=500, d=20") + 
  xlab(dict_axisLabs("interventionSize")) + ylab(dict_axisLabs(y))+ ylim(ylimits)+
  scale_fill_manual(values = colors_skel, breaks = aggrFcts)
fig1a_runtime = f1 + f2 + f3 + plot_layout(nrow = 1, guides="collect") & 
  theme(legend.position = "right", plot.title = element_text(face="bold"))
fig1a_runtime


# ---- skel: pooled data ----
# low dimensional
df = readRDS(paste(data_folder, "supp_pooled.Rds",sep=""))$df
df = df[df$interventionSize %in% c(1,6,10),]
df = pretty_names(df)
df$method = factor(df$method, levels = c("baseObs","ltest","mean","median","pooled"))

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}
fig_pooled = ggplot(df, aes(df[,x],SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs(x)) + labs(title="p=20, n=500, d=21",fill="Aggr. Fct.") +
  theme(legend.position = "right", plot.title = element_text(face="bold",hjust=0.5)) +
  scale_fill_manual(values = colors_skel) + 
  theme(strip.text.x = element_text(size = 12))
fig_pooled





# ---- skel: different intervention types ----
# load and process data
data_folder = 'data/'
nds_diff = readRDS("data_new/skeleton_recovery/skeleton_01.Rds")$df # ndatasets
smpl_diff = readRDS("data_new/skeleton_recovery/skeleton_02.Rds")$df # totalsamples
intvS_diff = readRDS("data_new/skeleton_recovery/skeleton_03.Rds")$df # intervSize
intvS_pfct = rbind(skeleton_03$df,skeleton_03_1$df,skeleton_03_2$df) # itnervSize perfect
intvS_pfct = intvS_pfct[intvS_pfct$kindOfIntervention == "perfect",]
smpl_pfct = rbind(skeleton_02$df,skeleton_02_1$df) # totalSamples perfect
smpl_pfct = smpl_pfct[smpl_pfct$kindOfIntervention == "perfect",]
nds_pfct = rbind(skeleton_01$df,skeleton_01_1$df,skeleton_01_2$df) # ndatasets
nds_pfct = nds_pfct[nds_pfct$kindOfIntervention == "perfect",]
nds = pretty_names(rbind(
  nds_diff,
  nds_pfct[nds_pfct$interventionSize == 10,]
))
smpl = pretty_names(rbind(
  smpl_diff[smpl_diff$totalSamples %in% c(500,1000,5000),],
  smpl_pfct[smpl_pfct$totalSamples %in% c(500,1000,5000) & smpl_pfct$interventionSize == 10,]
))
intvS = pretty_names(rbind(
  intvS_diff[intvS_diff$interventionSize %in% c(10,50,100),],
  intvS_pfct[intvS_pfct$interventionSize %in% c(10,50,100),]
))

# plotting
interventions = c("perfect","inhibitory","flipped")
titles = c("p=500, d=21, k=10","p=500, n=1000, k=10", "p=500, n=1000, d=21")
smpl_plot = smpl[smpl$kindOfIntervention %in% interventions,]
smpl_plot$kindOfIntervention = factor(smpl_plot$kindOfIntervention, levels=interventions)
nds_plot = nds[nds$kindOfIntervention %in% interventions,]
nds_plot$kindOfIntervention = factor(nds_plot$kindOfIntervention, levels=interventions)
intvS_plot = intvS[intvS$kindOfIntervention %in% interventions,]
intvS_plot$kindOfIntervention = factor(intvS_plot$kindOfIntervention, levels=interventions)
fig_smpls = ggplot(smpl_plot, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs("totalSamples")) + labs(title=titles[1],fill="Aggr. Fct.") +
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5)) +
  scale_fill_manual(values = colors_skel) + 
  theme(strip.text.x = element_text(size = 13))
fig_nds = ggplot(nds_plot, aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs("ndatasets")) + labs(title=titles[2],fill="Aggr. Fct.") +
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5)) +
  scale_fill_manual(values = colors_skel) + 
  theme(strip.text.x = element_text(size = 13))
fig_intvS = ggplot(intvS_plot, aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs("interventionSize")) + labs(title=titles[3],fill="Aggr. Fct.") +
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5)) +
  scale_fill_manual(values = colors_skel) + 
  theme(strip.text.x = element_text(size = 13))
fig_skel_intv = fig_smpls + fig_nds + fig_intvS + plot_layout(ncol=1,guides="collect") &
  theme(legend.position = "bottom")
fig_skel_intv



#---- ort: diff intervention types ----
# load low-dim data
titles = c("p=20, d=11, k=2","p=20, n=500, k=2", "p=20, n=500, d=11")
ort_nds = readRDS("data_new/Orientation/ld_03-01.Rds")$df; unq(ld1) # nds
ort_smpl = readRDS("data_new/Orientation/ld_03-02.Rds")$df; unq(ld2) # totalSmppl
ort_interv = readRDS("data_new/Orientation/ld_03-03.Rds")$df; unq(ld3) # interv

# load high-dim data
titles = c("p=500, d=21, k=10","p=500, n=1000, k=10", "p=500, n=1000, d=21")
ort_nds = readRDS("data_new/Orientation/03-01.Rds")$df #ndatasets
ort_smpl = readRDS("data_new/Orientation/03-02.Rds")$df # totalsamples
ort_interv = readRDS("data_new/Orientation/03-03.Rds")$df  # intervSize
ort_smpl = ort_smpl[ort_smpl$totalSamples %in% c(500,1000,5000),]
ort_interv = ort_interv[ort_interv$interventionSize %in% c(2,10,20),]

# processing of either LD or HD data (execute either load LD or load HD and then this section)
interventions = c("perfect","inhibitory","flipped")
smpl_plot = ort_smpl[ort_smpl$kindOfIntervention %in% interventions,]
smpl_plot$kindOfIntervention = factor(smpl_plot$kindOfIntervention, levels=interventions)
nds_plot = ort_nds[ort_nds$kindOfIntervention %in% interventions,]
nds_plot$kindOfIntervention = factor(nds_plot$kindOfIntervention, levels=interventions)
intvS_plot = ort_interv[ort_interv$kindOfIntervention %in% interventions,]
intvS_plot$kindOfIntervention = factor(intvS_plot$kindOfIntervention, levels=interventions)

nds_rand = mean(nds_plot[nds_plot$proc == "random","SHD"])
smpl_rand = mean(smpl_plot[smpl_plot$proc == "random","SHD"])
intvS_rand = mean(intvS_plot[intvS_plot$proc == "random","SHD"])

smpl_plot = pretty_names(smpl_plot[smpl_plot$proc != "random", ])
nds_plot = pretty_names(nds_plot[nds_plot$proc != "random",])
intvS_plot = pretty_names(intvS_plot[intvS_plot$proc != "random",])

fig_smpls = ggplot(smpl_plot, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs("totalSamples")) + labs(title=titles[1],fill="Method") +
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5)) +
  geom_hline(yintercept = smpl_rand, linetype = "dashed") +
  scale_fill_manual(values = colors_ort) + 
  theme(strip.text.x = element_text(size = 13))
fig_nds = ggplot(nds_plot, aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs("ndatasets")) + labs(title=titles[2],fill="Method") +
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5)) +
  geom_hline(yintercept = nds_rand, linetype = "dashed") +
  scale_fill_manual(values = colors_ort) + 
  theme(strip.text.x = element_text(size = 13))
fig_intvS = ggplot(intvS_plot, aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs("interventionSize")) + labs(title=titles[3],fill="Method") +
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5)) +
  geom_hline(yintercept = intvS_rand, linetype = "dashed") +
  scale_fill_manual(values = colors_ort) + 
  theme(strip.text.x = element_text(size = 13))
fig_ort_intv = fig_smpls + fig_nds + fig_intvS + plot_layout(ncol=1,guides="collect") &
  theme(legend.position = "bottom")
fig_ort_intv



#---- unbalanced sample sizes ----
load("data_new/unbalanced(2)/unbalanceddata.RData")

# Skeleton
skelplot = skeleton_04
skelplot$sdatasets = as.character(as.numeric(skelplot$sdatasets) / 1000)
p_s<-ggplot(skelplot,
            aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  xlab(expression(paste("w"[o])))+
  theme(legend.position = "bottom",legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(face="bold"))+
  guides(fill=guide_legend(title="")) +
  scale_fill_manual(values=colors_skel)

# Orientation
breaks = c("P.1, refined, BIC",  "P.1, refined, IRC", "P.1, simple, BIC", "P.1, simple, IRC", 
           "P.2, refined, IRC","P.2, refined, BIC", "P.2, simple, IRC",  "P.2, simple, BIC")
orplot = or
orplot$sdatasets = as.character(as.numeric(orplot$sdatasets) / 1000)
p_o<-ggplot(orplot[orplot$method!="random" & orplot$kindOfIntervention == "perfect",],
            aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  geom_hline(yintercept=mm, linetype="dotted")+xlab(expression(paste("w"[o])))+
  theme(legend.position = "bottom",legend.key.size = unit(0.4, 'cm'))+
  theme(axis.title.x = element_text(vjust=-1))+
  guides(fill=guide_legend(title=""))+ylab("") +
  scale_fill_manual(values=colors_ort, breaks=breaks)

# figure for skeleton and orientation
fig_unb_skel_ort = p_s + p_o + plot_layout(nrow=1) + 
  plot_annotation(title="p=500, n=1000, d=21, k=10", 
                  theme=theme(plot.title=element_text(face="bold",hjust=0.5)));  # 1000 by 500

# comparison to GIES
sdaplot = sda
sdaplot$sdatasets = as.character(as.numeric(sdaplot$sdatasets) / 1000)
fig_unb_GIES <- ggplot(sdaplot,aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  xlab(expression(paste("w"[o])))+ labs(title = "p=500, n=1000, d=21, k=10",fill="Method") +
  ylab("SHD")+geom_hline(yintercept=m_dag, linetype="dotted")+
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5))+
  scale_fill_manual(values = colors_ort, breaks=c("GIES","P.2, simple, IRC")) # 450 by 450

