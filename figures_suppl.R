library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

#---- helpful plot functions ----
# dictionary for converting axis labels
dict_axisLabs = function(x){
  # val = switch(x,
  #              "ndatasets" = "number of interv. datasets",
  #              "totalSamples" = "number of samples",
  #              "interventionSize" = "size of interventions",
  #              "time_s" = "time (s)"
  #              )
  val = switch(x,
               "ndatasets" = "d",
               "totalSamples" = "n",
               "interventionSize" = "k",
               "time_s" = "time (s)"
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




#---- explore data ----
load("data_new/many_data (1).RData")
load("data_new/skeleton_data.RData")
dev.new()

#ORIENTATION LD
plot_orientations(ld1, "totalSamples","SHD")
plot_orientations(ld2, "interventionSize","SHD")
plot_orientations(ld3, "ndatasets","SHD")

##Orientation HD
plot_orientations(l03_01, "ndatasets","time_s")
plot_orientations(l03_02, "totalSamples","SHD")
plot_orientations(l03_03, "interventionSize","SHD")


## full algo very HD (only mean,3simp,TEST)
plot_orientations(l05_01,"interventionSize","SHD") # 1000 nodes
plot_orientations(l05_02,"interventionSize","SHD") # 2000 nodes



#SKELETON LD
# base intervsize = 2
plot_skeleton(l_0, "totalSamples","SHD")
plot_skeleton(l_1, "interventionSize","SHD") 
plot_skeleton(l_2, "ndatasets","SHD")   

# intervsize = 4
plot_skeleton(l_3, "totalSamples","SHD") 

# intervSize = 10
plot_skeleton(l_4, "totalSamples","SHD")


##SKELETON HD
l = NULL; l$df<-rbind(skeleton_03$df,skeleton_03_1$df,skeleton_03_2$df); l$str = skeleton_03$str
plot_skeleton(l, "interventionSize","SHD")

l = NULL; l$df<-rbind(skeleton_02$df,skeleton_02_1$df); l$str = skeleton_02$str
plot_skeleton(l, "totalSamples","SHD")

l = NULL; l$df<-rbind(skeleton_01$df,skeleton_01_1$df,skeleton_01_2$df); l$str = skeleton_01$str
plot_skeleton(l, "interventionSize","SHD")



# increase parameters until GIES fails
NULL


#---- explore the skeleton recovery ----
# ld skel: don't think that I need this
ld_skl = readRDS("data_new/skeleton_recovery/ld_skeleton.Rds"); 
unq(ld_skl[[1]]$df)   # intervsize
unq(ld_skl[[2]]$df)   # totalSample
unq(ld_skl[[3]]$df)   # ndatasets

# high-dim
skl_01 = readRDS("data_new/skeleton_recovery/skeleton_01.Rds"); unq(skl_01$df) #ndatasets
skl_02 = readRDS("data_new/skeleton_recovery/skeleton_02.Rds"); unq(skl_02$df) #totalsamples


skl_03 = readRDS("data_new/skeleton_recovery/skeleton_03.Rds"); unq(skl_03$df) #intervSize
skl_04 = readRDS("data_new/skeleton_recovery/skeleton_04.Rds"); unq(skl_04$df) #perfect interv






#---- orientation analysis -----
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

ll1 = readRDS("data_new/Orientation/03-01.Rds"); unq(ll1$df) # ndatasets
ll2 = readRDS("data_new/Orientation/03-02.Rds"); unq(ll2$df) # totalsamples
ll3 = readRDS("data_new/Orientation/03-03.Rds"); unq(ll3$df) # intervSize






# ---- figure: different neighborhoods ----
load_dir = 'old_cluster_comp/'
l6 = readRDS(paste(load_dir,"05-06.Rds",sep=''))
l6_GIES = readRDS(paste(load_dir,"05-06_GIES.Rds",sep=''))
nbh1 = rbind(l6$df, l6_GIES$df)
nbh1 = rbind(nbh1)
l7 = readRDS(paste(load_dir,"05-07.Rds",sep=''))
l7_GIES = readRDS(paste(load_dir,"05-07_GIES.Rds",sep=''))
nbh5 = rbind(l7$df, l7_GIES$df)
nbh5 = rbind(nbh5)
l8 = readRDS(paste(load_dir,"05-08.Rds",sep=''))
nbh10 = rbind(l8$df)


dfAll = rbind(nbh1, nbh5, nbh10)
y = "SHD"
means = c("1" = 1000, "5" =2000, "10"=3000)
ggplot(dfAll, aes(totalSamples, dfAll[,y], fill=method)) + geom_boxplot() + ylab(y)+
  geom_hline(aes(yintercept = value), data=means) +  
  facet_grid(cols = ggplot2::vars(dag_nbh), labeller = labeller(.default = label_both))

ylims = c(min(dfAll$SHD), max(dfAll$SHD))
nbh1_rand = mean(nbh1[nbh1$method == "shd_rand","SHD"])
nbh1_plt = nbh1[ !(nbh1$method %in% c("shd_opt","shd_rand")),]
# nbh1_plt = pretty_names(nbh1_plt)
fnbh1 = ggplot(nbh1_plt, aes(totalSamples, nbh1_plt[,y], fill=factor(method))) + geom_boxplot()+
  xlab(dict_axisLabs("totalSamples")) + ylab(dict_axisLabs(y)) + ylim(ylims)+
  geom_hline(yintercept = nbh1_rand, linetype = "dotted") + 
  labs(title = "p=500, d=21, k=10, h=1")
nbh5_rand = mean(nbh5[nbh5$method == "shd_rand","SHD"])
nbh5_plt = nbh5[ !(nbh5$method %in% c("shd_opt","shd_rand")),]
# nbh5_plt = pretty_names(nbh5_plt)
fnbh5 = ggplot(nbh5_plt, aes(totalSamples, nbh5_plt[,y], fill=factor(method))) + geom_boxplot()+
  xlab(dict_axisLabs("totalSamples")) + ylab(dict_axisLabs(y)) +ylim(ylims)+
  geom_hline(yintercept = nbh5_rand, linetype = "dotted") + 
  labs(title = "p=500, d=21, k=10, h=5")
nbh10_rand = mean(nbh10[nbh10$method == "shd_rand","SHD"])
nbh10_plt = nbh10[ !(nbh10$method %in% c("shd_opt","shd_rand")),]
# nbh10_plt = pretty_names(nbh10_plt)
# nbh10_plt$method = factor(nbh10_plt$method, levels = c("GIES","P.1,simple,IRC","P.2,simple,IRC"))
fnbh10 = ggplot(nbh10_plt, aes(totalSamples, nbh10_plt[,y], fill=factor(method))) + geom_boxplot()+
  xlab(dict_axisLabs("totalSamples")) + ylab(dict_axisLabs(y)) +ylim(ylims)+
  geom_hline(yintercept = nbh10_rand, linetype = "dotted") + 
  labs(title = "p=500, d=21, k=10, h=10")
fnbhs = fnbh1 + fnbh5 + fnbh10 + plot_layout(nrow=1, guides="collect") &
  theme(plot.title = element_text(face="bold"))
fnbhs



# all mean,3simp,TEST
l = readRDS("data_new/Results/05-vlarge_011.Rds"); unq(l$df) # PT: 1000 nodes, totalSamples
l = readRDS("data_new/Results/05-vlarge_012.Rds"); unq(l$df) # DAG: 1000 nodes, nbh=1,5,10,20
  # vlarge for GIES did not finish
  #  ATM: we are running vlarge with GIES

l = readRDS("data_new/Results/05-large_011.Rds"); unq(l$df) # DAG, 1000 nodes, 20 nbh
# l = readRDS("data_new/Results/05-large_022.Rds"); unq(l$df) # DAG: 2000 nodes, 20 nbh


tmp = l$df
tmp$method = "SHD_rand"
tmp$SHD = tmp$shd_rand
tmp2 = l$df
tmp2$method = "SHD_opt"
tmp2$SHD = tmp2$SHD_optimal
dfAll = rbind(l$df, tmp, tmp2)
dfAll$method = factor(dfAll$method, levels=c("SHD_rand","mean,3simp,TEST","SHD_opt"))
ggplot(dfAll, aes(factor(dag_nbh), SHD, fill=factor(method))) + geom_boxplot()



# ---- runtime skeleton recovery ----
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


# ---- pooled data ----
# low dimensional
df = readRDS(paste(data_folder, "supp_pooled.Rds",sep=""))$df
df = df[df$interventionSize %in% c(1,6,10),]
df$method = factor(df$method, levels = c("baseObs","ltest","mean","median","baseAll"))

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}
fig_pooled = ggplot(df, aes(df[,x],SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=ggplot2::vars(kindOfIntervention), labeller=labeller(.default = capitalize))+
  xlab(dict_axisLabs(x)) + labs(title="p=20, n=500, n=21",fill="Aggr. Fct.") +
  theme(legend.position = "right", plot.title = element_text(face="bold")) +
  scale_fill_manual(values = colors_skel) + 
  theme(strip.text.x = element_text(size = 30))
fig_pooled

# x = "interventionSize"
# ylim_ld = df[df$kindOfIntervention %in% c("perfect","flipped"),"SHD"]
# ylim_ld = c(min(ylim_ld), max(ylim_ld))
# df_ld_pfct = df[df$kindOfIntervention == "perfect",]
# fpool_ld_pfct = ggplot(df_ld_pfct, aes(df_ld_pfct[,x],SHD,fill=factor(method))) + 
#   geom_boxplot() + ylim(ylim_ld) + 
#   xlab(dict_axisLabs(x)) + labs(title="perfect, p=20, n=500, n=21",fill="Aggr. Fct.") +
#   theme(plot.title = element_text(face="bold")) + scale_fill_manual(values = colors_skel)
# df_ld_flip = df[df$kindOfIntervention == "flipped",]
# fpool_ld_flip = ggplot(df_ld_flip, aes(df_ld_flip[,x],SHD,fill=factor(method))) + 
#   geom_boxplot() + ylim(ylim_ld) + ylab(NULL)+
#   xlab(dict_axisLabs(x)) + labs(title="flipped, p=20, n=500, n=21",fill="Aggr. Fct.") +
#   theme(plot.title = element_text(face="bold")) + scale_fill_manual(values = colors_skel)
# 
# # high-dim
# skl_03 = readRDS("data_new/skeleton_recovery/skeleton_03.Rds"); unq(skl_03$df) #intervSize
# # skl_04 = readRDS("data_new/skeleton_recovery/skeleton_04.Rds"); unq(skl_04$df) #perfect interv
# df_hd = rbind(skl_03$df, skeleton_03$df)
# df_hd = df_hd[df_hd$interventionSize %in% c(1,10,20,50,100,250),]
# df_hd$method = factor(df_hd$method, levels = c("baseObs","ltest","mean","median","baseAll"))
# ylim_hd = df_hd[df_hd$kindOfIntervention %in% c("perfect","flipped"),"SHD"]
# ylim_hd = c(min(ylim_hd), max(ylim_hd))
# df_hd_pfct = df_hd[df_hd$kindOfIntervention == "perfect",]
# fpool_hd_pfct = ggplot(df_hd_pfct, aes(df_hd_pfct[,x],SHD,fill=factor(method))) + 
#   geom_boxplot() + ylim(ylim_hd) +
#   xlab(dict_axisLabs(x)) + labs(title="perfect, p=500, n=1000, n=21",fill="Aggr. Fct.") +
#   theme(plot.title = element_text(face="bold")) + scale_fill_manual(values = colors_skel)
# df_hd_flip = df_hd[df_hd$kindOfIntervention == "flipped",]
# fpool_hd_flip = ggplot(df_hd_flip, aes(df_hd_flip[,x],SHD,fill=factor(method))) + 
#   geom_boxplot() + ylim(ylim_hd)+ ylab(NULL) +
#   xlab(dict_axisLabs(x)) + labs(title="flipped, p=500, n=1000, n=21",fill="Aggr. Fct.") +
#   theme(plot.title = element_text(face="bold")) + scale_fill_manual(values = colors_skel)
# 
# # merge
# fpool = fpool_ld_pfct + fpool_ld_flip + fpool_hd_pfct + fpool_hd_flip +
#   plot_layout(ncol = 2, guides="collect"); fpool



# ---- skel different intervention types ----
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

nds = rbind(
  nds_diff,
  nds_pfct[nds_pfct$interventionSize == 10,]
)
smpl = rbind(
  smpl_diff[smpl_diff$totalSamples %in% c(500,1000,5000),],
  smpl_pfct[smpl_pfct$totalSamples %in% c(500,1000,5000) & smpl_pfct$interventionSize == 10,]
)
intvS = rbind(
  intvS_diff[intvS_diff$interventionSize %in% c(10,50,100),],
  intvS_pfct[intvS_pfct$interventionSize %in% c(10,50,100),]
)


# plot skeleton data nicely
plot_skel = function(df, x, y, ylims, title_str){
  fig = ggplot(df, aes(df[,x], df[,y], fill=factor(method))) +
    geom_boxplot() + ylim(ylims) +
    xlab(dict_axisLabs(x)) + ylab(dict_axisLabs(y)) +
    labs(title = title_str, fill="Aggr. fct.")+
    scale_fill_manual(values = colors_skel, breaks = aggrFcts) +
    theme(plot.title = element_text(face = "bold"))
  return(fig)
}

# plotting
cols = c("perfect","inhibitory","random")
rows = c("totalSamples","ndatasets","interventionSize")
titles = c("p=500, d=21, k=10","p=500, n=1000, k=10", "p=500, n=1000, d=21")
aggrFcts = c("baseObs","ltest","mean","median","pooled")
rowDfs = list(pretty_names(smpl), pretty_names(nds), pretty_names(intvS))
plotDfs = list()
figs = vector("list",length=length(cols)*length(rows))
i = 1
for(r in 1:length(rows)){
  ylimits = c(min(rowDfs[[r]]$SHD), max(rowDfs[[r]]$SHD))
  for(col in 1:length(cols)){
    df = rowDfs[[r]][rowDfs[[r]]$kindOfIntervention == cols[col], ]
    figs[[i]] = plot_skel(df, rows[r], "SHD", ylimits, paste(cols[col],", ",titles[r],sep=""))
    if(col != 1){
      figs[[i]] = figs[[i]] + ylab(NULL) + theme(axis.title.y=element_blank(),
                                                 axis.text.y=element_blank(),
                                                 axis.ticks.y=element_blank())
    }
    i = i + 1
  }
}
fig = Reduce("+", figs[-1], figs[[1]]) + plot_layout(ncol = 3, guides="collect") &
  theme(legend.position = "bottom")
fig


fig1 = figs[[1]]+labs(title="perfect") + 
  figs[[2]] +labs(title="inhibitory") +
  figs[[3]] +labs(title="random") +
  plot_layout(nrow=1, guides="collect") +
  plot_annotation(title=titles[1], theme=theme(
            plot.title = element_text(hjust = 0.5, face="bold"))) &
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face="bold")); fig1



#---- ort: diff intervention types ----



