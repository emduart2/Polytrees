# library(RColorBrewer)
# 
l1 = readRDS("cluster_computation/01_results/l101.Rds")
l2 = readRDS("cluster_computation/01_results/l102.Rds")

ind = c(3,4,7,8)
reds = brewer.pal(n = 9, name = 'Reds')[c(3,4,7,8)]
blues = brewer.pal(n = 9, name = 'Greens')[c(3,4,7,8)]
greens = brewer.pal(n = 9, name = 'Blues')[c(3,4,7,8)]
pattern = c(4,3,2,1)
colors1 = c(reds[pattern], blues[pattern], greens[pattern])

# dict = function(x){return(switch(x,
#   "gtruth,1,BIC" = "P.1, BIC, refined",
#   "gtruth,1,TEST" = "P.1, TEST, refined",
#   "gtruth,1simp,BIC" = "P.1, BIC, simple",
#   "gtruth,1simp,TEST" = "P.1, TEST, simple",
#   "gtruth,2,BIC" = "P.2, BIC, refined",
#   "gtruth,2,TEST" = "P.2, TEST, refined",
#   "gtruth,2simp,BIC" = "P.2, BIC, simple",
#   "gtruth,2simp,TEST" = "P.2, TEST, simple",
#   "gtruth,3,BIC" = "P.3, BIC, refined",
#   "gtruth,3,TEST" = "P.3, TEST, refined",
#   "gtruth,3simp,BIC" = "P.3, BIC, simple",
#   "gtruth,3simp,TEST" = "P.3, TEST, simple",
# ))}
dict = function(x){return(switch(x,
   "gtruth,1,BIC" = "P.1, refined, BIC",
   "gtruth,1,TEST" = "P.1, refined, TEST",
   "gtruth,1simp,BIC" = "P.1, simple, BIC",
   "gtruth,1simp,TEST" = "P.1, simple, TEST",
   "gtruth,2,BIC" = "P.2, refined, BIC",
   "gtruth,2,TEST" = "P.2, refined, TEST",
   "gtruth,2simp,BIC" = "P.2, simple, BIC",
   "gtruth,2simp,TEST" = "P.2, simple, TEST",
   "gtruth,3,BIC" = "P.3, refined, BIC",
   "gtruth,3,TEST" = "P.3, refined, TEST",
   "gtruth,3simp,BIC" = "P.3, simple, BIC",
   "gtruth,3simp,TEST" = "P.3, simple, TEST",
))}
l1$df$method = sapply(l1$df$method, dict, USE.NAMES = FALSE)
l2$df$method = sapply(l2$df$method, dict, USE.NAMES = FALSE)
legend_str = "Method"

ggplot(l1$df, aes(ndatasets, l1$df[, y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+
  ggplot(l2$df, aes(totalSamples, l2$df[,y], fill=factor(method))) + geom_boxplot()+
  ylab(y) + scale_fill_manual(values=colors1) + labs(fill=legend_str)+
  plot_layout(ncol = 1, guides="collect") 

# size: 900 by 430



# OLD COLOR CODE
# display.brewer.pal(n = 9, name = 'YlOrRd')
# display.brewer.pal(n = 9, name = 'YlGnBu')
# display.brewer.pal(n = 9, name = 'YlGn')
# reds_yell = brewer.pal(n = 9, name = 'YlOrRd')[ind]
# blues_yell = brewer.pal(n = 9, name = 'YlGnBu')[ind]
# greens_yell = brewer.pal(n = 9, name = 'YlGn')[ind]
# colors2 = c(reds[4],reds_yell[3],reds[2],reds_yell[1], 
#             blues[4],blues_yell[3],blues[2],blues_yell[1],
#             greens[4],greens_yell[3],greens[2],greens_yell[1])
# colorsBorder = c("#000000","#FFCC00","#000000","#FFCC00",
#                  "#000000","#FFCC00","#000000","#FFCC00",
#                  "#000000","#FFCC00","#000000","#FFCC00")

# with different colors
# ggplot(l1$df, aes(ndatasets, SHD, fill=method, color=method)) + geom_boxplot() +
#   scale_fill_manual(values=colors1) + scale_color_manual(values=colorsBorder)

