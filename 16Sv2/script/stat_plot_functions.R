# 程序功能：R语言包安装、统计、绘图大全
# Functions: Frequency used packages,  statistics and plotting functions in R
# 主要内容 Main steps: 
# - 1. 常用包自动安装和加载
# - 2. 绘图参数和函数
# - 3. 统计函数



# 1. 常用包自动安装和加载

## 1.1 安装CRAN来源常用包

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("dplyr","vegan", "reshape2", "ggplot2", "devtools", "ggthemes", "agricolae","pheatmap","plyr") # ,"multcompView","ggpubr"
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# ## 1.2 安装bioconductor常用包
#  
# package_list = c("ggrepel")
# for(p in package_list){
#   if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
#     source("https://bioconductor.org/biocLite.R")
#     biocLite(p)
#     suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
#   }
# }
# 
# # 1.3 安装Github常用包
# 
# package_list = c("kassambara/ggpubr")
# for(p in package_list){
#   q=unlist(strsplit(p,split = "/"))[2]
#   if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
#     install_github(p)
#     suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
#   }
# }


  
# 2. 绘图参数和函数
  
  # Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.

## 2.1 设置颜色

# alpha  =  .7
# c_yellow  =           rgb(255 / 255, 255 / 255,   0 / 255, alpha)
# c_blue  =             rgb(  0 / 255, 000 / 255, 255 / 255, alpha)
# c_orange  =           rgb(255 / 255,  69 / 255,   0 / 255, alpha)
# c_green  =            rgb(  50/ 255, 220 / 255,  50 / 255, alpha)
# c_dark_green  =       rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
# c_very_dark_green  =  rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
# c_sea_green  =        rgb( 46 / 255, 129 / 255,  90 / 255, alpha)
# c_black  =            rgb(  0 / 255,   0 / 255,   0 / 255, alpha)
# c_grey  =             rgb(180 / 255, 180 / 255,  180 / 255, alpha)
# c_dark_brown  =       rgb(101 / 255,  67 / 255,  33 / 255, alpha)
# c_red  =              rgb(200 / 255,   0 / 255,   0 / 255, alpha)
# c_dark_red  =         rgb(255 / 255, 130 / 255,   0 / 255, alpha)

## 2.2 调置ggplot2主题

fontsize=8
# 主题为空白，轴线宽和颜色，刻度和文字，图例位置，正文文字
main_theme = theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(size = .5, colour = "black"),
    axis.line.y = element_line(size = .5, colour = "black"),
    axis.ticks = element_line(size = .5, color = "black"),
    axis.text = element_text(size = fontsize, color = "black"),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = fontsize),
    text = element_text(family = "sans", size = fontsize))




#  3. 统计函数

## 3.1 summarySE：计算样本均值、标准差、标准误和置信区间
summarySE = function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # 计算长度
  length2 = function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac = ddply(data, groupvars, .drop=.drop,.fun = function(xx, col) {
                  c(N    = length2(xx[[col]], na.rm=na.rm),
                    mean = mean   (xx[[col]], na.rm=na.rm),
                    sd   = sd     (xx[[col]], na.rm=na.rm))},
                    measurevar)
  # 重命名  
  datac = plyr::rename(datac, c("mean" = measurevar))
  # 计算标准偏差
  datac$se = datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult = qt(conf.interval/2 + .5, datac$N-1)
  datac$ci = datac$se * ciMult
  return(datac)
}


## 3.2 解析Vegan::cca结果
variability_table  =  function(cca){
  chi  =  c(cca$tot.chi,
           cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table  =  cbind(chi, chi/chi[1])
  colnames(variability_table)  =  c("inertia", "proportion")
  rownames(variability_table)  =  c("total", "constrained", "unconstrained")
  return(variability_table)
}

cap_var_props  =  function(cca){
  eig_tot  =  sum(cca$CCA$eig)
  var_propdf  =  cca$CCA$eig/eig_tot
  return(var_propdf)
}

pca_var_props  =  function(cca){
  eig_tot  =  sum(cca$CA$eig)
  var_propdf  =  cca$CA$eig/eig_tot
  return(var_propdf)
}

cca_ci  =  function(cca, permutations=5000){
  var_tbl  =  variability_table(cca)
  p  =  permutest(cca, permutations=permutations)
  ci  =  quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
  return(ci)
}



## 3.3 三元图绘制函数
tern_e = function (x, scale = 1, dimnames = NULL, dimnames_position = c("corner",
                                                                       "edge", "none"), dimnames_color = "black", id = NULL, id_color = "black",
                  coordinates = FALSE, grid = TRUE, grid_color = "gray", labels = c("inside",
                                                                                    "outside", "none"), labels_color = "darkgray", border = "black",
                  bg = "white", pch = 19, cex = 1, prop_size = FALSE, col = "red",
                  main = "ternary plot", newpage = TRUE, pop = TRUE, ...)
{
  labels  =  match.arg(labels)
  if (grid == TRUE)
    grid  =  "dotted"
  if (coordinates)
    id  =  paste("(", round(x[, 1] * scale, 1), ",", round(x[,
                                                            2] * scale, 1), ",", round(x[, 3] * scale, 1), ")",
                sep = "")
  dimnames_position  =  match.arg(dimnames_position)
  if (is.null(dimnames) && dimnames_position != "none")
    dimnames  =  colnames(x)
  if (is.logical(prop_size) && prop_size)
    prop_size  =  3
  if (ncol(x) != 3)
    stop("Need a matrix with 3 columns")
  if (any(x < 0))
    stop("X must be non-negative")
  s  =  rowSums(x)
  if (any(s <= 0))
    stop("each row of X must have a positive sum")
  x  =  x/s
  top  =  sqrt(3)/2
  if (newpage)
    grid.newpage()
  xlim  =  c(-0.03, 1.03)
  ylim  =  c(-1, top)
  pushViewport(viewport(width = unit(1, "snpc")))
  if (!is.null(main))
    grid.text(main, y = 0.9, gp = gpar(fontsize = 18, fontstyle = 1))
  pushViewport(viewport(width = 0.8, height = 0.8, xscale = xlim,
                        yscale = ylim, name = "plot"))
  eps  =  0.01
  grid.polygon(c(0, 0.5, 1), c(0, top, 0), gp = gpar(fill = bg,
                                                     col = border), ...)
  if (dimnames_position == "corner") {
    grid.text(x = c(0, 1, 0.5), y = c(-0.02, -0.02, top +
                                        0.02), label = dimnames, gp = gpar(fontsize = 12))
  }
  if (dimnames_position == "edge") {
    shift  =  eps * if (labels == "outside")
      8
    else 0
    grid.text(x = 0.25 - 2 * eps - shift, y = 0.5 * top +
                shift, label = dimnames[2], rot = 60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.75 + 3 * eps + shift, y = 0.5 * top +
                shift, label = dimnames[1], rot = -60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.5, y = -0.02 - shift, label = dimnames[3],
              gp = gpar(col = dimnames_color))
  }
  if (is.character(grid))
    for (i in 1:4 * 0.2) {
      grid.lines(c(1 - i, (1 - i)/2), c(0, 1 - i) * top,
                 gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(1 - i, 1 - i + i/2), c(0, i) * top,
                 gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(i/2, 1 - i + i/2), c(i, i) * top, gp = gpar(lty = grid,
                                                               col = grid_color))
      if (labels == "inside") {
        grid.text(x = (1 - i) * 3/4 - eps, y = (1 - i)/2 *
                    top, label = i * scale, gp = gpar(col = labels_color),
                  rot = 120)
        grid.text(x = 1 - i + i/4 + eps, y = i/2 * top -
                    eps, label = (1 - i) * scale, gp = gpar(col = labels_color),
                  rot = -120)
        grid.text(x = 0.5, y = i * top + eps, label = i *
                    scale, gp = gpar(col = labels_color))
      }
      if (labels == "outside") {
        grid.text(x = (1 - i)/2 - 6 * eps, y = (1 - i) *
                    top, label = (1 - i) * scale, gp = gpar(col = labels_color))
        grid.text(x = 1 - (1 - i)/2 + 3 * eps, y = (1 -
                                                      i) * top + 5 * eps, label = i * scale, rot = -120,
                  gp = gpar(col = labels_color))
        grid.text(x = i + eps, y = -0.05, label = (1 -
                                                     i) * scale, vjust = 1, rot = 120, gp = gpar(col = labels_color))
      }
    }
  xp  =  x[, 2] + x[, 3]/2
  yp  =  x[, 3] * top
  size = unit(if (prop_size)
    #emiel inserted this code. x are proportions per row.  x*s is original data matrix. s = rowsums of original data matrix (x*s)
    prop_size * rowSums(x*x*s) / max(  rowSums(x*x*s) )
    #prop_size * rowSums(    (x*s) * ((x*s)/s)) / max(  rowSums(    (x*s) * ((x*s)/s)) )
    else cex, "lines")
  grid.points(xp, yp, pch = pch, gp = gpar(col = col), default.units = "snpc",
              size = size, ...)
  if (!is.null(id))
    grid.text(x = xp, y = unit(yp - 0.015, "snpc") - 0.5 *
                size, label = as.character(id), gp = gpar(col = id_color,
                                                          cex = cex))
  if (pop)
    popViewport(2)
  else upViewport(2)
}



## 3.4 da_adonis：距离矩阵adonis组间差异统计
# Compare each group distance matrix by vegan adonis in bray_curtis
da_adonis = function(sampleV){
  sampleA = as.matrix(sampleV$sampA)
  sampleB = as.matrix(sampleV$sampB)
  design2 = subset(sub_design, group %in% c(sampleA,sampleB))
  if (length(unique(design2$group))>1) {
    sub_dis_table = dis_table[rownames(design2),rownames(design2)]
    sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
    adonis_table = adonis(sub_dis_table ~ group, data = design2, permutations = 1000) 
    adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
    print(paste("In", m, "pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
    adonis_pvalue = paste(m, sampleA, sampleB, adonis_pvalue, sep="\t")
    return(adonis_pvalue)
  }
}



## 3.5 相关评估r2
lm_eqn  =  function(df){
  m  =  lm(y ~ x, df);
  eq  =  substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



## 3.6 样品稀释曲线抽样统计函数 sample rarefracation curve

# 传入参数为OTU表

sample_rare  =  function(df, count_cutoff =3, length = 30, rep = 30){  
  # 小于一定频率是认为是噪间 set threshold to distinguish noise, default cutoff = 3, min 1, max 30
  # count_cutoff = 3
  # loop for 1 to n samples, set length
  # length = 30
  # if combo > 30, only sample 30 times to reduce calculate
  # rep = 30
  
  # Set result as final result
  result = data.frame(sample = c(0), richness = c(0))
  # Set inital value of otu table
  count=df
  # Set otu table to binary for easy calculate richness
  
  count[count < count_cutoff] = 0
  count[count >= count_cutoff] = 1
  # Sample number
  n = dim(count)[2]
  
  x = unique(as.integer(seq(1, n, length.out = length)))
  for (i in x){
    # i = 1
    # get combination, need transposition and format to df
    # combn组合太多，时间过长，如combn(30,15)的数量为choose(30,15)即1.5亿行
    # list = as.data.frame(t(combn(n, i)))

    # choose calcuate combn, then sample get list for each rep
    m = choose(n, i)
    if (m > rep){m = rep}
    
    # loop list and calculate richness
    for(j in 1:m) {
      idx = sample(n, i)
      temp = count[, idx, drop = F]
      # subset non-zero row
      temp1=temp[rowSums(temp)>0, , drop=F]
      # row number is accumulative OTUs
      result = rbind(result, c(i, dim(temp1)[1]))
    }
  }
  # remove start 0,0
  result = result[-1,]
  # factor draw as each box
  result$sample = factor(result$sample, levels = unique(result$sample))
  return(result)
}



## 3.7 基于数据框，默认的group列和批量的m列绘制箱线图+统计

boxplot_anova  =  function(index, m){  
  # boxplot stat code
model = aov(index[[m]] ~ group, data=index)
Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_table = as.data.frame(Tukey_HSD$group) 
# write.table(paste(m, "\n\t", sep=""), file=paste(wd, b, m,  ".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# suppressWarnings(write.table(Tukey_HSD_table, file=paste(wd, b, m,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

out = LSD.test(model,"group", p.adj="none") # alternative fdr
stat = out$groups
index$stat=stat[as.character(index$group),]$groups
max=max(index[,c(m)])
min=min(index[,c(m)])
x = index[,c("group",m)]
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',m,')',sep=""))
y=as.data.frame(y)
rownames(y)=y$group
index$y=y[as.character(index$group),]$Max + (max-min)*0.05

p = ggplot(index, aes(x=group, y=index[[m]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste(m, "")) + theme_classic() + main_theme +
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+theme(legend.position="none")
if (length(unique(index$group))>5){	p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
p
}


## 3.8 单列数据绘制饼形图，要求列V1为名称，V2为数值

plot_pie  =  function(df, radius, threshold){  
  # 绘制柱状图、饼图代码
  data = df
  width = radius
  colnames(data)=c("class", "CK")
  # 生成自定义图例标签，在标签上显示百分比；
  lab1 =as.vector(data$class)
  # lab1
  lab1 =paste(lab1, "(", round(data$CK, 2), ")",sep = "") # /sum(data$CK)*100 # "%)"
  # lab1
  # 生成仅有“百分比”的数据标签，并将数值小于threshold，如3%的标签设置为“空”，避免占比较小的区域标签互相重叠；
  lab2 = round(data$CK/sum(data$CK)*100,1)
  # lab2
  n = length(lab2)
  # n
  for(i in 1: n){
    if(as.numeric(lab2[i])< threshold )
      lab2[i] = ""
    else
      lab2[i] = paste(lab2[i],"%",sep= "")
  }
  lab2
  # 绘制堆叠条形图，设置条形图的“边框”为白色；
  # p1 = ggplot(data=data,aes(x="",y=data$CK,fill=(data$class)))+geom_bar(stat="identity",width=0.4,color="white",linetype=1,size=1)
  # p1
  # #添加百分比标签；
  # p2 = p1+geom_text(aes(x=1,label=lab2),position= position_stack(reverse =F,vjust=0.5),size=6)
  # p2
  # coord_polar()将直角坐标系转为极坐标系；
  p3 = ggplot(data=data,aes(x="",y=data$CK,fill=data$class))+ 
    geom_bar(stat="identity",width = width,color="white",linetype=1,size=1)+
    coord_polar(theta="y", start = pi/2) + labs(x="",y="",title="")
  p3 = p3 + geom_text(aes(x=1.25,label=lab2),position = position_stack(reverse =F, vjust=0.5), size=4)
  # vjust的值取 0位于底部，0.5中间， 1 (the default) 在上部；
  # p3
  # 更改图表的主题，去掉横纵坐标标题，添加图片标题，将背景改为白色，图列位置（这里仍保持右侧）；
  # 更新为自定义的图例标签，使标签显示百分比；
  p4 = p3 + theme_bw() + 
    theme(axis.text=element_blank(),panel.border=element_blank(),axis.ticks=element_blank(), panel.grid=element_blank(),legend.title= element_blank(), legend.position = "right") +
    scale_fill_discrete(breaks= data$class, labels = lab1)
  p4
}

## 3.9 Tukey检验转换字母
# 将Tukey检验结果P值转换为显著字母分组
generate_label_df = function(TUKEY, variable){
  library(multcompView)
  # 转换P值为字母分组
  ## 提取图基检验中分组子表的第4列P adjust值
  Tukey.levels = TUKEY[[variable]][,4]
  ## multcompLetters函数将两两p值转换为字母，data.frame并生成列名为Letters的数据框
  Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
  # 按分组名字母顺序
  ## 提取字母分组行名为group组名
  Tukey.labels$group = rownames(Tukey.labels)
  # 按组名的字线顺序排列，默认的Levels
  Tukey.labels=Tukey.labels[order(Tukey.labels$group) , ]
  return(Tukey.labels)
}


## 双Y轴折线图

# ##使用gtable + grid 包中的一些命令，组合ggplot2作图结果，得到双坐标y轴图
# library(ggplot2)
# library(gtable)
# library(grid) 

##定义组合函数
ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  library(gtable)
  library(grid)
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

# #使用的 R 自带的 mtcars 数据集，分别绘制“mpg-disp”与“mpg-drat”折线图
# g1 <- ggplot(mtcars, aes(mpg, disp)) + geom_line(colour = 'blue') + labs(y = 'disp') + theme_bw()
# g2 <- ggplot(mtcars, aes(mpg, drat)) + geom_line(colour = 'red') + labs(y = 'drat') + theme_bw() %+replace% 
#   theme(panel.background = element_rect(fill = NA))
# 
# #使用 ggplot2.two_y_axis() 整合结果，并保存图片
# ggsave('g1.pdf', g1)
# ggsave('g2.pdf', g2)
# ggsave('g12.pdf', ggplot2.two_y_axis(g1, g2))