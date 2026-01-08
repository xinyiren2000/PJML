rm(list = ls())
setwd("C:/Users/Administrator/Desktop/PJML/r")

source("gendata.r")
source("pjml.r")
source("others.r")
source("kernel.r")
source("numfac.r")
source("evaluation.r")

load_pkgs <- c("doParallel", "foreach")
sapply(load_pkgs, require, character = TRUE)

## returns names of functions loaded in current session
.export <- unclass(lsf.str())
.packages <- c("MASS", "mvtnorm", "base", "PDSCE")

sapply(.packages, require, character = TRUE)



#read data--------------------------------------------------------------------
datainf <- read.csv("real data/us_macroeconometrics/old/datainf.csv")[1:4]
DI <- subset(datainf,E.F..==1,select=c(Mnemonic,Trans..Code))

dat1 <- read.csv("real data/us_macroeconometrics/old/dat1.csv")
data1 <- (dat1[1:574,-1]+dat1[2:575,-1]+dat1[3:576,-1])/3
data1 <- data1[seq(1,574,3),]

dat2 <- read.csv("real data/us_macroeconometrics/old/dat2.csv")
dat <- cbind(dat2[-1],data1)


tf <- function(x,type){
  l = nrow(x)
  if(type==1){
    y = x[3:l,]
  }else if(type==2){
    y = x[3:l,]-x[2:(l-1),]
  }else if(type==3){
    del1 = x[2:l,]-x[1:(l-1),]
    y = del1[2:(l-1)]-del1[1:(l-2)]
  }else if(type==4){
    y = log(x[3:l,])
  }else if(type==5){
    y = log(x[3:l,])-log(x[2:(l-1),])
  }else if(type==6){
    logdel1 = log(x[2:l,])-log(x[1:(l-1),])
    y = logdel1[2:(l-1)]-logdel1[1:(l-2)]
  }else{
    y = 'wrong'
  }
  return(y)
}

data <- data.frame(matrix(NA,190,109))
colnames(data) <- DI$Mnemonic
rownames(data) <- dat2[c(-1,-2),1]
for(i in 1:nrow(DI)){
  data[,i] <- tf(dat[DI$Mnemonic[i]],DI$Trans..Code[i])
}

Ts = nrow(data)
p = ncol(data)
h = (2.35/sqrt(12)) * Ts^(-1/5) * p^(-1/10)

Ydat = scale(as.matrix(data))

# time-varying factor model estimation ------------------------------------------------------------------------------------------
qopt = select_q(Ydat)

wls = local_pca(Y=Ydat,q=qopt)
lams.tv = CVPJML.tv(Y=Ydat,q=qopt,lam1.max=0.4)
apjml.tv = PJML.tv(Y=Ydat,q=qopt,lam1=lams.tv$lam1,lam2=lams.tv$lam2,ap=TRUE)
pjml.tv = PJML.tv(Y=Ydat,q=qopt,lam1=lams.tv$lam1,lam2=lams.tv$lam2,ap=FALSE,maxit=5)


# regression -------------------------------------------------------------------------------------------------------------------
gdp = scale(tf(dat2[,2,drop=F],5))

model_wls0 = lm(gdp ~ 0 + wls$pre_F[,1]+wls$pre_F[,2])
summary(model_wls0)

model_pjml0 = lm(gdp ~ 0 + pjml.tv$Fhat[,1]+pjml.tv$Fhat[,2])
summary(model_pjml0)

model_apjml0 = lm(gdp ~ 0 + apjml.tv$Fhat[,1]+apjml.tv$Fhat[,2])
summary(model_apjml0)




# heatmap ----------------------------------------------------------------------------------------------------------------------
# 安装并加载必要的包
#install.packages("ggplot2")
#install.packages("reshape2")
library(ggplot2)
library(reshape2)
library(patchwork)
# library(cowplot)
# library(gridExtra)
#library(latex2exp)
library(extrafont)
font_import()
loadfonts(device='win')


# 获取t=1时的图像
i = 1 #129-166

# 创建示例数据
data <- t(wls$B_tilde[,,i])#matrix(rnorm(100), nrow=10
rownames(data) <- c("f1", "f2")
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pwls = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==1)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'), plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x = paste0("t=",i), fill = "Value")

# 创建示例数据
data <- t(pjml.tv$Bhat[,,i])#matrix(rnorm(100), nrow=10)
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pjtv = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==1)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'), plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 创建示例数据
data <- t(apjml.tv$Bhat[,,i])#matrix(rnorm(100), nrow=10)
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pajtv = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==1)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'),  plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 不显示图例和横坐标
pwls1 = pwls + theme(legend.position = 'none', axis.title.y = element_blank()) 
pjtv1 = pjtv + theme(legend.position = 'none', axis.title.y = element_blank()) 
pajtv1 = pajtv + theme(legend.position = 'none', axis.title.y = element_blank())



# 获取t=96时的图像
i = 96 #129-166

# 创建示例数据
data <- t(wls$B_tilde[,,i])#matrix(rnorm(100), nrow=10
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pwls = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==96)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'),  plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 创建示例数据
data <- t(pjml.tv$Bhat[,,i])#matrix(rnorm(100), nrow=10)
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pjtv = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==96)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'),  plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 创建示例数据
data <- t(apjml.tv$Bhat[,,i])#matrix(rnorm(100), nrow=10)
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pajtv = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==96)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'),  plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 不显示图例和横坐标
pwls96 = pwls + theme(legend.position = 'none', axis.title.y = element_blank()) + ggtitle('WLS')
pjtv96 = pjtv + theme(legend.position = 'none', axis.title.y = element_blank()) + ggtitle('PJML')
pajtv96 = pajtv + theme(legend.position = 'none', axis.title.y = element_blank()) + ggtitle('APJML')


# 获取t=190时的图像
i = 190 #129-166

# 创建示例数据
data <- t(wls$B_tilde[,,i])#matrix(rnorm(100), nrow=10
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pwls = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==190)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'),  plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 创建示例数据
data <- t(pjml.tv$Bhat[,,i])#matrix(rnorm(100), nrow=10)
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pjtv = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==190)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'),  plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 创建示例数据
data <- t(apjml.tv$Bhat[,,i])#matrix(rnorm(100), nrow=10)
rownames(data) <- c('f1', 'f2')
# colnames(data) <- colnames(Ydat)
# 转换数据为数据框
data_melt <- melt(data)
# 自定义颜色渐变，调整颜色范围
my_colors <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2))
# 绘制热力图
pajtv = ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  my_colors +
  scale_x_discrete(labels=expression(italic(f)[1],italic(f)[2]), name=expression(italic(t)==190)) +
  theme_minimal() +
  theme(text = element_text(family = 'Times New Roman'),  plot.margin = margin(t=0,b=0),
        axis.text.y = element_blank(),  # 省略行坐标
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(angle = 0, vjust = 0.5, size=10),
        #legend.position = 'top',legend.direction = 'horizontal',
        plot.title = element_text(margin = margin(b=0), size=10),
        legend.key.size = unit(0.5, "cm"),  # 调整图例大小
        legend.text = element_text(size = 7),  # 调整图例文字大小
        legend.title = element_text(size = 10)) +  # 调整图例标题大小
  labs(y = "U.S. macroeconomic data", x =  paste0("t=",i), fill = "Value")

# 不显示图例和横坐标
pwls190 = pwls + theme(legend.position = 'none', axis.title.y = element_blank())
pjtv190 = pjtv + theme(legend.position = 'none', axis.title.y = element_blank())
pajtv190 = pajtv + theme(legend.position = 'none', axis.title.y = element_blank())


# 组合热力图
com_plot = pwls1 | pwls96 | pwls190 | 
  pjtv1 | pjtv96 | pjtv190 | 
  pajtv1 | pajtv96 | pajtv190 + 
  plot_layout(guides = 'collect') & theme(legend.position = 'right') 
#final_plot = com_plot + plot_annotation()
ggsave(filename = 'real data/heatmap.png',plot=com_plot,device='png',width=6,height=7.5,units='in',dpi=600,bg='white')


