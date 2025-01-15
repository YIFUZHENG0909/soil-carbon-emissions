#安装 plspm 包
#install.packages('plspm')
#devtools::install_github('gastonstat/plspm')

#加载 plspm 包
library(plspm)
# 设置工作目录
setwd("D:/碳通量数据/data上传/Fig5")
#读取数据
dat <- read.delim('data2.txt', sep = '\t')
#指定潜变量，在 R 中以列表（list）存储变量和潜变量的关系
###############################6个地区
###############################6个地区
###############################6个地区
#您可以直接指定列名称，或者指定列的下标都可以，我个人习惯指定列名称
dat_blocks <- list(
    Factory = 'Distance',
    soil = c('TC',	'CN'),
    Actinobacteria =  'Actinobacteria',
    shannon =  'Shannon',
    Gene_CD =  c('K01179',	'K01178',	'K01176',	'K01209',	'K01191',	'K06113',	'K01205'),
    cf = "cf"
)
dat_blocks


#通过 0-1 矩阵描述潜变量之间的关联，其中 0 代表变量间没有关联，1 代表有关联
########88888888888888888888888888888888888888888888888888888
Factory<- c(0,0,0,0,0,0)
soil<- c(1,0,0,0,0,0)
Actinobacteria<- c(1,1,0,0,0,0)
shannon<- c(0,1,1,0,0,0)
Gene_CD<- c(0,1,1,1,0,0)
cf<- c(0,0,0,0,1,0)
########66666666666666666666666666666666666666666666

#dat_path <- rbind(TC,CN,pH,shannon,beta,Streptomyces,K01187,K02274,K02276,K05343,K02275,K00121,cf)
dat_path <- rbind(Factory,soil,Actinobacteria,shannon,Gene_CD,cf)
colnames(dat_path) <- rownames(dat_path)
dat_path


#指定因果关系，可选 A（代表列是行的因） 或 B（代表行是列的因）
dat_modes <- rep('A', 6)
dat_modes

##一个简单的 PLS-PM，更多参数详情 ?plspm
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes)
dat_pls
summary(dat_pls)

#结果内容比较多，细节部分还需自行参阅 plspm 包的用户手册：
#完整版手册，235页：https://www.gastonsanchez.com/PLS_Path_Modeling_with_R.pdf
#简版手册，10页：https://rdrr.io/cran/plspm/f/inst/doc/plspm_introduction.pdf

#以下仅展示了一部分相对重要的内容

#查看路径系数的参数估计值，以及相关的统计信息
dat_pls$path_coefs
dat_pls$inner_model

#查看因果关系的路径图，详情 ?innerplot
innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)

#查看作为外源潜变量和内源潜变量的状态
dat_pls$inner_summary

#查看变量间的影响状态
dat_pls$effects

#查看观测变量和潜变量关系，可通过 outerplot() 画图展示类似路径图的结构，详情 ?outerplot
dat_pls$outer_model
outerplot(dat_pls, what = 'loadings', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')
outerplot(dat_pls, what = 'weights', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')

#goodness-of-fit 值可以帮助评估模型优度
dat_pls$gof

#查看潜变量得分，可以理解为标准化后的潜变量的值
dat_pls$scores

#输出潜变量的值
#latent <- data.frame(dat_pls$scores)
#latent <- cbind(dat$site, latent)
#write.csv(latent, 'latent.csv')

