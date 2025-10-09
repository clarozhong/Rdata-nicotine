#关注基因箱型图----
library(ggplot2)
library(reshape2)
library(ggpubr)
library(stringr)
library(reshape2)
g<-c("AKT1", "JUN", "ALB","TP53", "STAT3", 
"CASP3", 
"IL1B", 
"TGFB1", 
"FOS", 
"TNF", 
"PTGS2", 
"CTNNB1" )

exp<-`TCGA-LUAD_tpm`
#exp$symbol<-rownames(exp)
#exp$symbol <- trimws(exp$symbol)
IF_exp <- exp[rownames(exp) %in% g,]


group<- `TCGA-LUAD_group`

rownames(group)<-colnames(IF_exp)
rt<-t(IF_exp)
rt1 <- cbind(group, rt)
rt1<-rt1[,-1]
rt_long=melt(rt1)     #melt表示拆分数据
colnames(rt_long)=c('Group','Gene','Expression') #设置拆分后列名
head(rt_long)

pdf("表达箱线图.pdf",width = 10,height = 6)
ggplot(rt_long, aes(x = Gene, y = Expression, fill = Group)) +
  #geom_violin(trim = FALSE)  +  # 使用箱型图
  geom_boxplot(width = 0.8, outlier.shape = NA) +  # 可选：在小提琴图上叠加一个小型箱型图
  labs(x = "TCGA-LUAD", y = "Expression",title = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # 去除网格线
    axis.line = element_line(color = "black"),  # 显示坐标轴线
    axis.title.x = element_text(face = "bold"),  # x轴标题样式
    axis.title.y = element_text(face = "bold")   # y轴标题样式
  ) +  
  scale_fill_manual(values = c("yellowgreen", "violetred1"))+
  #scale_fill_manual(values = c("#1F78B4", "#C51B7D"))+
  stat_compare_means(label.y = 14,
                     method = "t.test",
                     label="p.signif")
dev.off()

#12个基因GSEA----
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)
library(ggplot2)
library(ggridges)
dat_res_diff <- readRDS('diff_res.Rds')
dat_res_diff <- dat_res_diff[rownames(dat_res_diff) %in% g,]
# 将基因名与logFC对应起来
logFC_all_gene <- rlang::set_names(dat_res_diff$log2FoldChange, rownames(dat_res_diff))
# 按照logFC对基因进行排序
logFC_all_gene <- sort(logFC_all_gene, decreasing = T)
# 基因集（GSEA一般用C2）
geneset_c2 <- read.gmt('c2.all.v2024.1.Hs.symbols.gmt')
geneset_c2 <- read.gmt('h.all.v2024.1.Hs.symbols.gmt')
# GSEA分析
res_gsea <- GSEA(
  logFC_all_gene,
  # 基因集
  TERM2GENE = geneset_c2,
  seed = 2024,
  # p值矫正方法
  pAdjustMethod = 'BH',
  # p值阈值
  pvalueCutoff = 0.05
)
# 提取GSEA分析结果
dat_res_gsea <- as.data.frame(res_gsea)

#LASSO筛选----
library(pROC)
# 假设train_data是训练数据集，group是目标变量
set.seed(123)
rt1<-train
# 随机分配数据到训练集和测试集，这里以80%训练，20%测试为例

sampleIndex <- sample(1:nrow(rt1), 0.7 * nrow(rt1))

# 训练集
train_data <- rt1[sampleIndex, ]

# 测试集
test_data <- rt1[-sampleIndex, ]
#运行LASSO
library(glmnet)#加载glmnet包

#colnames(mydata[,1:17])#查看前17列的列名（根据自己数据调整）
y <- as.matrix(train_data[, 1])  # 提取第1列作为结局（建议放在第一列）
x <- as.matrix(train_data[, 2:13])  # 第2至第17列为自变量
t_y<-as.matrix(test_data[, 1])
t_x<-as.matrix(test_data[, 2:13])
#x<- as.numeric(x)
#后边的代码除了s值基本不需更改
lasso_model <- glmnet(x, y, family = "binomial",
                      alpha = 1) # 表示采用L1正则化，即Lasso回归。
max(lasso_model$lambda)
print(lasso_model) 
#绘制LASSO图
pdf("LASSO系数路径图.pdf",width = 8,height = 6)
plot(lasso_model,
     xvar = "lambda")
dev.off()

#交叉验证并绘制可视化结果
pdf("LASSO交叉验证图.pdf",width = 8,height = 6)
cv_model <- cv.glmnet(x, y, family = "binomial",alpha = 1,nfolds = 10)
plot(cv_model)
dev.off()
#根据交叉验证结果，选择lambda值，lambda.min或lambda.1se。
lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se

#s为Lambda大小，Lambda越大表示模型的正则化强度越大，选择的自变量也越少。
#这里选择的是刚刚得到的lambda_min的值
coef_lasso <- coef(lasso_model,
                   s =  lambda_1se)
coef_lasso
#结果显示后边带有数值的变量为筛选得到的变量
# 使用最优lambda值重新训练模型
lasso_model_final <- glmnet(x, y, family = "binomial", alpha = 1, lambda = lambda_1se)

# 预测测试集
predictions <- predict(lasso_model_final, newx = t_x, s = lambda_1se, type = "response")

# 计算ROC曲线和AUC值
roc_obj <- roc(response = t_y, predictor = predictions)
#### ROC曲线的绘制----
pdf("LASSO_ROC.pdf",width = 8,height = 6)
plot(roc_obj,
     print.auc=TRUE, #设置是否添加AUC值标签
     auc.polygon=T, #设置是否添加AUC值面积多边形
     grid=c(0.1, 0.2), #设置是否添加网格线
     grid.col=c("#B6D7CE","#E8BBC5"), ##设置网格线颜色
     col="red",#设置ROC曲线颜色
     max.auc.polygon=TRUE, #设置是否添加最大AUC值面积多边形
     legacy.axes=F, #x轴格式更改
     auc.polygon.col="skyblue", #设置AUC值面积多边形的填充色
     print.thres=T, #是否添加截点和95%CI，
     print.thres.col="black",#设置阈值的颜色
     main="ROC curve")#图的主标题
dev.off()
#GSVA----
library(GSVA)
library(tidyverse)
library(limma)
library(ggplot2)
library(pheatmap)
exp<-`TCGA-LUAD_tpm`
dat_group <- `TCGA-LUAD_group`
c<-c("AKT1", "JUN", "ALB","TP53", "STAT3", "CASP3", "IL1B",  "TGFB1", "FOS", "TNF", "PTGS2",  "CTNNB1",
     "PHLPP1", "PHLPP2", "RGCC", "APPL1", "FOXO4", "RPS6KB2", "AKT2", "PTEN", "THEM4", "RICTOR",
     "SLCO1B3", "SLCO1B1", "SLCO1A2", "SLC10A1", "LCAT", "ABCC3", "FABP6", "HPX", "AFM", "GC",
     "CAPG", "ACIN1", "APPL1", "DFFA", "STK3", "GAS2", "DFFB", "DIABLO", "STK26", "DSG1",
     "LEF1", "CTNNA1", "CTNNBIP1", "LRRFIP1", "TCF7L2", "CDH1", "APC", "SOX1", "SOX17", "CHD8",
     "JUN", "IGFBP7", "JUNB", "KDM6B", "ATF4", "ELK1", "FOSB", "FOSL2", "RELA", "DUSP1",
     "MAP3K3", "IL1R1", "IL1R2", "A2M", "IL1RAP", "CASP1", "SQSTM1", "IL1A", "IRAK3", "MAPK8IP2",
     "FOS", "ATF3", "MAPK8", "IGFBP7", "MAPK9", "ATF4", "CREB5", "KDM6B", "ATF2", "SENP2",
     "PTGIS", "TBXAS1", "PTGDS", "PTGS1", "HTR1B", "ALOX5AP", "CAV1", "LTC4S", "NUCB1", "ELAVL1",
     "EGFR", "NFKBIZ", "STAT5B", "PTK2B", "DUT", "PIAS3", "STAT1", "IL17F", "PDIA3", "IL6ST",
     "TGFBR3", "TGFBR2", "LTBP1", "LTBP4", "MMP9", "TGFBR1", "NOS2", "LRRC32", "HIF1A", "CEBPA",
     "TNFRSF1A", "TRADD", "TNFRSF1B", "CEBPA", "RIPK1", "TRAF2", "ADAM17", "CHUK", "TNFAIP3", "HNRNPA1",
     "MDM2", "MYB", "COP1", "SIN3A", "TP53BP1", "TP53BP2", "MDM4", "TP63", "MAP3K5", "PMAIP1")

c<-c( "JUN", "CASP3", "IL1B",  "FOS", 
     "CAPG", "ACIN1", "APPL1", "DFFA", "STK3", "GAS2", "DFFB", "DIABLO", "STK26", "DSG1",
     "JUN", "IGFBP7", "JUNB", "KDM6B", "ATF4", "ELK1", "FOSB", "FOSL2", "RELA", "DUSP1",
     "MAP3K3", "IL1R1", "IL1R2", "A2M", "IL1RAP", "CASP1", "SQSTM1", "IL1A", "IRAK3", "MAPK8IP2",
     "FOS", "ATF3", "MAPK8", "IGFBP7", "MAPK9", "ATF4", "CREB5", "KDM6B", "ATF2", "SENP2")

dat_expr <- exp[rownames(exp) %in% c,]

# GSVA分析 ---------------------------------------------------------------

# 基因集（GSVA一般用HALLMARK）
geneset_h <- GSEABase::getGmt('h.all.v2024.1.Hs.symbols.gmt')
# GSVA分析
dat_res_gsva <- gsvaParam(as.matrix(dat_expr),
                          # 基因集
                          geneset_h,
                          # 富集得分计算方式为最大正负随机游走偏差之间的幅度差
                          maxDiff = T,
                          # 当输入的表达谱为Count（整数）时选择'Poisson'，其余情况均用'Gaussian'
                          kcdf = 'Gaussian') %>% gsva() %>% as.data.frame()


# GSVA差异分析 ---------------------------------------------------------------

# 创建模型矩阵
design <- model.matrix( ~ 0 + factor(dat_group$group))
colnames(design) <- levels(factor(dat_group$group))
rownames(design) <- dat_group$sample

# 使用线性模型拟合富集得分数据
fit <- lmFit(dat_res_gsva, design)

# 为组间差异创建对照
contrasts <- c('LUAD-Control')
# 一定要疾病组在前，对照组在后

# 创建线性对比矩阵
cont_matrix <- makeContrasts(contrasts = contrasts, levels = design)

# 对拟合的线性模型进行再次拟合
fit2 <- contrasts.fit(fit, cont_matrix)

# Bayes调整
fit2 <- eBayes(fit2)

# 提取差异分析结果
dat_res_diff_gsva <- topTable(fit2, coef = contrasts, n = Inf) %>% na.omit()

# 将为0的p值改成除0外最小的p值，以防止p值为0取不到对数
dat_res_diff_gsva$P.Value[dat_res_diff_gsva$P.Value == 0] <- min(dat_res_diff_gsva$P.Value[dat_res_diff_gsva$P.Value != 0])
dat_res_diff_gsva$adj.P.Val[dat_res_diff_gsva$adj.P.Val == 0] <- min(dat_res_diff_gsva$adj.P.Val[dat_res_diff_gsva$adj.P.Val != 0])

# 筛选差异显著的通路
dat_res_diff_gsva_sig <- dat_res_diff_gsva %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::arrange(desc(t))
dat_res_diff_gsva_sig$t_group <- ifelse(dat_res_diff_gsva_sig$t < 0, 'Down', 'Up')
# GSVA双向柱状图 ---------------------------------------------------------------

# 双向柱状图数据
dat_res_diff_gsva_sig$pathway <- factor(rownames(dat_res_diff_gsva_sig),
                                        levels = rownames(dat_res_diff_gsva_sig))

# 双向柱状图
p <- ggplot(dat_res_diff_gsva_sig, aes(t, pathway, fill = t_group)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  guides(fill = 'none') +
  labs(y = element_blank()) +
  scale_fill_manual(values = c('#4DBBD5', '#E64B35'))

ggsave(file = '39geneGSVA_barplot.pdf', p, width = 8, height = 8)
rm(list=ls())
library(gwasglue)
library(ieugwasr)
library(dplyr)
library(tidyr)
library(CMplot)
library(TwoSampleMR)
library(gwasvcf)
library(VariantAnnotation)
##暴露因素
setwd("~/LXY/吸烟")
exposure<-readVcf("ieu-b-142.vcf.gz")
exposure<-gwasvcf_to_TwoSampleMR(vcf=exposure)##OK
exposure<-exposure[,c(1:9,11,13,16)]
colnames(exposure)<-c("chr.exposure","pos.exposure","other_allele.exposure","effect_allele.exposure",
                     "beta.exposure","se.exposure","pval.exposure","eaf.exposure","samplesize.exposure","SNP","exposure","id.exposure")
#去连锁不平衡
e_snp<-exposure %>%dplyr::select(rsid=SNP, pval=pval.exposure)
setwd("~/LXY")
e_clumped <- ieugwasr::ld_clump(e_snp,
                                clump_kb = 10000,
                                clump_r2 = 0.001,
                                clump_p = 1,
                                plink_bin = plinkbinr::get_plink_exe(),
                                bfile = "/home/xiyouyun-lqh/LXY/EUR/EUR",
                                pop = "EUR")
exposure<-exposure[exposure$SNP%in%e_clumped$rsid,]
###筛选强相关工具变量
exposure1<-subset(exposure,pval.exposure<5e-8)
###根据F筛选强工具变量
exposure1$MAF<-ifelse(exposure1$eaf.exposure<0.5,exposure1$eaf.exposure,1-exposure1$eaf.exposure)
exposure1$R2<-2*exposure1$beta.exposure^2*exposure1$MAF*(1-exposure1$MAF)/((exposure1$se.exposure)^2*exposure1$samplesize.exposure)
exposure1$F<-exposure1$R2*(1-exposure1$R2)*(exposure1$samplesize.exposure-2)
exposure2<-subset(exposure1,F>10&MAF>0.05)
##结局
##输入在IEU网站申请的令牌
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxMjE5OTg3MjQwQHFxLmNvbSIsImlhdCI6MTc0Mzc2MTM1NiwiZXhwIjoxNzQ0OTcwOTU2fQ.Dn7nwLKR9MU6EZVIRrxEGixeP0Qi71hVF2510ASGYvErK4a6ekTYECE41vnmprnBQuKfWjxSpCJ1WsYy8aaNW8IVNQ4KxYNxDPquFzi9RaaWDw5kfAkqs1u_dqBvD6ZZAQfQ8zZRNi09f-r7nmkSnzRVX38IrilAqijT7l8Pf54tre3dKYzrNq2P0kM2uYUJj8o4A0fzcUh4YEvgi2uwO-p8tQ4xlPkzLvUpDqIUX1w5XgZ8uPhUQNQGcnytdVW7L2fN-enCtFMIYr8Sk-Sb9In2FXOolKERaYg9yF3bhoEzEYMnD4Ir8W4_R2AkMTvvXRjWHYkMNTEhgx4jkT5FhQ")
##肺腺癌的结局数据编号ebi-a-GCST004744
outcome<-extract_outcome_data(snps=exposure2$SNP,outcomes ="ebi-a-GCST004744",proxies = F)
##匹配,寻找共有SNP
mrdata<-harmonise_data(exposure_dat = exposure,
                       outcome_dat = outcome,
                       action = 2)
{mrdata[,6]<-as.numeric(mrdata[,6])
  mrdata[,7]<-as.numeric(mrdata[,7])
  mrdata[,16]<-as.numeric(mrdata[,16])
  mrdata[,18]<-as.numeric(mrdata[,18])
  mrdata[,26]<-as.numeric(mrdata[,26])
  mrdata[,27]<-as.numeric(mrdata[,27])}
##steiger过滤，确保因果方向，若访问得到的数据没有samplesize则需要在数据库里自行查看并在下面输入
#mrdata$samplesize.outcome<-xxx
mrdata<-steiger_filtering(mrdata)
mrdata<-mrdata[mrdata$steiger_dir==T,]
#####再GWAS catalog里搜索各个snp对应的性状是否为可能的混杂因素并将之除去
mrdata<-mrdata[c(4,9,11,12,13),]##选定的SNP仅被报道和吸烟情况有关
#####正式mr分析
res<-mr(mrdata)
#异质性分析
het<-mr_heterogeneity(mrdata,method_list = c('mr_ivw'))##pval<0.05,IVW法需要换用随机效应模型
res1<-mr(mrdata,method_list = "mr_ivw_mre")
#多效性检验
Pleio<-mr_pleiotropy_test(mrdata)
save(res,res1,Pleio,file="吸烟对肺腺癌的因果.RData")

load("~/LXY/吸烟对肺腺癌的因果(1).RData")
####单SNP森林图
res1<-mr_singlesnp(mrdata)##单个snp的mr数据
res1$lo_ci<-res1$b-1.96*res1$se
res1$up_ci<-res1$b+1.96*res1$se
res1$b_0.95a<-ifelse(is.na(res1$b),"",
                     sprintf("%.3f(%.3f to %.3f)",res1$b,res1$lo_ci,res1$up_ci))
res1$` `<-paste(rep(" ",40),collapse=" ")
colnames(res1)[6]<-""
###画法
tm<-forest_theme(base_size=10,##图层大小
                 refline_col="red",
                 core = list(
                   bg_params = list(fill = "white"),  # 背景色
                   fg_params = list( hjust = 0.5,x = 0.5)),
                 ci_alpha = 0.7,#透明度
                 ci_pch = 22,##箱心形状,15~25时才可以填充
                 ci_col = "#4575b4",##箱边框颜色
                 ci_fill = "#9975b4",##箱填充颜色
                 vertline_lwd=0.6,#额外垂线宽度
                 #vertline_lty=，#额外垂线类型
                 vertline_col="#4575b4",#额外垂线颜色
                 ci_lty=1,
                 ci_lwd=2,
                 ci_Theight = 0.2,
                 arrow_type="closed",
                 footnote_gp=gpar(col="blue",cex=1))
p<-forest(res1[,c(6,13)],
          est=res1$b,
          lower=res1$lo_ci,
          upper=res1$up_ci,
          sizes=0.5,
          ci_column=2,
          ref_line=0,
          xlim=c(-2.5,2.5),
          ticks_a=c(-2,-1,0,1,2),
          nudge_y = 0.05,##间距
          #footnote="P<0.05 was considered\nstatistically significant",
          title=" Cigarettes smoked per day associations to Lung adenocarcinoma",
          theme=tm)
p
{#为标题部分的顶部添加下划线，分隔不同部分
  g<-edit_plot(p,row=1:23,which="background",gp=gpar(fill="#ffffff"))
  g<-add_border(g,part="header",row=2,col=2,where="top")
  #为标题部分的底部添加下划线
  g<-add_border(g,part="header",row=24,col=2,where="bottom")
  #在第9行的指定列添加上边框,可以多加
  g<-add_border(g,row=1:23,col=2,where="left")
  g<-add_border(g,row=1:23,col=2,where="right")
  g<-edit_plot(g,row = 22:23,col = 2,which="ci",
               gp=gpar(col="#aa2266",fill="#112264"))
  g<-add_border(g,row=22,col=2,where="top",gp=gpar(col="grey"))
}
g
pdf("SNP森林图.pdf",width = 8, height = 10 )
par(mar=c(1, 1, 1, 1)) 
plot(g)##每一条曲线代表一个变量的回归过程
dev.off()
####留一法森林图
res2 <- mr_leaveoneout(mrdata)
res2$lo_ci<-res2$b-1.96*res2$se
res2$up_ci<-res2$b+1.96*res2$se
res2$b_0.95a<-ifelse(is.na(res2$b),"",
                     sprintf("%.3f(%.3f to %.3f)",res2$b,res2$lo_ci,res2$up_ci))
res2$`  `<-paste(rep(" ",40),collapse=" ")
res2$SNP<-paste0("Without ",res2$SNP)
res2$SNP[22]<-"All"
colnames(res2)[6]<-""
###画法
tm<-forest_theme(base_size=10,##图层大小
                 refline_col="red",
                 core = list(
                   bg_params = list(fill = "white"),  # 背景色
                   fg_params = list( hjust = 0.5,x = 0.5)),
                 ci_alpha = 0.7,#透明度
                 ci_pch = 22,##箱心形状,15~25时才可以填充
                 ci_col = "#45aa77",##箱边框颜色
                 ci_fill = "#9975b4",##箱填充颜色
                 vertline_lwd=0.6,#额外垂线宽度
                 #vertline_lty=，#额外垂线类型
                 vertline_col="#4575b4",#额外垂线颜色
                 ci_lty=1,
                 ci_lwd=2,
                 ci_Theight = 0.2,
                 arrow_type="closed",
                 footnote_gp=gpar(col="blue",cex=1))
p<-forest(res2[,c(6,13)],
          est=res2$b,
          lower=res2$lo_ci,
          upper=res2$up_ci,
          sizes=0.5,
          ci_column=2,
          ref_line=0,
          xlim=c(-0.1,0.9),
          ticks_a=c(0,0.25,0.50,0.75),
          nudge_y = 0.05,##间距
          #footnote="P<0.05 was considered\nstatistically significant",
          title="                                        Leave one out plot",
          theme=tm)
p
{#为标题部分的顶部添加下划线，分隔不同部分
  g<-edit_plot(p,row=1:22,which="background",gp=gpar(fill="#ffffff"))
  g<-add_border(g,part="header",row=2,col=2,where="top")
  #为标题部分的底部添加下划线
  g<-add_border(g,part="header",row=23,col=2,where="bottom")
  #在第9行的指定列添加上边框,可以多加
  g<-add_border(g,row=1:22,col=2,where="left")
  g<-add_border(g,row=1:22,col=2,where="right")
  g<-edit_plot(g,row = 22,col = 2,which="ci",
               gp=gpar(col="#aa2266",fill="#112264"))
  g<-add_border(g,row=22,col=2,where="top",gp=gpar(col="grey"))
}
g
pdf("去一法森林图.pdf",width = 8, height = 10 )
par(mar=c(1, 1, 1, 1)) 
plot(g)##每一条曲线代表一个变量的回归过程
dev.off()



load("~/LXY/吸烟对肺腺癌的因果(1).RData")
res_single<-mr_singlesnp(mrdata)
p<-mr_scatter_plot(res,mrdata) 
# 散点图
data<-p[["RF2Ajh.ebi-a-GCST004744"]][["data"]]
p1<-p[[1]] + 
  theme_bw() +  # 白色背景 + 灰色网格线
  theme(
    panel.background = element_rect(fill = "white"),  # 背景色
    panel.grid.major = element_line(color = "grey90"), # 主网格线颜色
    panel.grid.minor = element_blank()                # 移除次要网格线
  )
data2<-p1[["layers"]][[4]][["data"]]
p<-p1 + geom_point(color = "#ff55aa",fill="#226699" ,shape = 15, size = 1)+  # 蓝色三角形点，大小=3
  geom_errorbarh(data = data,
    aes(xmin = beta.exposure-se.exposure, xmax = beta.exposure+se.exposure),  # 替换为你的数据列名（如 lci, uci）
    color = "#9999cc",alpha=0.7)+
  geom_errorbar(data =data,
    aes(ymin = beta.outcome-se.outcome, ymax = beta.outcome+se.outcome),  # 替换为你的数据列名（如 lci, uci）
    color = "#9999cc",alpha=0.7)+
  geom_abline(data = data2,
    aes(intercept = a,slope = b,color = method),
    linewidth = 1.2)+
  scale_color_manual(values = c("#99dddd","#005500","#44c4aa","#77cc66","#ffdd44"))+
  labs(x = "SNP effect on daily smoking volumn",             # X轴标题
    y = "SNP effect on lung adenocarchinoma")+
  theme(axis.title = element_text(size = 14))
pdf("散点图.pdf",width = 10, height = 8 )
par(mar=c(1, 1, 3, 3)) 
plot(p)##每一条曲线代表一个变量的回归过程
dev.off()
#漏斗图
p4 <- mr_funnel_plot(res_single)
data<-p4[["RF2Ajh.ebi-a-GCST004744"]][["data"]]
p<-ggplot(data,aes(x=b,y=se))+
  geom_point(color ="#ff1111" ,shape = 4, size = 3, stroke = 1.2)+
  theme_bw() +  # 白色背景 + 灰色网格线
  theme(
    panel.background = element_rect(fill = "white"),  # 背景色
    panel.grid.major = element_line(color = "grey90"), # 主网格线颜色
    panel.grid.minor = element_blank()                # 移除次要网格线
  )+
  geom_vline(p4[["RF2Ajh.ebi-a-GCST004744"]][["layers"]][[2]][["data"]][1,],
             mapping=aes(xintercept = b,color = SNP),            # 颜色
             linewidth = 1,linetype = "dashed")+
  scale_color_manual(values = c("#9999cc","#99dddd"))+
  labs(
    x = "Effect Size",             # X轴标题
    y = "Standard Error",   # Y轴标题
    title = "Funnel Plot(MR)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "plain"), # 标题居中加粗
    axis.title = element_text(size = 14),
    legend.position = c(0.99, 0.25),  # 左上角坐标
    legend.justification = c(1, 1),   # 对齐方式
    legend.background = element_rect(fill = "white", color = "grey80"))   # 坐标轴标题大小
p
pdf("漏斗图.pdf",width = 7, height = 7 )
par(mar=c(1, 1, 1, 1)) 
plot(p)##每一条曲线代表一个变量的回归过程
dev.off()
