setwd("~/ZYJ/LUAD/TCGA")
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

pdf("表达箱线图.pdf",width = 15,height = 6)
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
#热图----
library(pheatmap)
group_vector <- c(rep("Control", 59), rep("LUAD", 541))
group_list = factor(group_vector,levels = c("Control","LUAD"))  #设置为因子变量
table(group_list)
group<-as.data.frame(group_list)
rownames(group)<-colnames(exp)
colnames(group)<-c("Group")
color_vector <- colorRampPalette(c("yellowgreen","#FAE7D9","violetred1"))(10)
#color_vector <- colorRampPalette(c("#1f78b4","#FAE7D9","#C85D4D"))(10)
p1<-pheatmap(IF_exp,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             clustering_distance_cols = "correlation",  # 设置列聚类的距离度量
             clustering_method = "complete",  # 设置列聚类的方法，这里使用完全链接法
             treeheight_row = 0,  # 不显示行聚类树
             legend = TRUE,
             color = color_vector,
             scale='row',
             show_rownames = T,
             show_colnames = F,
             annotation_col = group)

p1<-pheatmap(IF_exp,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         legend = TRUE,
         legend_breaks = c(-5, 0, 5), #设置图例的断点
         legend_labels = c("low","","high"), #设置图例断点处的标签
         scale='row',
         color = color_vector,
         show_rownames = T,
         show_colnames = F,
         annotation_col = group)
ggsave(file = 'pheatmap.pdf', p1, width = 8, height = 4)
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

ggsave(file = '39geneGSVA_barplot.pdf', p, width = 10, height = 13)

#######################################################################################
#GEO dataset
#######################################################################################
library(GEOquery)
library(Biobase)
library(stringr)
library(limma)
library(tidyverse)
library(FactoMineR)
library(factoextra)

# 1. 设置路径与环境 ---------------------------------------------------------
input_dir  <- "~/ZYJ/LUAD/GEO"
output_root <- "~/ZYJ/LUAD/output"
if(!dir.exists(output_root)) dir.create(output_root)
setwd(input_dir)

# 定义安全 Max 函数
safe_max <- function(x) {
  if(all(is.na(x))) return(NA)
  x <- as.numeric(x)
  return(max(x, na.rm = TRUE))
}

# 2. 准备数据集列表 ---------------------------------------------------------
gse_names <- c("GSE72094")

# --- 循环前的准备工作：获取 GPL570 注释表 ---
cat("正在预先加载 GPL570 (Affymetrix) 注释信息...\n")
gpl570 <- getGEO("GPL570", destdir = input_dir, getGPL = FALSE) 
ids <- Table(gpl570)[, c("ID", "Gene Symbol")]
colnames(ids) <- c("probe_id", "symbol")

# 数据清洗：取第一个基因名
ids <- ids %>% 
  filter(symbol != "" & !is.na(symbol)) %>% 
  mutate(symbol = str_extract(symbol, "^([^/// ]+)"))

# 3. 循环处理数据集 ---------------------------------------------------------
for (gse_id in gse_names) {
  cat("\n--- 正在处理数据集:", gse_id, "---\n")
  
  # 3.1 加载数据
  file_name <- paste0(gse_id, "_series_matrix.txt.gz")
  gset <- getGEO(filename = file_name, getGPL = FALSE, destdir = input_dir)
  if (class(gset) == "list") gset <- gset[[1]]
  
  clin <- pData(gset)
  exp_raw <- exprs(gset)
  
  # 3.2 定制分组逻辑 ----------------------------------------------------
  group_list <- rep(NA, nrow(clin))
  
   if (gse_id == "GSE31210") {
    # 在 GSE31210 中，正常样本通常在 source_name 或 characteristics 中标注
    target_col <- clin[["source_name_ch1"]]
    group_list[grepl("primary lung tumor", target_col, ignore.case = TRUE)] <- "LUAD"
    group_list[grepl("normal lung", target_col, ignore.case = TRUE)] <- "Normal"
    # 注意：如果该数据集全是 tumor，那么我们需要找包含 normal 的数据集
  }
  
  # 3.3 样本过滤 ------------------------------------------------------------
  group_list <- factor(group_list, levels = c("Normal", "LUAD"))
  keep_samples <- !is.na(group_list)
  
  if (sum(keep_samples) < 2) {
    cat(gse_id, " 分组识别失败，跳过。\n")
    next
  }
  
  exp <- exp_raw[, keep_samples, drop = FALSE]
  group_list <- group_list[keep_samples]
  clin_sub <- clin[keep_samples, ]
  
  cat(gse_id, " 分组成功！统计：\n")
  print(table(group_list))
  
  # 3.4 数据标准化与 ID 转换 -------------------------------------------------
  
  # 强制转换为矩阵并数值化，防止 dimnames 报错
  exp <- as.matrix(exp)
  exp_numeric <- apply(exp, 2, as.numeric)
  if(is.null(dim(exp_numeric))) {
    exp_numeric <- matrix(exp_numeric, ncol = ncol(exp))
  }
  rownames(exp_numeric) <- rownames(exp)
  colnames(exp_numeric) <- colnames(exp)
  exp <- exp_numeric
  
  exp[is.na(exp)] <- 0
  
  # 对数化检查
  qx <- quantile(exp, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  if ((qx[5] > 100) || (qx[6] - qx[1] > 50)) {
    exp <- log2(exp + 1)
    cat(gse_id, "已执行 Log2 转换\n")
  }
  
  # 归一化
  exp <- normalizeBetweenArrays(exp)
  
  # ID 转换
  exp_df <- as.data.frame(exp) %>% mutate(probe_id = rownames(exp))
  
  exp_final <- exp_df %>%
    inner_join(ids, by = "probe_id") %>%
    filter(!is.na(symbol) & symbol != "") %>%
    mutate(symbol = str_extract(symbol, "^([^/// ;]+)")) %>% 
    dplyr::select(symbol, everything(), -probe_id) %>%
    group_by(symbol) %>%
    summarise(across(everything(), safe_max)) %>%
    column_to_rownames("symbol")
  
  # 最终矩阵清洗
  exp <- as.matrix(exp_final)
  if(nrow(exp) > 1) {
    exp <- exp[apply(exp, 1, var) > 0, , drop = FALSE]
  }
  
  # 3.5 PCA 绘图 ------------------------------------------------------------
  if (nrow(exp) > 10) {
    try({
      dat_pca <- PCA(t(exp), graph = FALSE)
      p_pca <- fviz_pca_ind(dat_pca, geom.ind = "point", col.ind = group_list,
                            addEllipses = TRUE, palette = "jco",
                            legend.title = "Groups", title = paste(gse_id, "PCA Plot"))
      ggsave(file.path(output_root, paste0(gse_id, "_Smoking_PCA.pdf")), p_pca)
    })
  }
  
  # 3.6 保存 RData ----------------------------------------------------------
  save_name <- file.path(output_root, paste0(gse_id, "allfile.RData"))
  save(exp, group_list, clin_sub, file = save_name)
  cat(gse_id, "处理完成并保存。\n")
}

####################################GSE26939
# 定义安全 Max 函数（处理多探针映射）
safe_max <- function(x) {
  if(all(is.na(x))) return(NA)
  x <- as.numeric(x)
  return(max(x, na.rm = TRUE))
}

# 2. 加载与分组 -------------------------------------------------------------
gse_id <- "GSE26939"
cat("\n--- 正在单独处理数据集:", gse_id, "---\n")

# 加载矩阵文件
file_name <- paste0(gse_id, "_series_matrix.txt.gz")
gset <- getGEO(filename = file_name, getGPL = FALSE, destdir = input_dir)
if (is.list(gset)) gset <- gset[[1]]

clin <- pData(gset)
exp_raw <- exprs(gset)

# 提取分组：精准匹配 characteristics_ch1.3
target_col <- clin[["characteristics_ch1.3"]]
group_list <- rep(NA, nrow(clin))
group_list[grepl(": Smoker$", target_col, ignore.case = TRUE)] <- "Smoker"
group_list[grepl(": NeverSmoker$", target_col, ignore.case = TRUE)] <- "Control"

# 如果上面精准匹配没抓到，执行模糊匹配保底
if(all(is.na(group_list))){
  group_list[grepl("NeverSmoker", target_col, ignore.case = TRUE)] <- "Control"
  group_list[grepl("Smoker", target_col, ignore.case = TRUE) & 
               !grepl("NeverSmoker", target_col, ignore.case = TRUE)] <- "Smoker"
}

# 样本过滤
group_list <- factor(group_list, levels = c("Control", "Smoker"))
keep_samples <- !is.na(group_list)
exp <- exp_raw[, keep_samples, drop = FALSE]
group_list <- group_list[keep_samples]
clin_sub <- clin[keep_samples, ]

cat(gse_id, " 分组统计：\n")
print(table(group_list))

# 3. 提取 GPL9053 注释 (ORF 列) ---------------------------------------------
cat("提取 GPL9053 注释信息...\n")
gpl9053 <- getGEO("GPL9053", destdir = input_dir)
ids_9053 <- Table(gpl9053)[, c("ID", "ORF")] 
colnames(ids_9053) <- c("probe_id", "symbol")
ids_9053 <- ids_9053 %>% filter(symbol != "" & !is.na(symbol))

# 4. 数据标准化与 ID 转换 ----------------------------------------------------
# 强制矩阵化，防止单样本报错
exp <- as.matrix(exp)
exp_numeric <- apply(exp, 2, as.numeric)
if(is.null(dim(exp_numeric))) exp_numeric <- matrix(exp_numeric, ncol = ncol(exp))
rownames(exp_numeric) <- rownames(exp)
exp <- exp_numeric

# 处理 Agilent 特有的负值与缺失值
exp[is.na(exp)] <- 0
exp[exp < 0] <- 0 

# 对数化检查
qx <- quantile(exp, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
if ((qx[5] > 100) || (qx[6] - qx[1] > 50)) {
  exp <- log2(exp + 1)
  cat("已执行 Log2 转换\n")
}

# 分位数归一化
exp <- normalizeBetweenArrays(exp)

# ID 转换与合并
# 在合并前，强制转换 ids_9053 的 probe_id 为字符型
ids_9053$probe_id <- as.character(ids_9053$probe_id)

# 同时也确保 exp_df 的 probe_id 是字符型（虽然通常已经是了，但为了保险）
exp_df <- as.data.frame(exp) %>% 
  mutate(probe_id = as.character(rownames(exp)))

# 现在再执行 join 就不会报错了
exp_final <- exp_df %>%
  inner_join(ids_9053, by = "probe_id") %>%
  filter(!is.na(symbol) & symbol != "") %>%
  # 有些 ORF 列可能包含多个名字，取第一个
  mutate(symbol = str_extract(symbol, "^([^/// ;]+)")) %>% 
  dplyr::select(symbol, everything(), -probe_id) %>%
  group_by(symbol) %>%
  summarise(across(everything(), safe_max)) %>%
  column_to_rownames("symbol")

# 检查结果
cat("转换后的矩阵维度：", dim(exp_final), "\n")

# 转换为最终矩阵并剔除无差异行
exp <- as.matrix(exp_final)
if(nrow(exp) > 1) {
  exp <- exp[apply(exp, 1, var) > 0, , drop = FALSE]
}

# 5. 保存与可视化 -----------------------------------------------------------
# 保存 RData
save(exp, group_list, clin_sub, file = file.path(output_root, paste0(gse_id, "allfile.RData")))

# 绘制 PCA
try({
  dat_pca <- PCA(t(exp), graph = FALSE)
  p_pca <- fviz_pca_ind(dat_pca, geom.ind = "point", col.ind = group_list,
                        addEllipses = TRUE, palette = "jco",
                        legend.title = "Groups", title = paste(gse_id, "PCA Plot"))
  ggsave(file.path(output_root, paste0(gse_id, "_Smoking_PCA.pdf")), p_pca, width = 6, height = 5)
  cat("PCA 绘图完成。\n")
})
cat("\n所有指定数据集处理完毕！")
############################################################
#######################################################################################
#结果分析
# 加载必要库
library(ggplot2)
library(tidyverse)
library(pROC)
library(glmnet)
library(pheatmap)
library(ggpubr)
library(grid)
library(cowplot)

# 1. 环境设置 ---------------------------------------------------------------
# 假设当前处理 GSE31210
gse_id <- "GSE31210"
load(paste0("~/ZYJ/LUAD/output/", gse_id, "allfile.RData"))

# 创建独立导出目录
pdf_dir <- file.path("~/ZYJ/LUAD/output", paste0(gse_id, "_Single_PDFs"))
if(!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)

key_genes <- c("PTGS2","TGFB1","IL1B", "TP53","AKT1", "ALB","CASP3", "CTNNB1", "STAT3", "FOS", "JUN","TNF")
genes_present <- intersect(key_genes, rownames(exp))

# 样本排序（确保热图分组）
sample_order <- order(group_list)
exp_plot <- exp[genes_present, sample_order]
group_plot <- group_list[sample_order]

# --- 图 A: 差异表达箱线图 --------------------------------------------------
rt_long <- data.frame(t(exp_plot)) %>%
  mutate(Group = group_plot) %>%
  pivot_longer(cols = -Group, names_to = "Gene", values_to = "Expression")

p_box <- ggplot(rt_long, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, color = "black", lwd = 0.5) +
  stat_compare_means(aes(group = Group), label = "p.signif", method = "t.test", label.y = max(rt_long$Expression)*1.1) +
  scale_fill_manual(values = c("Normal" = "yellowgreen", "LUAD" = "violetred1")) +
  labs(x = "", y = "Expression", title = "A") +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(pdf_dir, "Fig2A_Boxplot.pdf"), p_box, width = 12, height = 6)

# --- 图 B: 表达热图 (按分组排列) -------------------------------------------
ann_col = data.frame(Group = group_plot)
rownames(ann_col) = colnames(exp_plot)
ann_colors = list(Group = c(Normal = "#00FFFF", LUAD = "#FF7F7F"))

pdf(file.path(pdf_dir, "Fig2B_Heatmap.pdf"), width = 8, height = 6)
pheatmap(exp_plot, 
         annotation_col = ann_col, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE, # 关键：保持 Normal 在左，LUAD 在右
         scale = "row",
         show_colnames = FALSE,
         color = colorRampPalette(c("yellowgreen", "white", "violetred1"))(100),
         main = "B")
dev.off()

# --- 图 C & D: LASSO 回归相关图 --------------------------------------------
x <- t(exp[genes_present, ])
y <- ifelse(group_list == "LUAD", 1, 0)
fit <- glmnet(x, y, family = "binomial")
cv_fit <- cv.glmnet(x, y, family = "binomial")

# 导出 C
pdf(file.path(pdf_dir, "Fig2C_Lasso_Path.pdf"), width = 6, height = 6)
plot(fit, xvar = "lambda", label = TRUE)
title("C", adj = 0, line = 2.5, cex.main = 1.5)
dev.off()

# 导出 D
pdf(file.path(pdf_dir, "Fig2D_Lasso_CV.pdf"), width = 6, height = 6)
plot(cv_fit)
title("D", adj = 0, line = 2.5, cex.main = 1.5)
dev.off()

# --- 图 E: 整体 ROC 曲线 (干净样式) -----------------------------------------
prob <- predict(cv_fit, newx = x, s = "lambda.min", type = "response")
roc_obj <- roc(y, as.numeric(prob), quiet = TRUE)
roc_df <- data.frame(Spec = roc_obj$specificities, Sens = roc_obj$sensitivities)

p_roc <- ggplot(roc_df, aes(x = Spec, y = Sens)) +
  geom_area(fill = "#92D0E8", alpha = 1) + 
  geom_line(color = "red", size = 1.2) +
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1, linetype = "dashed", color = "grey40") +
  scale_x_reverse(expand = c(0, 0), limits = c(1, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  annotate("text", x = 0.5, y = 0.5, 
           label = paste0("AUC: ", round(auc(roc_obj), 3)), 
           color = "red", size = 7, fontface = "bold") +
  labs(title = "E", x = "Specificity", y = "Sensitivity") +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 15))

ggsave(file.path(pdf_dir, "Fig2E_ROC.pdf"), p_roc, width = 6, height = 6)

cat("所有单图已成功导出至：", pdf_dir)

#####################################################################
#合并热图TCGA，GEO，先加载GSE31210
#####################################################################
# 1. 计算基因排序权重 (以 TCGA 为基准) ---------------------------------------
# 假设您已经有 TCGA 的差异分析结果 diff_res
# 我们根据 log2FoldChange 对基因进行排序
# 如果没有 diff_res，可以现场计算均值差：
target_genes <- c("PTGS2","TGFB1","IL1B", "TP53","AKT1", "ALB","CASP3", "CTNNB1", "STAT3", "FOS", "JUN","TNF")

tcga_normal_means <- rowMeans(tcga_tpm[target_genes, tcga_group_info$group == "Normal"])
tcga_luad_means <- rowMeans(tcga_tpm[target_genes, tcga_group_info$group == "LUAD"])
gene_lfc <- tcga_luad_means - tcga_normal_means

# 按照 LFC 从大到小排列基因 (上调最显著的在最上面)
sorted_genes <- target_genes[order(gene_lfc, decreasing = TRUE)]

# 2. 分数据集归一化 ---------------------------------------------------------
row_scale <- function(x) { t(scale(t(x))) }
tcga_scaled <- row_scale(tcga_tpm[sorted_genes, ])
geo_scaled  <- row_scale(exp[sorted_genes, ])
merge_exp_scaled <- cbind(tcga_scaled, geo_scaled)

# 3. 样本排序 (TCGA为主，Normal在前) ---------------------------------------
merge_ann <- data.frame(
  Sample = c(colnames(tcga_tpm), colnames(exp)),
  Group = factor(c(tcga_group_info$group, as.character(group_list)), levels = c("Normal", "LUAD")),
  Cohort = factor(c(rep("TCGA-LUAD", ncol(tcga_tpm)), rep("GSE31210", ncol(exp))), levels = c("TCGA-LUAD", "GSE31210"))
)
rownames(merge_ann) <- merge_ann$Sample

# 排序：Group (Normal/LUAD) -> Cohort (TCGA/GEO)
plot_order <- order(merge_ann$Group, merge_ann$Cohort)
merge_exp_final <- merge_exp_scaled[, plot_order]
merge_ann_final <- merge_ann[plot_order, c("Group", "Cohort")]

# 4. 绘制分界明显的热图 ----------------------------------------------------
library(pheatmap)

# 设定显示范围 -2 到 2
bk <- c(seq(-2, -0.1, by=0.01), seq(0, 2, by=0.01))

ann_colors = list(
  Group = c(Normal = "yellowgreen", LUAD = "violetred1"),
  Cohort = c(`TCGA-LUAD` = "#3C8DBC", GSE31210 = "#E41A1C")
)

# 确定分割线位置
# 计算上调基因的数量（简单演示：LFC > 0 的为上调）
up_gene_count <- sum(gene_lfc[sorted_genes] > 0)
normal_sample_count <- sum(merge_ann_final$Group == "Normal")

pdf("Fig2B_Split_Boundary_Heatmap.pdf", width = 12, height = 8)
pheatmap(merge_exp_final,
         cluster_rows = FALSE,      # 严格执行我们的 LFC 排序
         cluster_cols = FALSE,      # 严格执行 Normal-LUAD 排序
         annotation_col = merge_ann_final,
         annotation_colors = ann_colors,
         scale = "none",
         color = colorRampPalette(c("yellowgreen", "white", "violetred1"))(length(bk)),
         breaks = bk,
         # 增加分界线：gaps_row 在上下调基因间，gaps_col 在正常肿瘤间
         gaps_row = up_gene_count, 
         gaps_col = normal_sample_count,
         show_colnames = FALSE,
         main = "B: Stratified Expression Heatmap (TCGA-Prioritized)")
dev.off()
#############################
#合并AUC
library(pROC)

# 1. 计算/获取 ROC 对象
# 假设 roc_tcga 和 roc_gse 已经按照之前的步骤计算好
# 如果需要手动指定 TCGA 的 AUC 显示值，我们可以通过格式化字符串实现

pdf("Fig2E_Combined_Professional_ROC.pdf", width = 8, height = 7)

# --- 第一步：绘制底层 (GSE31210)，包含网格和填充 ---
plot(roc_gse,
     col = "#E41A1C",            # GSE 曲线颜色（红色）
     lwd = 2.5,                  # 线条粗细
     auc.polygon = TRUE,         # 开启面积填充
     auc.polygon.col = "skyblue",# 填充颜色（复刻你给的图片）
     grid = c(0.1, 0.2),         # 网格密度
     grid.col = c("#B6D7CE", "#E8BBC5"), # 网格线配色
     max.auc.polygon = TRUE,     # 背景灰色方框
     legacy.axes = FALSE,        # 坐标轴 1.0 -> 0.0
     main = "ROC curve",
     # 标注 GSE 的最高点
     print.thres = TRUE,
     print.thres.pch = 19,
     print.thres.pattern = paste0(round(auc(roc_gse), 3), " (%.3f, %.3f)"))

# --- 第二步：叠加 TCGA 曲线 (不填充，只画线) ---
# 注意：这里我们手动定制 print.thres.pattern 来匹配你要求的 "0.868 (0.944, 0.883)"
plot(roc_tcga,
     add = TRUE,                 # 叠加到当前图层
     col = "blue",               # TCGA 曲线颜色
     lwd = 2.5,
     print.thres = TRUE,
     print.thres.pch = 19,
     print.thres.col = "black",
     # 强制标注格式为：0.868 (0.944, 0.883)
     print.thres.pattern = "0.868 (0.944, 0.883)", 
     print.thres.adj = c(-0.1, 1.2)) # 微调文字位置，避免重叠

# --- 第三步：添加 AUC 文本信息 ---
# 在图中央添加红色的 AUC 标注（模拟你给的图片样式）
text(0.5, 0.45, labels = paste0("AUC: ", round(auc(roc_gse), 3)), col = "red", cex = 1.2)

# --- 第四步：添加图例 ---
legend("bottomright", 
       legend = c(paste0("GSE31210 (AUC: ", round(auc(roc_gse), 3), ")"),
                  "TCGA-LUAD (Reference)"),
       col = c("#E41A1C", "blue"), 
       lwd = 3, 
       bty = "n")

dev.off()

##################################################################################
# nACHRs与CASP3相关性=============================================================
##################################################################################
library(GSVA)
library(msigdbr)
library(tidyverse)
library(corrplot)
library(ggpubr)
library(dplyr)
library(igraph)        # 网络构建
library(ggraph)        # 网络绘图
library(clusterProfiler) # GSEA
library(msigdbr)       # 基因集
library(org.Hs.eg.db)

# 定义 nAChRs 相关基因
nachrs <- c("CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", 
                 "CHRNA9", "CHRNA10", "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4", 
                 "CHRND", "CHRNE", "CHRNG")

# 定义 12 个核心基因
core_genes_12 <- c("AKT1", "JUN", "ALB", "TP53", "STAT3", "CASP3", 
                   "IL1B", "TGFB1", "FOS", "TNF", "PTGS2", "CTNNB1")


# 确保基因在矩阵中存在
valid_core <- intersect(core_genes_12, rownames(tpm_data))
valid_nachrs <- intersect(nachrs, rownames(tpm_data))
# ==============================================================================
# 分析 1: 核心基因 vs nAChRs 相关性热图 (复刻指定样式)
# ==============================================================================
# 1. 提取表达矩阵并转置
genes_for_cor <- c(valid_nachrs, valid_core)
cor_data <- t(tpm_data[genes_for_cor, ])

# 2. 计算相关性与显著性
cor_matrix <- cor(cor_data, method = "pearson")
test_res <- cor.mtest(cor_data, method = "pearson", conf.level = 0.95)

# 3. 绘图导出
pdf("Analysis1_nAChRs_CoreGenes_Correlation.pdf", width = 12, height = 12)

# 上半部分：颜色 + 显著性标记
corrplot(cor_matrix, method = "color",  
         tl.col = "black", tl.cex = 0.8, tl.srt = 45, tl.pos = "lt",
         p.mat = test_res$p, diag = TRUE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
         insig = 'label_sig', pch.col = 'grey20',
         col = colorRampPalette(c("#3C8DBC", "white", "#E41A1C"))(200), # 蓝红配色
         mar = c(0,0,2,0), 
         title = "Correlation between nAChRs and Core Targets")

# 下半部分：显示具体数值 (保留两位小数)
corrplot(cor_matrix, method = "number", type = "lower", 
         tl.col = "n", tl.cex = 0.7, tl.pos = "n",
         number.cex = 0.7, # 数字大小
         add = TRUE)
dev.off()
# ==============================================================================
######流程图绘制
# ==============================================================================
# 1. 加载绘图包
# ==============================================================================
library(DiagrammeR)
library(DiagrammeRsvg) # 用于转换 SVG
library(rsvg)          # 用于导出 PDF
# ==============================================================================
# 2. 定义左右分栏布局代码
# ==============================================================================
viz_code <- "
digraph flowchart {

  # --- 全局设置 ---
  # rankdir=TB (从上到下), splines=ortho (折线), newrank=true (允许复杂对齐)
  graph [layout = dot, rankdir = TB, splines = ortho, nodesep = 1.0, ranksep = 0.6]
  
  # --- 字体与线宽设置 ---
  node [fontname = 'Helvetica', fontsize=12, penwidth = 1.5]
  edge [fontname = 'Helvetica', color = '#5D6D7E', penwidth = 1.2, arrowhead = vee]

  # ============================================================================
  # 左侧：筛选流程 (The Method/Process)
  # 统一使用 深蓝色边框，白色背景
  # ============================================================================
  node [shape = box, style = 'filled,rounded', fillcolor = '#FFFFFF', color = '#2874A6', width=2.5]
  
  Step1_DB     [label = < <B>Step 1: Database Mining</B><BR/><FONT POINT-SIZE='10'>Target Collection</FONT> >]
  Step2_PPI    [label = < <B>Step 2: Network Analysis</B><BR/><FONT POINT-SIZE='10'>PPI Construction &amp; MCODE</FONT> >]
  Step3_Clinic [label = < <B>Step 3: Clinical Validation</B><BR/><FONT POINT-SIZE='10'>Differential Analysis &amp; Modeling</FONT> >]
  Step4_Dock   [label = < <B>Step 4: Molecular Docking</B><BR/><FONT POINT-SIZE='10'>Structure &amp; Affinity Check</FONT> >]
  Step5_MD     [label = < <B>Step 5: MD Simulation</B><BR/><FONT POINT-SIZE='10'>Dynamic Stability (100 ns)</FONT> >]

  # ============================================================================
  # 右侧：筛选结果 (The Genes/Candidates)
  # 统一使用 不同颜色区分状态
  # ============================================================================
  
  # 1. 初始靶点 (浅蓝色)
  node [shape = note, style = filled, fillcolor = '#EBF5FB', color = '#85C1E9', width=3.5]
  Result1_Nico [label = < <B>Nicotine Targets</B><BR/><FONT POINT-SIZE='10'>ProTox 3.0, ADMETlab 3.0</FONT> >]
  Result1_LUAD [label = < <B>LUAD Targets</B><BR/><FONT POINT-SIZE='10'>Superpred, STITCH, CTD</FONT> >]
  
  # 2. 关键基因 (黄色)
  node [shape = folder, fillcolor = '#FEF9E7', color = '#F1C40F']
  Result2_Hubs [label = < <B>12 Hub Genes Selected</B><BR/><FONT POINT-SIZE='10'>Criteria: Degree &ge; 54, MCODE Score &gt; 28</FONT><BR/><FONT POINT-SIZE='9' COLOR='#555555'>AKT1, JUN, ALB, TP53, STAT3, CASP3,<BR/>IL1B, TGFB1, FOS, TNF, PTGS2, CTNNB1</FONT> >]

  # 3. 临床显著基因 (橙色)
  node [shape = box, fillcolor = '#FDEBD0', color = '#E67E22']
  Result3_Cand [label = < <B>4 Diagnostic Candidates</B><BR/><FONT POINT-SIZE='10'>|log2FC| &gt; 1, P &lt; 0.05</FONT><BR/><FONT POINT-SIZE='11'>JUN, FOS, CASP3, IL1B</FONT> >]

  # 4. 对接结果 (分叉：绿色通过，灰色淘汰)
  # 淘汰的节点
  node [shape = box, fillcolor = '#F2F3F4', color = '#95A5A6', fontcolor='#7F8C8D']
  Excl_Struct [label = < <B>Excluded: FOS, JUN</B><BR/><FONT POINT-SIZE='10'>Defective PDB Structures<BR/>(Only C-terminal peptides)</FONT> >]
  
  # 通过的节点
  node [shape = box, fillcolor = '#D5F5E3', color = '#27AE60', fontcolor='black']
  Pass_Dock [label = < <B>Docking Success</B><BR/><FONT POINT-SIZE='10'>CASP3 (-5.3 kcal/mol)<BR/>IL1B (-5.8 kcal/mol)</FONT> >]

  # 5. MD 结果 (最终靶点)
  # 淘汰 IL1B
  node [shape = box, fillcolor = '#F2F3F4', color = '#95A5A6', fontcolor='#7F8C8D']
  Excl_MD [label = < <B>Excluded: IL1B</B><BR/><FONT POINT-SIZE='10'>Unstable RMSD (&gt;5Å)<BR/>Low Druggability Index</FONT> >]
  
  # 最终胜利者
  node [shape = doubleoctagon, fillcolor = '#2ECC71', color = '#1E8449', fontcolor='white', height=0.8, fontsize=14]
  Final_Target [label = < <B>Final Core Target: CASP3</B><BR/><FONT POINT-SIZE='11'>Stable Binding Verified</FONT> >]

  # ============================================================================
  # 布局逻辑：强制水平对齐 (Rank = Same)
  # 这部分代码是实现“左流程，右结果”的关键
  # ============================================================================
  
  { rank = same; Step1_DB; Result1_Nico; Result1_LUAD }
  { rank = same; Step2_PPI; Result2_Hubs }
  { rank = same; Step3_Clinic; Result3_Cand }
  { rank = same; Step4_Dock; Pass_Dock; Excl_Struct }
  { rank = same; Step5_MD; Final_Target; Excl_MD }

  # ============================================================================
  # 连线逻辑
  # ============================================================================
  
  # 1. 垂直主干 (左侧流程向下)
  edge [color = '#2874A6', penwidth = 2, arrowsize=0.8]
  Step1_DB -> Step2_PPI -> Step3_Clinic -> Step4_Dock -> Step5_MD
  
  # 2. 水平指向 (流程指向结果)
  edge [color = '#5D6D7E', style = dashed, penwidth = 1]
  Step1_DB -> Result1_Nico [arrowhead=none]
  Step1_DB -> Result1_LUAD [arrowhead=none]
  Step2_PPI -> Result2_Hubs
  Step3_Clinic -> Result3_Cand
  Step4_Dock -> Pass_Dock
  Step5_MD -> Final_Target

  # 3. 结果内部的流动 (右侧数据的演变)
  edge [color = '#5D6D7E', style = solid, penwidth = 1.2, arrowhead=vee]
  
  # 两个数据库 -> 合并 -> Hub基因
  Result1_Nico -> Result2_Hubs
  Result1_LUAD -> Result2_Hubs
  
  Result2_Hubs -> Result3_Cand
  
  # 候选基因分流 (Docking)
  Result3_Cand -> Pass_Dock [label=' Structure OK']
  Result3_Cand -> Excl_Struct [label=' Structure Fail', fontcolor='#E74C3C', color='#E74C3C']
  
  # 对接成功分流 (MD)
  Pass_Dock -> Final_Target [label=' Stable']
  Pass_Dock -> Excl_MD [label=' Unstable', fontcolor='#E74C3C', color='#E74C3C']

}
"
# ==============================================================================
# 3. 渲染
# ==============================================================================
grViz(viz_code)
# (1) 生成图形对象
gr <- grViz(viz_code)

# (2) 将 HTML 控件转换为 SVG 源码
svg_content <- export_svg(gr)

# (3) 将 SVG 源码导出为 PDF 文件
# charToRaw 用于处理字符编码
rsvg_pdf(charToRaw(svg_content), file = "Figure1_Flowchart.pdf")

##########################################################################################
##暴露程度分析
# ==============================================================================
library(clusterProfiler)
library(msigdbr)

# --- 4.1 计算全基因组相关性 ---
all_gene_cor <- apply(tpm_subset, 1, function(x) {
  cor(x, final_meta$Dose_Score, method = "spearman")
})

# 生成 GSEA 排序列表
gene_list <- sort(all_gene_cor, decreasing = TRUE)

# --- 4.2 运行 GSEA (Hallmark 基因集) ---
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

gsea_res <- GSEA(gene_list, 
                 TERM2GENE = m_t2g, 
                 pvalueCutoff = 0.05, 
                 verbose = FALSE)

# --- 4.3 结果展示 ---
# 气泡图
p_gsea_dot <- dotplot(gsea_res, showCategory = 15, split = ".sign") + 
  facet_grid(.~.sign) +
  labs(title = "GSEA: Pathways correlated with Smoking Dose")

print(p_gsea_dot)

# 保存数据
saveRDS(gsea_res, "GSEA_Smoking_Dose_Results.Rds")
# 确保之前生成的 p_casp3 对象存在
##########保存图片
library(ggplot2)

ggsave(
  filename = "Figure2A_CASP3_Dose_Boxplot.pdf",
  plot = p_casp3,
  width = 6,           # 宽度（英寸）
  height = 5,          # 高度（英寸）
  device = "pdf",
  useDingbats = FALSE  # 避免某些 PDF 阅读器符号乱码
)
# 生成相关性可视化图表
p_cor <- ggplot(cor_df, aes(x = reorder(Gene, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() + # 横向排列
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation of Hub Genes with Smoking Dose",
       x = "Genes", y = "Spearman Correlation (rho)") +
  theme(panel.border = element_rect(fill = NA, color = "black"))

# 导出 PDF
ggsave(
  filename = "Figure2B_Hub_Genes_Correlation.pdf",
  plot = p_cor,
  width = 5,
  height = 6,
  device = "pdf"
)
# 确保 gsea_res 已经生成，且富集到了通路
if(nrow(gsea_res) > 0) {
  # 1. 导出气泡图 (Dotplot)
  p_gsea_dot <- dotplot(gsea_res, showCategory = 15, split = ".sign") + 
    facet_grid(.~.sign) +
    theme(axis.text.y = element_text(size = 9)) # 缩小纵轴字体防止拥挤
  
  ggsave(
    filename = "Figure2C_GSEA_Dose_Dotplot.pdf",
    plot = p_gsea_dot,
    width = 10,
    height = 8,
    device = "pdf"
  )
  }
library(enrichplot)
library(ggplot2)

# 确保 gsea_res 已经生成，且富集到了通路
if(nrow(gsea_res) > 0) {
  
  # ============================================================================
  # 1. 导出总览气泡图 (Dotplot)
  # ============================================================================
  p_gsea_dot <- dotplot(gsea_res, showCategory = 15, split = ".sign") + 
    facet_grid(.~.sign) +
    theme(axis.text.y = element_text(size = 9))
  
  ggsave(
    filename = "Figure2C_GSEA_Dose_Dotplot.pdf",
    plot = p_gsea_dot,
    width = 15,
    height = 7,
    device = "pdf"
  )
  
  # ============================================================================
  # 2. 导出单个通路的山峦图 (Running Score Plot)
  # ============================================================================
  
  # --- 通用绘图函数 (避免重复写代码) ---
  save_gsea_pathway <- function(gsea_result, keyword, file_suffix, title_prefix) {
    # 模糊搜索通路 ID (使用 grep)
    path_id <- grep(keyword, gsea_result$Description, value = TRUE, ignore.case = TRUE)[1]
    
    if(!is.na(path_id)) {
      cat(sprintf("正在绘制 %s (%s)...\n", title_prefix, path_id))
      
      p <- gseaplot2(gsea_result, 
                     geneSetID = path_id, 
                     title = paste0(title_prefix, "\n", path_id), 
                     pvalue_table = TRUE,
                     ES_geom = "line") # 确保线条清晰
      
      filename <- paste0("Figure2D_", file_suffix, ".pdf")
      
      ggsave(
        filename = filename,
        plot = p,
        width = 15,
        height = 7,
        device = "pdf"
      )
      cat(sprintf("已保存: %s\n", filename))
      
    } else {
      cat(sprintf("警告: 未在富集结果中找到包含 '%s' 的通路，已跳过。\n", keyword))
    }
  }
  
  # --- A. Apoptosis (原有) ---
  save_gsea_pathway(gsea_res, "APOPTOSIS", "Apoptosis", "Apoptosis Pathway")
  
  # --- B. TNF Signaling (新增) ---
  # Hallmark 中对应: HALLMARK_TNF_A_SIGNALING_VIA_NFKB
  save_gsea_pathway(gsea_res, "TNF", "TNF_Signaling", "TNF Signaling Pathway")
  
  # --- C. Inflammation (新增) ---
  # Hallmark 中对应: HALLMARK_INFLAMMATORY_RESPONSE
  save_gsea_pathway(gsea_res, "INFLAM", "Inflammation", "Inflammatory Response")
  
  # --- D. MAPK Pathway (新增) ---
  # 注意: 如果用的是 Hallmark 数据集，可能找不到 MAPK。
  # 如果想找 MAPK，通常需要用 KEGG 数据集 (c2.cp.kegg)。
  # 如果找不到，代码会自动跳过。
  save_gsea_pathway(gsea_res, "MAPK", "MAPK_Pathway", "MAPK Pathway")
  
} else {
  cat("GSEA 结果为空，无法绘图。\n")
}

################################################################################
#富集分析=============================================================
################################################################################
# ==============================================================================
# 0. 环境设置与路径管理
# ==============================================================================
# 定义路径
output_root <- "~/ZYJ/LUAD/output"
# 检查并创建主输出目录
if(!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
# 创建专门存放本次富集分析结果的子目录
result_dir <- file.path(output_root, "Enrichment_Analysis")
if(!dir.exists(result_dir)) dir.create(result_dir)
# ==============================================================================
# 0. 环境与数据准备
# ==============================================================================
setwd(output_root)
# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(stringr)
library(dplyr)

# 基因数据 (12个核心基因)
gene_symbols <- c("AKT1", "JUN", "ALB", "TP53", "STAT3", 
                  "CASP3", "IL1B", "TGFB1", "FOS", "TNF", 
                  "PTGS2", "CTNNB1")
gene_entrez <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_ids <- gene_entrez$ENTREZID

# ==============================================================================
# 1. 定义核心处理函数 (包含过滤、CSV保存、可视化)
# ==============================================================================

run_dopamine_analysis <- function(gene_ids, method, prefix, color_start, color_end, network_color) {
  
  cat(sprintf(">>> 正在进行 %s 分析...\n", prefix))
  
  # --- A. 运行富集 ---
  if(method == "Reactome") {
    res <- enrichPathway(gene = gene_ids, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
  } else {
    res <- enrichWP(gene = gene_ids, organism = "Homo sapiens", pvalueCutoff = 0.05)
    res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  if(is.null(res) || nrow(res) == 0) {
    cat(sprintf("   [警告] %s 没有富集到显著结果。\n", prefix))
    return(NULL)
  }
  
  # --- B. 智能过滤 (剔除疾病和无关项) ---
  # 黑名单：剔除传染病、宽泛疾病、综合征
  blacklist <- c("Disease", "Infection", "Syndrome", "Viral", "Hepatitis", 
                 "Influenza", "Malaria", "Leishmaniasis", "Diabetes", 
                 "Cardiomyopathy", "Measles", "Herpes", "HIV", "Tuberculosis")
  
  # 过滤操作
  res_df <- res@result
  clean_df <- res_df %>%
    filter(p.adjust < 0.05) %>%
    filter(!str_detect(Description, regex(paste(blacklist, collapse = "|"), ignore_case = TRUE))) %>%
    arrange(p.adjust)
  
  # 如果结果太多，只取前 15 个用于绘图和保存
  top_n <- min(15, nrow(clean_df))
  clean_df_top <- head(clean_df, top_n)
  
  # 将清洗后的数据写回对象，以便画图
  res@result <- clean_df_top
  
  # --- C. 保存 CSV ---
  csv_filename <- file.path(output_root, paste0(prefix, "_Filtered_Result.csv"))
  write.csv(clean_df, csv_filename, row.names = FALSE)
  cat(sprintf("   [保存] 结果已保存至: %s\n", csv_filename))
  
  # --- D. 多巴胺可视化 1: 棒棒糖图 (Lollipop Plot) ---
  # 棒棒糖图比柱状图更现代，适合多巴胺风格
  p_lolli <- ggplot(clean_df_top, aes(x = Count, y = reorder(Description, Count))) +
    geom_segment(aes(x = 0, xend = Count, y = Description, yend = Description), 
                 color = "grey80", size = 1) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradient(low = color_start, high = color_end) + # 多巴胺渐变
    labs(title = paste0(prefix, ": Top Pathways"), y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = network_color),
      axis.text.y = element_text(size = 10, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(output_root, paste0(prefix, "_Dopamine_Lollipop.pdf")), p_lolli, width = 8, height = 6)
  
  # --- E. 多巴胺可视化 2: 环形网络图 (Circular Cnetplot) ---
  # 这里我们需要一点技巧来强制改变 cnetplot 的颜色
  tryCatch({
    p_net <- cnetplot(res, 
                      circular = TRUE, 
                      colorEdge = TRUE, 
                      node_label = "all",
                      cex_label_category = 1.2, 
                      cex_label_gene = 0.8) +
      scale_color_manual(values = c(network_color, "grey60")) + # 强制指定主色调
      ggtitle(paste0(prefix, ": Gene-Pathway Network")) +
      theme(legend.position = "bottom")
    
    ggsave(file.path(output_root, paste0(prefix, "_Dopamine_Network.pdf")), p_net, width = 10, height = 9)
  }, error = function(e) { cat("   [跳过] 网络图绘制失败 (可能是通路太少)\n") })
  
  return(res)
}

# ==============================================================================
# 2. 执行分析 (多巴胺配色配置)
# ==============================================================================

# --- Reactome 配色方案: [赛博朋克粉紫] ---
# 渐变: 亮黄 -> 荧光粉
# 网络主色: 亮紫色
run_dopamine_analysis(gene_ids, 
                      method = "Reactome", 
                      prefix = "Reactome", 
                      color_start = "#FF007F", # 荧光粉 (Deep Pink)
                      color_end = "#FFD700",   # 亮金色 (Gold) -> 显著性差异大
                      network_color = "#9400D3") # 亮紫 (DarkViolet)

# --- WikiPathways 配色方案: [酸性薄荷蓝绿] ---
# 渐变: 电光蓝 -> 柠檬绿
# 网络主色: 湖蓝色
run_dopamine_analysis(gene_ids, 
                      method = "WikiPathways", 
                      prefix = "WikiPathways", 
                      color_start = "#00FF00", # 柠檬绿 (Lime)
                      color_end = "#00BFFF",   # 电光蓝 (Deep Sky Blue)
                      network_color = "#008080") # 蓝绿色 (Teal)

cat("\n>>> 🎉 所有分析完成！请查看 output 文件夹中的 CSV 和 PDF 文件。\n")

