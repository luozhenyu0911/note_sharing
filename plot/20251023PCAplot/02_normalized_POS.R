rm(list = ls())
setwd('E:\\09_metabolome\\02.cleandata')

# 加载必要包
library(ggplot2)
library(ggrepel)  # 用于避免标签重叠
library(readxl)
library(DESeq2)
library(openxlsx)
library(dplyr)
library(limma)
library(sva)
library(tidyr)

data <- read_excel('../raw_data_from_xueling/BIGCS_dr2505/LCQEPOSfiltered.dropna.xlsx')
rowname_data <- data$`MS2 name`
data <- data[, 2:ncol(data)]
rownames(data)<- rowname_data

random_sample_distribution <- function(data, number){
  set.seed(123) # 设置随机种子保证可重复性
  selected_cols <- sample(colnames(data), number) # 随机选3列
  df_long <- data %>%
  # 保留化学物质名称（行名转为列）
  tibble::rownames_to_column("Compound") %>%
  # 只选择需要的列
  select(Compound, all_of(selected_cols)) %>%
  # 转换为长格式
  pivot_longer(
    cols = -Compound,
    names_to = "Sample",
    values_to = "Value"
  )

  # 绘制密度图
  ggplot(df_long, aes(x = Value, fill = Sample)) +
    geom_density(alpha = 0.5) + 
    labs(title = "Density Plot of Randomly Selected Columns",
         x = "Measurement Value", 
         y = "Density") +
    theme_minimal()
}

random_sample_distribution(data, 10)

# data=data.frame(normalizeBetweenArrays(log(data+1)))
data=log(data+1)
random_sample_distribution(data, 10)

### loading pheno data
coldata <- read.table('../01.sampleinfo/dr2505_idmatch_en_addQC_year.txt', sep = '\t', header = 1)
colnames(coldata)
table(coldata$type_year)

PCA_plot <- function(data, dataname, coldata, figname, figname2){
  data <- na.omit(as.data.frame(data) ) 
  write.xlsx(data, dataname,rowNames = TRUE)
  pca_result <- prcomp(t(data), scale= T)  # scale.=TRUE确保标准化
  summary(pca_result)  # 查看主成分方差贡献率
  
  # 提取主成分得分和方差贡献率
  pca_scores <- as.data.frame(pca_result$x)
  percentage <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)  # 计算各PC方差百分比
  
  Year = as.factor(coldata$Year)
  Member = as.factor(coldata$Sample_type)
  
  p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Year,shape = Member, label = rownames(pca_scores))) +
    geom_point(size = 1) +
    stat_ellipse(aes(fill = Year), geom = "polygon", alpha = 0.2, level = 0.95) +
    # 添加参考线
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # Y=0的线
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # X=0的线
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    labs(
      # title = "PCA Plot with Confidence Ellipses",
      x = paste0("PC1 (", percentage[1], "%)"),
      y = paste0("PC2 (", percentage[2], "%)")
    ) +
    theme_minimal()+
    guides(
      color = guide_legend(override.aes = list(size = 3)),  # 颜色图例符号大小
      shape = guide_legend(override.aes = list(size = 3))   # 形状图例符号大小
    ) 
  p2 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Member,shape = Member, label = rownames(pca_scores))) +
    geom_point(size = 1) +
    stat_ellipse(aes(fill = Member), geom = "polygon", alpha = 0.2, level = 0.95) +
    # 添加参考线
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # Y=0的线
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # X=0的线
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    labs(
      # title = "PCA Plot with Confidence Ellipses",
      x = paste0("PC1 (", percentage[1], "%)"),
      y = paste0("PC2 (", percentage[2], "%)")
    ) +
    theme_minimal()+
    guides(
      color = guide_legend(override.aes = list(size = 3)),  # 颜色图例符号大小
      shape = guide_legend(override.aes = list(size = 3))   # 形状图例符号大小
    ) 
  ggsave(figname, plot = p, width = 6, height = 4, units = "in")
  ggsave(figname2, plot = p2, width = 6, height = 4, units = "in")
  
}


#  13  14  15  16  17  QC

#################  raw pca plot  ###########################
pca_scores <- PCA_plot(data, 'LCQEPOSfiltered.dropna.xlsx', coldata, 
                       'LCQEPOSfiltered.dropna.PCA_year.pdf', 'LCQEPOSfiltered.dropna.PCA_member.pdf')

###################  combat 去批次 #############################
# # 设置生物学分类，告诉函数不要把生物学差异整没了 
coldata$Sample_type <- factor(coldata$Sample_type, levels = c("Cord blood", "Maternal blood", "QC"))
batch <- coldata$Year # 批次

mod <- model.matrix(~as.factor(Sample_type), data=coldata) 

expr_combat <- ComBat(dat = data, batch = batch)# , mod=mod

PCA_plot(expr_combat, "LCQEPOSfiltered.dropna.combat.xlsx", coldata , 
         'LCQEPOSfiltered.dropna.combat.PCA_year.pdf', 'LCQEPOSfiltered.dropna.combat.PCA_member.pdf')



#############  2. 组织类型分离，再根据年份进行提取子集，再对子集的每一行进行Z-score

# "Maternal blood"
table(coldata$Sample_type)
MB_pheno <- coldata %>% filter(Sample_type =="Maternal blood")
MB_data <- subset(data,  select = c(MB_pheno$Meta_ID))

#################  raw pca plot  ###########################
PCA_plot(MB_data, 'LCQEPOSfiltered.dropna.MB.xlsx', MB_pheno, 
         'LCQEPOSfiltered.dropna.MB.PCA_year.pdf', 'LCQEPOSfiltered.dropna.MB.PCA_member.pdf')

###################  combat 去批次 #############################
table(MB_pheno$Year)
# # 设置生物学分类，告诉函数不要把生物学差异整没了 
batch <- MB_pheno$Year # 批次

# mod <- model.matrix(~as.factor(Sample_type), data=MB_pheno) # mod=mod

MB_expr_combat <- as.data.frame( ComBat(dat = MB_data, batch = batch))

PCA_plot(MB_expr_combat, "LCQEPOSfiltered.dropna.MB.combat.xlsx", MB_pheno , 
         'LCQEPOSfiltered.dropna.MB.combat.PCA_year.pdf', 'LCQEPOSfiltered.dropna.MB.combat.PCA_member.pdf')

# "Cord blood"
table(coldata$Sample_type)
CB_pheno <- coldata %>% filter(Sample_type =="Cord blood")
CB_data <- subset(data,  select = c(CB_pheno$Meta_ID))

#################  raw pca plot  ###########################
PCA_plot(CB_data, 'LCQEPOSfiltered.dropna.CB.xlsx', CB_pheno, 
         'LCQEPOSfiltered.dropna.CB.PCA_year.pdf', 'LCQEPOSfiltered.dropna.CB.PCA_member.pdf')

###################  combat 去批次 #############################
table(CB_pheno$Year)
# # 设置生物学分类，告诉函数不要把生物学差异整没了 
batch <- CB_pheno$Year # 批次

mod <- model.matrix(~as.factor(Sample_type), data=CB_pheno) # mod=mod

CB_expr_combat <- as.data.frame( ComBat(dat = CB_data, batch = batch))

PCA_plot(CB_expr_combat, "LCQEPOSfiltered.dropna.CB.combat.xlsx", CB_pheno , 
         'LCQEPOSfiltered.dropna.CB.combat.PCA_year.pdf', 'LCQEPOSfiltered.dropna.CB.combat.PCA_member.pdf')


########################################   NEG   ###############################################
data <- read_excel('../raw_data_from_xueling/BIGCS_dr2505/LCQENEGfiltered.dropna.xlsx')
rowname_data <- data$`MS2 name`
data <- data[, 2:ncol(data)]
rownames(data)<- rowname_data

random_sample_distribution <- function(data, number){
  set.seed(123) # 设置随机种子保证可重复性
  selected_cols <- sample(colnames(data), number) # 随机选3列
  df_long <- data %>%
    # 保留化学物质名称（行名转为列）
    tibble::rownames_to_column("Compound") %>%
    # 只选择需要的列
    select(Compound, all_of(selected_cols)) %>%
    # 转换为长格式
    pivot_longer(
      cols = -Compound,
      names_to = "Sample",
      values_to = "Value"
    )
  
  # 绘制密度图
  ggplot(df_long, aes(x = Value, fill = Sample)) +
    geom_density(alpha = 0.5) + 
    labs(title = "Density Plot of Randomly Selected Columns",
         x = "Measurement Value", 
         y = "Density") +
    theme_minimal()
}

random_sample_distribution(data, 10)

# data=data.frame(normalizeBetweenArrays(log(data+1)))
data=log(data+1)
random_sample_distribution(data, 10)

### loading pheno data
coldata <- read.table('../01.sampleinfo/dr2505_idmatch_en_addQC_year.txt', sep = '\t', header = 1)
colnames(coldata)
table(coldata$type_year)

PCA_plot <- function(data, dataname, coldata, figname, figname2){
  data <- na.omit(as.data.frame(data) ) 
  write.xlsx(data, dataname,rowNames = TRUE)
  pca_result <- prcomp(t(data), scale= T)  # scale.=TRUE确保标准化
  summary(pca_result)  # 查看主成分方差贡献率
  
  # 提取主成分得分和方差贡献率
  pca_scores <- as.data.frame(pca_result$x)
  percentage <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)  # 计算各PC方差百分比
  
  Year = as.factor(coldata$Year)
  Member = as.factor(coldata$Sample_type)
  
  p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Year,shape = Member, label = rownames(pca_scores))) +
    geom_point(size = 1) +
    stat_ellipse(aes(fill = Year), geom = "polygon", alpha = 0.2, level = 0.95) +
    # 添加参考线
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # Y=0的线
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # X=0的线
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    labs(
      # title = "PCA Plot with Confidence Ellipses",
      x = paste0("PC1 (", percentage[1], "%)"),
      y = paste0("PC2 (", percentage[2], "%)")
    ) +
    theme_minimal()+
    guides(
      color = guide_legend(override.aes = list(size = 3)),  # 颜色图例符号大小
      shape = guide_legend(override.aes = list(size = 3))   # 形状图例符号大小
    ) 
  p2 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Member,shape = Member, label = rownames(pca_scores))) +
    geom_point(size = 1) +
    stat_ellipse(aes(fill = Member), geom = "polygon", alpha = 0.2, level = 0.95) +
    # 添加参考线
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # Y=0的线
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +  # X=0的线
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#CC79A7","#CC0033", "#4A4D9A")) +
    labs(
      # title = "PCA Plot with Confidence Ellipses",
      x = paste0("PC1 (", percentage[1], "%)"),
      y = paste0("PC2 (", percentage[2], "%)")
    ) +
    theme_minimal()+
    guides(
      color = guide_legend(override.aes = list(size = 3)),  # 颜色图例符号大小
      shape = guide_legend(override.aes = list(size = 3))   # 形状图例符号大小
    ) 
  ggsave(figname, plot = p, width = 6, height = 4, units = "in")
  ggsave(figname2, plot = p2, width = 6, height = 4, units = "in")
  
}


#  13  14  15  16  17  QC

#################  raw pca plot  ###########################
PCA_plot(data, 'LCQENEGfiltered.dropna.xlsx', coldata, 
         'LCQENEGfiltered.dropna.PCA_year.pdf', 'LCQENEGfiltered.dropna.PCA_member.pdf')

###################  combat 去批次 #############################
# # 设置生物学分类，告诉函数不要把生物学差异整没了 
coldata$Sample_type <- factor(coldata$Sample_type, levels = c("Cord blood", "Maternal blood", "QC"))
batch <- coldata$Year # 批次

mod <- model.matrix(~as.factor(Sample_type), data=coldata) 

expr_combat <- ComBat(dat = data, batch = batch)# , mod=mod

PCA_plot(expr_combat, "LCQENEGfiltered.dropna.combat.xlsx", coldata , 
         'LCQENEGfiltered.dropna.combat.PCA_year.pdf', 'LCQENEGfiltered.dropna.combat.PCA_member.pdf')



#############  2. 组织类型分离，再根据年份进行提取子集，再对子集的每一行进行Z-score

# "Maternal blood"
table(coldata$Sample_type)
MB_pheno <- coldata %>% filter(Sample_type =="Maternal blood")
MB_data <- subset(data,  select = c(MB_pheno$Meta_ID))

#################  raw pca plot  ###########################
PCA_plot(MB_data, 'LCQENEGfiltered.dropna.MB.xlsx', MB_pheno, 
         'LCQENEGfiltered.dropna.MB.PCA_year.pdf', 'LCQENEGfiltered.dropna.MB.PCA_member.pdf')

###################  combat 去批次 #############################
table(MB_pheno$Year)
# # 设置生物学分类，告诉函数不要把生物学差异整没了 
batch <- MB_pheno$Year # 批次

# mod <- model.matrix(~as.factor(Sample_type), data=MB_pheno) # mod=mod

MB_expr_combat <- as.data.frame( ComBat(dat = MB_data, batch = batch))

PCA_plot(MB_expr_combat, "LCQENEGfiltered.dropna.MB.combat.xlsx", MB_pheno , 
         'LCQENEGfiltered.dropna.MB.combat.PCA_year.pdf', 'LCQENEGfiltered.dropna.MB.combat.PCA_member.pdf')

# "Cord blood"
table(coldata$Sample_type)
CB_pheno <- coldata %>% filter(Sample_type =="Cord blood")
CB_data <- subset(data,  select = c(CB_pheno$Meta_ID))

#################  raw pca plot  ###########################
PCA_plot(CB_data, 'LCQENEGfiltered.dropna.CB.xlsx', CB_pheno, 
         'LCQENEGfiltered.dropna.CB.PCA_year.pdf', 'LCQENEGfiltered.dropna.CB.PCA_member.pdf')

###################  combat 去批次 #############################
table(CB_pheno$Year)
# # 设置生物学分类，告诉函数不要把生物学差异整没了 
batch <- CB_pheno$Year # 批次

mod <- model.matrix(~as.factor(Sample_type), data=CB_pheno) # mod=mod

CB_expr_combat <- as.data.frame( ComBat(dat = CB_data, batch = batch))

PCA_plot(CB_expr_combat, "LCQENEGfiltered.dropna.CB.combat.xlsx", CB_pheno , 
         'LCQENEGfiltered.dropna.CB.combat.PCA_year.pdf', 'LCQENEGfiltered.dropna.CB.combat.PCA_member.pdf')








