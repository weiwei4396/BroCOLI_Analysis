# 2025.11.21
# Review and organize the SIRV assessment code;
# Author: wei pan;

##########################

# 这一部分是SIRV的定量; ONT cDNA

##########################

# 1.读取标准答案;
library(readxl)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(reshape2)
Answer_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/SIRV_Set1_True.xlsx"
Answer <- read_excel(Answer_path, sheet = "Sheet2", col_names = FALSE)
colnames(Answer) <- c("Transcript", "Truth")
Answer$Truth_xlabel <- 0
for ( i in c(1:dim(Answer)[1]) ) {
  if (Answer$Truth[i] == 0.03125) {
    Answer$Truth_xlabel[i] = "1/32"
  } else if (Answer$Truth[i] == 0.25) {
    Answer$Truth_xlabel[i] = "1/4"
  } else if (Answer$Truth[i] == 0.5) {
    Answer$Truth_xlabel[i] = "1/2"
  } else if (Answer$Truth[i] == 1) {
    Answer$Truth_xlabel[i] = "1"
  } else {
    Answer$Truth_xlabel[i] = "4"
  }
}
Answer$Truth_xlabel <- factor(Answer$Truth_xlabel, levels=c("1/32","1/4","1/2","1","4"))
rm(Answer_path, i)


# 2.获得每个转录本的长度并合并;
# 2.1 获得每个转录本的长度;
sirv_gtf_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/SIRV_100.gtf"
sirv_gtf <- read.delim(sirv_gtf_path, header=FALSE, comment.char="#")
sirv_gtf <- tidyr::separate(sirv_gtf,V9,into = c("b","gene_ID","c","transcript_ID","d","exon_number"),sep = ' ')
sirv_gtf <- tidyr::separate(sirv_gtf,transcript_ID,into = c("transcript_ID","c"),sep = ';')
sirv_gtf <- tidyr::separate(sirv_gtf,gene_ID,into = c("gene_ID","d"),sep = ';')
sirv_gtf <- subset(sirv_gtf,select=c("V1","V2","V3","V4","V5","V7","transcript_ID","gene_ID"))
colnames(sirv_gtf) <-c("chr","seqname","type","start","end","strand","transcript_ID","gene_ID")
rm(sirv_gtf_path)
sirv_transcript_name <- unique(sirv_gtf$transcript_ID)
sirv_transcript_length <- c()
for (namei in sirv_transcript_name) {
  temp_sirv <- subset(sirv_gtf, transcript_ID == namei)
  lengthAll = 0
  for ( i in c(1:dim(temp_sirv)[1]) ) {
    lengthAll = lengthAll + (temp_sirv[i,5] - temp_sirv[i,4] + 1)
  }
  sirv_transcript_length <- c(sirv_transcript_length, lengthAll)
}
sirv2length <- data.frame(Transcript=sirv_transcript_name, Length=sirv_transcript_length)
rm(sirv_transcript_name, sirv_transcript_length, temp_sirv, lengthAll, namei, i)
rm(sirv_gtf)

# 2.2 合并长度和答案;
Answer <- merge(Answer, sirv2length, by="Transcript")
Answer$len_size <- 0
for (i in c(1:dim(Answer)[1])) {
  if (Answer$Length[i] < 500) {
    Answer$len_size[i] <- 1
  } else if (Answer$Length[i] < 1000) {
    Answer$len_size[i] <- 2
  } else if (Answer$Length[i] < 1500) {
    Answer$len_size[i] <- 3
  } else if (Answer$Length[i] < 2000) {
    Answer$len_size[i] <- 4
  } else if (Answer$Length[i] < 2500) {
    Answer$len_size[i] <- 5
  }
}
Answer$len_size <- factor(Answer$len_size, levels = c(1,2,3,4,5))
rm(sirv2length, i)


# 3.读取每个软件;
# 3.1 BroCOLI
# ONT cDNA;
BroCOLI_transcript_counts_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/Quantification/BroCOLI/counts_transcript.txt"
BroCOLI_predict <- read.table(file=BroCOLI_transcript_counts_path, header = TRUE)
rm(BroCOLI_transcript_counts_path)
colnames(BroCOLI_predict) <- c("Transcript", "Gene", "Predict")
data_plot <- merge(Answer, BroCOLI_predict, by=c("Transcript"))
data_plot$Predict[is.na(data_plot$Predict)] <- 0
data_plot$Predict <- data_plot$Predict + 0.1
pcc <- cor(data_plot$Truth, data_plot$Predict, method = "pearson")
pcc
scc <- cor(data_plot$Truth, data_plot$Predict, method = "spearman")
scc
data_plot$Truth <- as.factor(data_plot$Truth)

pdf("C:/Users/Lenovo/Desktop/BroCOLI_Quant.pdf", width=5, height=5)
ggplot(data_plot,aes(x=Truth_xlabel, y=Predict)) + 
  geom_boxplot(size=0.5, outlier.fill="white", outlier.shape=NA, width = 0.4) + 
  geom_jitter(aes(fill=Truth, color=Truth, size=len_size), width = 0, shape = 21, alpha=0.75)+
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_blank(),
        text = element_text(size = 21),
        plot.title = element_text(size=24,hjust=0.5),
        legend.position = "bottom", 
        legend.text = element_text(size = 24))+
  labs(x="SIRV concentration (fmol/ul)", y="Isoform abundance (read count)", title="BroCOLI") +
  scale_y_log10(
    breaks = c(10^(-1), 10^0, 10^1, 10^2, 10^3, 10^4),
    labels = c("0", expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
    limits = c(0.1,10^4)) +
  scale_size_manual(values = c(2, 3, 4, 5, 6), 
                    name = "Isoform length (nt)", 
                    labels = c("500", "1000", "1500", "2000", "2500") ) +
  guides(fill=FALSE, color=FALSE) +
  annotate("text", label = bquote("Pearson'" * r * " = " * .(round(pcc, 3))), x="1/32", y=10000, size=6, hjust = 0) +
  annotate("text", label = bquote("Spearman'" * rho * " = " * .(round(scc, 3))), x="1/32", y=3000, size=6, hjust = 0)
dev.off()

# 以上是定量图;
# 下面是计算Annotation百分百时候的sensitivity和precision
# 预测的\真实的   Positive   Negative   真实的
#      Positive     TP         FP
#      Negative     FN         TN
# 首先将预测的大于0的保留下来;
BroCOLI_predict <- BroCOLI_predict %>% filter(Predict > 0)
TP_FP <- length(BroCOLI_predict$Transcript)
TP_FN <- length(Answer$Transcript)
TP <- length(intersect(BroCOLI_predict$Transcript, Answer$Transcript))
precision <- TP/TP_FP
precision
sensitivity <- TP/TP_FN
sensitivity


# 3.2 Bambu 最佳NDR = 0.039
# ONT cDNA;
Bambu_transcript_counts_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/Quantification/Bambu/counts_transcript.txt"
Bambu_predict <- read.table(file=Bambu_transcript_counts_path, header = TRUE)
colnames(Bambu_predict) <- c("Transcript", "Gene", "Predict")
data_plot <- merge(Answer, Bambu_predict, by=c("Transcript"), all.x = TRUE)
data_plot$Predict[is.na(data_plot$Predict)] <- 0
data_plot$Predict <- data_plot$Predict + 0.1

pcc <- cor(data_plot$Truth, data_plot$Predict, method = "pearson")
pcc
scc <- cor(data_plot$Truth, data_plot$Predict, method = "spearman")
scc
data_plot$Truth <- as.factor(data_plot$Truth)


pdf("C:/Users/Lenovo/Desktop/Bambu_Quant.pdf", width=5, height=5)
ggplot(data_plot,aes(x=Truth_xlabel, y=Predict)) + 
  geom_boxplot(size=0.5, outlier.fill="white", outlier.shape=NA, width = 0.4) + 
  geom_jitter(aes(fill=Truth, color=Truth, size=len_size), width = 0, shape = 21, alpha=0.75)+
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_blank(),
        text = element_text(size = 24),
        plot.title = element_text(size=24,hjust=0.5),
        legend.position = "bottom", 
        legend.text = element_text(size = 24))+
  labs(x="SIRV concentration (fmol/ul)", y="Isoform abundance (read count)", title="Bambu") +
  scale_y_log10(
    breaks = c(10^(-1), 10^0, 10^1, 10^2, 10^3, 10^4),
    labels = c("0", expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
    limits = c(0.1,10^4)) +
  scale_size_manual(values = c(2, 3, 4, 5, 6), 
                    name = "Isoform length (nt)", 
                    labels = c("500", "1000", "1500", "2000", "2500") ) +
  guides(fill=FALSE, color=FALSE) +
  annotate("text", label = bquote("Pearson'" * r * " = " * .(round(pcc, 3))), x="1/32", y=10000, size=6, hjust = 0) +
  annotate("text", label = bquote("Spearman'" * rho * " = " * .(round(scc, 3))), x="1/32", y=3000, size=6, hjust = 0)
dev.off()
# 以上是定量图;
# 下面是计算Annotation百分百时候的sensitivity和precision
# 预测的\真实的   Positive   Negative   真实的
#      Positive     TP         FP
#      Negative     FN         TN
# 首先将预测的大于0的保留下来;
Bambu_predict <- Bambu_predict %>% filter(Predict > 0)
TP_FP <- length(Bambu_predict$Transcript)
TP_FN <- length(Answer$Transcript)
TP <- length(intersect(Bambu_predict$Transcript, Answer$Transcript))
precision <- TP/TP_FP
precision
sensitivity <- TP/TP_FN
sensitivity



# 3.3 Isoquant
# ONT cDNA;
# only annotated 
# Isoquant_transcript_counts_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/Quantification/Isoquant/OUT/OUT.transcript_counts.tsv"
# model
Isoquant_transcript_counts_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/Quantification/Isoquant/OUT/OUT.transcript_model_counts.tsv"
Isoquant_predict <- read.table(file=Isoquant_transcript_counts_path, header = FALSE)
colnames(Isoquant_predict) <- c("Transcript", "Predict")
data_plot <- merge(Answer, Isoquant_predict, by=c("Transcript"), all.x = TRUE)
data_plot$Predict[is.na(data_plot$Predict)] <- 0
data_plot$Predict <- data_plot$Predict + 0.1

pcc <- cor(data_plot$Truth, data_plot$Predict, method = "pearson")
pcc
scc <- cor(data_plot$Truth, data_plot$Predict, method = "spearman")
scc
data_plot$Truth <- as.factor(data_plot$Truth)

pdf("C:/Users/Lenovo/Desktop/IsoQuant_Quant.pdf", width=5, height=5)
ggplot(data_plot,aes(x=Truth_xlabel, y=Predict)) + 
  geom_boxplot(size=0.5, outlier.fill="white", outlier.shape=NA, width = 0.4) + 
  geom_jitter(aes(fill=Truth, color=Truth, size=len_size), width = 0, shape = 21, alpha=0.75)+
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_blank(),
        text = element_text(size = 24),
        plot.title = element_text(size=24,hjust=0.5),
        legend.position = "bottom", 
        legend.text = element_text(size = 24))+
  labs(x="SIRV concentration (fmol/ul)", y="Isoform abundance (read count)", title="IsoQuant") +
  scale_y_log10(
    breaks = c(10^(-1), 10^0, 10^1, 10^2, 10^3, 10^4),
    labels = c("0", expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
    limits = c(0.1,10^4)) +
  scale_size_manual(values = c(2, 3, 4, 5, 6), 
                    name = "Isoform length (nt)", 
                    labels = c("500", "1000", "1500", "2000", "2500") ) +
  guides(fill=FALSE, color=FALSE) +
  annotate("text", label = bquote("Pearson'" * r * " = " * .(round(pcc, 3))), x="1/32", y=10000, size=6, hjust = 0) +
  annotate("text", label = bquote("Spearman'" * rho * " = " * .(round(scc, 3))), x="1/32", y=3000, size=6, hjust = 0)
dev.off()
# 以上是定量图;
# 下面是计算Annotation百分百时候的sensitivity和precision
# 预测的\真实的   Positive   Negative   真实的
#      Positive     TP         FP
#      Negative     FN         TN
# 首先将预测的大于0的保留下来;
# Isoquant比别人特殊, 需要先删除67到69行, 可以仔细看一下;
Isoquant_predict <- Isoquant_predict[-((nrow(Isoquant_predict)-2):nrow(Isoquant_predict)), ]
Isoquant_predict <- Isoquant_predict %>% filter(Predict > 0)
TP_FP <- length(Isoquant_predict$Transcript)
TP_FN <- length(Answer$Transcript)
TP <- length(intersect(Isoquant_predict$Transcript, Answer$Transcript))
precision <- TP/TP_FP
precision
sensitivity <- TP/TP_FN
sensitivity


# 3.4 Isosceles
# ONT cDNA;
# strict
Isosceles_transcript_counts_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/Quantification/Isosceles/counts_transcript.txt"
# de novo loose
# Isosceles_transcript_counts_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/Quantification/Isosceles/counts_transcript_de_nove_loose.txt"
Isosceles_predict <- read.table(file=Isosceles_transcript_counts_path, header = TRUE)
colnames(Isosceles_predict) <- c("Others", "Predict", "Transcript")
data_plot <- merge(Answer, Isosceles_predict, by=c("Transcript"), all.x = TRUE)
data_plot$Predict[is.na(data_plot$Predict)] <- 0
data_plot$Predict <- data_plot$Predict + 0.1

pcc <- cor(data_plot$Truth, data_plot$Predict, method = "pearson")
pcc
scc <- cor(data_plot$Truth, data_plot$Predict, method = "spearman")
scc
data_plot$Truth <- as.factor(data_plot$Truth)

pdf("C:/Users/Lenovo/Desktop/Isosceles_Quant.pdf", width=5, height=5)
ggplot(data_plot,aes(x=Truth_xlabel, y=Predict)) + 
  geom_boxplot(size=0.5, outlier.fill="white", outlier.shape=NA, width = 0.4) + 
  geom_jitter(aes(fill=Truth, color=Truth, size=len_size), width = 0, shape = 21, alpha=0.75)+
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_blank(),
        text = element_text(size = 24),
        plot.title = element_text(size=24,hjust=0.5),
        legend.position = "bottom", 
        legend.text = element_text(size = 24))+
  labs(x="SIRV concentration (fmol/ul)", y="Isoform abundance (read count)", title="Isosceles") +
  scale_y_log10(
    breaks = c(10^(-1), 10^0, 10^1, 10^2, 10^3, 10^4),
    labels = c("0", expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
    limits = c(0.1,10^4)) +
  scale_size_manual(values = c(2, 3, 4, 5, 6), 
                    name = "Isoform length (nt)", 
                    labels = c("500", "1000", "1500", "2000", "2500") ) +
  guides(fill=FALSE, color=FALSE) +
  annotate("text", label = bquote("Pearson'" * r * " = " * .(round(pcc, 3))), x="1/32", y=10000, size=6, hjust = 0) +
  annotate("text", label = bquote("Spearman'" * rho * " = " * .(round(scc, 3))), x="1/32", y=3000, size=6, hjust = 0)
dev.off()
# 首先将预测的大于0的保留下来;
# 这个也需要做特殊处理;
Isosceles_predict <- Isosceles_predict %>% filter(Predict > 0)
for (i in c(1:dim(Isosceles_predict)[1])) {
  if (is.na(Isosceles_predict[i,3])) {
    Isosceles_predict[i,3] <- Isosceles_predict[i,1]
  }
}
TP_FP <- length(Isosceles_predict$Transcript)
TP_FN <- length(Answer$Transcript)
TP <- length(intersect(Isosceles_predict$Transcript, Answer$Transcript))
precision <- TP/TP_FP
precision
sensitivity <- TP/TP_FN
sensitivity


# 3.5 StringTie
# ONT cDNA;
StringTie_transcript_counts_path <- "D:/AAAAABroCOLICode/SIRV_Sampled_Annotation_Fig2/Quantification/StringTie/counts_transcript.txt"
StringTie_predict <- read.table(file=StringTie_transcript_counts_path, sep='\t', header = FALSE)
StringTie_predict <- StringTie_predict[grepl("transcript", StringTie_predict$V3),]
StringTie_predict <- tidyr::separate(StringTie_predict,V9,into = c("gene_id","gene_ID","transcript_id","Transcript_id", "cov", "Transcript", "FPKM", "FPKM_value", "TPM", "Predict"),sep = ' ')
StringTie_predict_novel <- StringTie_predict[grepl("cov", StringTie_predict$cov),]
StringTie_predict <- StringTie_predict[!grepl("cov", StringTie_predict$cov),]

StringTie_predict <- tidyr::separate(StringTie_predict,Predict,into = c("Predict","a"),sep = ';')
StringTie_predict <- tidyr::separate(StringTie_predict, Transcript, into = c("Transcript","q"),sep = ';')
StringTie_predict <- StringTie_predict[,c("Transcript", "Predict")]
StringTie_predict$Predict <- as.numeric(StringTie_predict$Predict)
rm(StringTie_transcript_counts_path)
data_plot <- merge(Answer, StringTie_predict, by=c("Transcript"), all.x = TRUE)
data_plot$Predict[is.na(data_plot$Predict)] <- 0
data_plot$Predict <- data_plot$Predict + 0.1

pcc <- cor(data_plot$Truth, data_plot$Predict, method = "pearson")
pcc
scc <- cor(data_plot$Truth, data_plot$Predict, method = "spearman")
scc
data_plot$Truth <- as.factor(data_plot$Truth)

pdf("C:/Users/Lenovo/Desktop/StringTie2_Quant.pdf", width=5, height=5)
ggplot(data_plot,aes(x=Truth_xlabel, y=Predict)) + 
  geom_boxplot(size=0.5, outlier.fill="white", outlier.shape=NA, width = 0.4) + 
  geom_jitter(aes(fill=Truth, color=Truth, size=len_size), width = 0, shape = 21, alpha=0.6)+
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_blank(),
        text = element_text(size = 24),
        plot.title = element_text(size=24,hjust=0.5),
        legend.position = "bottom", 
        legend.text = element_text(size = 24))+
  labs(x="SIRV concentration (fmol/ul)", y="Isoform abundance (read count)", title="StringTie2") +
  scale_y_log10(
    breaks = c(10^(-1), 10^0, 10^1, 10^2, 10^3, 10^4),
    labels = c("0", expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4)),
    limits = c(0.1,10^4)) +
  scale_size_manual(values = c(2, 3, 4, 5, 6), 
                    name = "Isoform length (nt)", 
                    labels = c("500", "1000", "1500", "2000", "2500") ) +
  guides(fill=FALSE, color=FALSE) +
  annotate("text", label = bquote("Pearson'" * r * " = " * .(round(pcc, 3))), x="1/32", y=10000, size=6, hjust = 0) +
  annotate("text", label = bquote("Spearman'" * rho * " = " * .(round(scc, 3))), x="1/32", y=3000, size=6, hjust = 0)
dev.off()
# 首先将预测的大于0的保留下来;
StringTie_predict <- StringTie_predict %>% filter(Predict > 0)
TP_FP <- length(StringTie_predict$Transcript) + dim(StringTie_predict_novel)[1]
TP_FN <- length(Answer$Transcript)
TP <- length(intersect(StringTie_predict$Transcript, Answer$Transcript))
precision <- TP/TP_FP
precision
sensitivity <- TP/TP_FN
sensitivity







library(ggplot2)
library(patchwork)
custom_colors <- c("#f4d3d7", "#e8acb3", "#ef9c84", "#dc566f", "#851321")
custom_colors <- c("#d9e6ef", "#c5c1d9", "#8bd0dd", "#66a6c3", "#415883")
custom_colors <- c("#a7e8e2", "#84bec8", "#659cb2", "#478074", "#185f59")
custom_colors <- c("#f5e0ba", "#f5cd8d", "#e19863", "#955a4c", "#783124")


# BroCOLI
X <- c(rep("IGROV-1",3), rep("SK-OV-3",3), rep("COV504",3))
Y <- c(rep(c("IGROV-1", "COV504", "SK-OV-3"),3))
# Correlation <- c(0.62, 0.24, 0.3, 0.28, 0.28, 0.63, 0.25, 0.64, 0.3) 
Correlation <- c(0.62, 0.24, 0.3, 0.28, 0.28, 0.63, 0.25, 0.64, 0.3) 
df <- data.frame(X=X, Y=Y, Correlation=Correlation)
df$Correlation <- round(df$Correlation,2)
df$X <- factor(df$X, levels = c("IGROV-1","SK-OV-3","COV504"))
df$Y <- factor(df$Y, levels = c("IGROV-1", "SK-OV-3", "COV504"))


p1 <- ggplot(df, aes(x = X, y = Y, fill = Correlation)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(
                       colours = custom_colors, 
                       limits=c(0.2,0.7)
                       ) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 1.5,
                                barheight = 8)) +
  geom_text(aes(label=Correlation), color='black', size = 7) + 
  xlab("Single cell cluster") +
  ylab("Bulk cell line") + 
  labs(title = "BroCOLI novel transcript") +
  theme(
    axis.text.x.bottom = element_text(size=16, angle = 0, hjust = 0.5, vjust= 1),
    axis.text.y.left = element_text(size = 16), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    axis.title.x =element_text(size=24), # x轴title字体大小
    axis.title.y=element_text(size=24), # y轴title字体大小
    title = element_text(size=24, hjust=0.5),   # title字体大小
    panel.background = element_blank(),
    panel.border = element_blank(), # 去掉外边框线
    legend.position = 'none'  # 去除legend
  )
# pdf("C:/Users/Lenovo/Desktop/Fig76.pdf", width=6, height=5)
# p1
# dev.off()





X <- c(rep("IGROV-1",3), rep("SK-OV-3",3), rep("COV504",3))
Y <- c(rep(c("IGROV-1", "COV504", "SK-OV-3"),3))
Correlation <- c(0.53, 0.25, 0.27, 0.3, 0.29, 0.52, 0.28, 0.52, 0.33) # Isoquant

df <- data.frame(X=X, Y=Y, Correlation=Correlation)
df$Correlation <- round(df$Correlation,2)
df$X <- factor(df$X, levels = c("IGROV-1","SK-OV-3","COV504"))
df$Y <- factor(df$Y, levels = c("IGROV-1", "SK-OV-3", "COV504"))

p2 <- ggplot(df, aes(x = X, y = Y, fill = Correlation)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(
    colours = custom_colors, 
    limits=c(0.2,0.7)
  ) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 1.5,
                                barheight = 8)) +
  geom_text(aes(label=Correlation), color='black', size = 7) + 
  xlab("Single cell cluster") +
  ylab("Bulk cell line") + 
  labs(title = "Isoquant novel transcript") +
  theme(
    axis.text.x.bottom = element_text(size=16,angle = 0,hjust = 0.5, vjust=1),
    axis.text.y.left = element_text(size = 16), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    axis.title.x =element_text(size=24), # x轴title字体大小
    axis.title.y=element_text(size=24), # y轴title字体大小
    title = element_text(size=24, hjust=0.5),   # title字体大小
    panel.background = element_blank(),
    panel.border = element_blank(), # 去掉外边框线
    # legend.position = 'none'  # 去除legend
  )

# pdf("C:/Users/Lenovo/Desktop/Fig77.pdf", width=6, height=5)
# p2
# dev.off()

# Isosceles
X <- c(rep("IGROV-1",3), rep("SK-OV-3",3), rep("COV504",3))
Y <- c(rep(c("IGROV-1", "COV504", "SK-OV-3"),3))
Correlation <- c(0.43, 0.23, 0.23, 0.26, 0.22, 0.46, 0.26, 0.48, 0.22) # Isosceles

df <- data.frame(X=X, Y=Y, Correlation=Correlation)
df$Correlation <- round(df$Correlation,2)
df$X <- factor(df$X, levels = c("IGROV-1","SK-OV-3","COV504"))
df$Y <- factor(df$Y, levels = c("IGROV-1", "SK-OV-3", "COV504"))

p3 <- ggplot(df, aes(x = X, y = Y, fill = Correlation)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(
    colours = custom_colors, 
    limits=c(0.2,0.7)
  ) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 1.5,
                                barheight = 8)) +
  geom_text(aes(label=Correlation), color='black', size = 7) + 
  xlab("Single cell cluster") +
  ylab("Bulk cell line") + 
  labs(title = "Isosceles novel transcript") +
  theme(
    axis.text.x.bottom = element_text(size=16,angle = 0,hjust = 0.5, vjust=1),
    axis.text.y.left = element_text(size = 16), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    axis.title.x =element_text(size=24), # x轴title字体大小
    axis.title.y=element_text(size=24), # y轴title字体大小
    title = element_text(size=24, hjust=0.5),   # title字体大小
    panel.background = element_blank(),
    panel.border = element_blank(), # 去掉外边框线
    # legend.position = 'none'  # 去除legend
  )






pdf("C:/Users/Lenovo/Desktop/Fig78.pdf", width=7, height=21)
(p1 / p2 / p3)  +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.text = element_text(size = 16))
dev.off()











