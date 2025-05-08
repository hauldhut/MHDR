#### read stored plot and plot then save to pdf file (For MH and H Methods)
library(ggplot2)
library('cowplot')

#Figure 3
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_M1-M2-M3-M12-M13-M23-M123_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_M1-M2-M3-M12-M13-M23-M123_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure 4

p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_H1-H2-H3-MH12-MH13-MH23-MH123_PREDICT_gamma_0.5_delta_0.5.rdata")
p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_H1-H2-H3-MH12-MH13-MH23-MH123_CHEM_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure4_AUROC.pdf", width = 14, height = 7)

p3 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_H1-H2-H3-MH12-MH13-MH23-MH123_PREDICT_gamma_0.5_delta_0.5.rdata")
p4 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_H1-H2-H3-MH12-MH13-MH23-MH123_CHEM_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, p3, p4, labels=c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2)

ggsave("~/Manuscripts/99MHDR/Figures/Figure4_AUROCnAUPRC.pdf", width = 14, height = 14)


#Figure Mono_k_simThes
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_M1-M1_10-M1_15-M1_03_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMono_k_simThes_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_M1-M1_10-M1_15-M1_03_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMono_k_simThes_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure Multi_k_simThes
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_M123-M1_1023-M1_1523-M1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMulti_k_simThes_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_M123-M1_1023-M1_1523-M1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMulti_k_simThes_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure H_k_simThes --> S3
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_H1-H1_10-H1_15-H1_03_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureH_k_simThes_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_H1-H1_10-H1_15-H1_03_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureH_k_simThes_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure H_k_simThes (CHEM)
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_H1-H1_10-H1_15-H1_03_CHEM_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureH_CHEM_k_simThes_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_H1-H1_10-H1_15-H1_03_CHEM_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureH_CHEM_k_simThes_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure MH_k_simThes --> S4
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_MH123-MH1_1023-MH1_1523-MH1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMH_k_simThes_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_MH123-MH1_1023-MH1_1523-MH1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMH_k_simThes_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure MH_k_simThes (CHEM)
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_MH123-MH1_1023-MH1_1523-MH1_0323_CHEM_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMH_CHEM_k_simThes_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_MH123-MH1_1023-MH1_1523-MH1_0323_CHEM_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/FigureMH_CHEM_k_simThes_AUROCnAUPRC.pdf", width = 14, height = 7)


#Figure 3 with k=10
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_M1_10-M1_102-M1_103-M1_1023_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_k10_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_M1_10-M1_102-M1_103-M1_1023_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_k10_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure 3 with k=15
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_M1_15-M1_152-M1_153-M1_1523_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_k15_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_M1_15-M1_152-M1_153-M1_1523_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_k15_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure 3 with sim>=0.3
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_M1_03-M1_032-M1_033-M1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_sim03_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_M1_03-M1_032-M1_033-M1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure3_sim03_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure 4 with k=10
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_H1_10-MH1_102-MH1_103-MH1_1023_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure4_k10_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_H1_10-MH1_102-MH1_103-MH1_1023_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure4_k10_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure 4 with k=15
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_H1_15-MH1_152-MH1_153-MH1_1523_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure4_k15_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_H1_15-MH1_152-MH1_153-MH1_1523_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure4_k15_AUROCnAUPRC.pdf", width = 14, height = 7)

#Figure 4 with sim>=0.3
p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUROC_H1_03-MH1_032-MH1_033-MH1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

print(p1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure4_sim03_AUROC.pdf", width = 7, height = 7)

p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_AUPRC_H1_03-MH1_032-MH1_033-MH1_0323_PREDICT_gamma_0.5_delta_0.5.rdata")

plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
ggsave("~/Manuscripts/99MHDR/Figures/Figure4_sim03_AUROCnAUPRC.pdf", width = 14, height = 7)
