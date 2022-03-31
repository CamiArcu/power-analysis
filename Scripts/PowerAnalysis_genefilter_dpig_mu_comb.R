library(POWSC) # Requires R > =  4.1
library(tidyverse)
library(purrr)
library(viridis)
library(scran)
library(DESeq2)
library(zinbwave)
library(BiocParallel)
library(furrr)
library(RhpcBLASctl)
library(gridExtra)
library(pROC)
library(limma)
library(splatter)
plan(multisession, workers = 9)
  
source("Scripts/runPOWSC_genefilter_furrr_dpig_mu_comb.R")
setwd("/home/usuario/Documentos/Cami/PowerAnalysis/")
sce <- readRDS("Data/sce_Dentate10x_r.rds")
sce <- sce[, colSums(counts(sce))>2000]
colSums(counts(sce)) %>% summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2002    2272    2535    2618    2862    5152 
saveRDS(sce, "Outputs/sce_orig_220129.rds")
est_Paras <- Est2Phase(sce, low.prob = 0.99)


# simData = Simulate2SCE_explicada_Zex2(mu = 2, sdP = 0.5, dpi.g = 0.3, n= 250, estParas1 = est_Paras, estParas2 = est_Paras)
# de = runDE_explicada(simData$sce, DE_Method = "DESeq2zw")
# pow1_rslt = Pow_Disc_different_DEgenes(de, simData = simData)
# mu = 2
# sdP = 0.5
# dpi.g = 0.3
# n= 250
# estParas1 = est_Paras
# estParas2 = est_Paras
# perDE = 0.05


args =  commandArgs(trailingOnly=TRUE)
contador = args[1]
print(args)
class(as.numeric(args[3]))
print(as.numeric(args[3]))
sim_size = eval(parse(text=args[2])) %>% as.list()
sim_mu = as.numeric(args[3])
sim_sd = as.numeric(args[4])
sim_dpi.g = as.numeric(args[5])
start <- Sys.time()
pow_rslt_t_DESeq <- future_map(sim_size, runPOWSC_genefilter_furrr, sim_mu, sim_sd, sim_dpi.g,
                               est_Paras = est_Paras,
                               per_DE=0.05, # Porcentaje de genes diferencialmente expresados
                               DE_Method = "DESeq2zw",
                               Cell_Type = "PW")
names(pow_rslt_t_DESeq) <- c(as.character(unlist(sim_size)))
end <- Sys.time()
print(paste0(round(end - start, 0), " ", units(end - start)))
size = paste(unlist(sim_size), collapse = "_")
saveRDS(pow_rslt_t_DESeq, paste0("Outputs/pow_scDESeq_220202_size_", size, "_mu_", sim_mu , "_dpig_", sim_dpi.g, "_round_", contador, ".rds"))
rm(pow_rslt_t_DESeq)

# #With POWSC
# sim_size = 500
# sim_mu = 0
# sim_sd = 0
# sim_dpi.g = 0
# simData = Simulate2SCE_explicada_Zex2(mu =sim_mu, sdP = sim_sd, dpi.g = sim_dpi.g, n=sim_size, estParas1 = est_Paras, estParas2 = est_Paras)
# saveRDS(simData, "Outputs/sim_POWSC_220129.rds")
# 
# # With splatter
# params <- splatEstimate(sce)
# sim <- splatSimulate(params, batchCells = 500)
# saveRDS(sim, "Outputs/sim_splat_220129.rds")
# 
# # With zinb
# params.zinb <- zinbEstimate(sce)
# sim.zinb <- zinbSimulate(params.zinb, model = zinbwave::zinbModel(n=500, J= 14545))
# saveRDS(sim.zinb, "Outputs/sce_zinb_220129.rds")
# 
# zinb <- readRDS("Outputs/sce_zinb_220129.rds")
# sce <- readRDS("Outputs/sce_orig_220129.rds")
# sce<- sce[, sample(500)]
# POWSC <- readRDS("Outputs/sim_POWSC_220129.rds")$sce
# sim <- readRDS("Outputs/sim_splat_220129.rds")
# sces <- list("orig" = sce, "splat" = sim, "POWSC" = POWSC, "zinb-wave" = zinb)
# compare <- compareSCEs(sces)
# p <- makeCompPanel(compare)
# ggsave("Outputs/Compare_simulators_220129.png", p)
# differences <- diffSCEs(sces, "orig")
# p <- makeDiffPanel(differences)
# ggsave("Outputs/Difference_simulators_220129.png", p)
# 
# ObtainKSd <- function(sceList, ref){
#   sceRef <- sceList[[ref]]
#   sceList <- sceList[-which(names(sceList) == ref)]
#   ks = do.call(rbind, map2(sceList, 1:length(sceList), function(x, y){
#     mu = ks.test(rowMeans(counts(x)), rowMeans(counts(sceRef)))$statistic
#     var = ks.test(rowVars(counts(x)), rowVars(counts(sceRef)))$statistic
#     pig = ks.test(rowMeans(counts(x)>0), rowMeans(counts(sceRef)>0))$statistic
#     l = ks.test(colSums(counts(x)), colSums(counts(sceRef)))$statistic
#     p0 = ks.test(colMeans(counts(x)>0), colMeans(counts(sceRef)>0))$statistic
#     data.frame(mu, var, pig, l, p0, sim.type = names(sceList)[y])
#   }))
#   KSdata = pivot_longer(ks, c(mu, var, pig, l, p0), names_to = "variable", values_to = "KSd")
#   ggplot(KSdata, aes(x = as.factor(variable), y = KSd, fill = sim.type)) + 
#     geom_col(position = position_dodge())
# }
# 
# p <- ObtainKSd(sces, "orig")
# ggsave("Outputs/KSdistance_simulators_220129.png", p)
# 
# #Two groups median differences equal 8. Splat
# sim.groups <- splatSimulate(params, batchCells = 500, group.prob = c(0.5, 0.5), method = "groups",
#                             verbose = FALSE, de.prob = 0.05, de.downProb = 0.5, de.facLoc = 5)
# 
# saveRDS(sim.groups, "Outputs/sim.groups_splat_220129.rds")
#         
# sceb <- readRDS("Outputs/sim_splat_220129.rds")
# sceb <- sim.groups
# blas_set_num_threads(5)
# nms <- c("counts", setdiff(assayNames(sceb), "counts"))
# assays(sceb) <- assays(sceb)[nms]
# #assay(sceb) <- as.matrix(assay(sceb))
# # Estimates DE genes
# # Seteo a la variable de los 'grupos' o 'cellTypes' para que esté en formato factor
# #colData(sceb)$cellTypes <- factor(x = colData(sceb)$cellTypes)
# colData(sceb)$Group <- factor(x = colData(sceb)$Group)
# # Seteo los rownames de colData y los colnames de assay. Si todo está bien, deberían ser los mismos
# rownames(colData(sceb)) <- 1:length(colData(sceb)[,1])
# colnames(assay(sceb)) <- 1:length(assay(sceb)[1,])
# sceb <-sceb[rowSums(assay(sceb))>0,]
# sceb <- zinbwave::zinbwave(sceb, K=0, observationalWeights=TRUE, epsilon=1e12, BPPARAM = BiocParallel::MulticoreParam(5))
# print("ya pase por zinwave")
# sceb <- scran::convertTo(sceb, type = "DESeq2")
# sceb <- scran::computeSumFactors(sceb)
# #design(sceb) <- formula(~ cellTypes)
# design(sceb) <- formula(~ Group)
# sd_vst <- DESeq2::varianceStabilizingTransformation(sceb, fitType = "local", blind = TRUE) %>%
#   assay() %>% apply(MARGIN = 1, sd) %>% sort(decreasing = T)
# #VariableFeatures <- names(which(sd_vst>1))
# VariableFeatures <- names(sort(sd_vst, decreasing = T)[1:2000])
# sceb <- DESeq2::estimateDispersions(object = sceb, fitType = "local")
# sceb <- sceb[VariableFeatures, ]
# sceb <- DESeq2::nbinomLRT(object = sceb, reduced= ~1, minmu = 1e-6)
# sceb <-  DESeq2::replaceOutliers(sceb, minReplicates = Inf)
# #res_DE <- DESeq2::results(object = sceb, contrast = c("cellTypes", "celltype1", "celltype2"), alpha = 0.05, independentFiltering=FALSE)
# res_DE <- DESeq2::results(object = sceb, contrast = c("Group", "Group1", "Group2"), alpha = 0.05, independentFiltering=FALSE)
# print("ya pase por results de DESeq2")
# res_DE <- data.frame(res_DE) %>% mutate(geneIndex = str_split_fixed(rownames(res_DE), "g", 2)[,2]) %>%
#   mutate(geneIndex = as.integer(geneIndex))
# res_DE_p_vals <- dplyr::select(res_DE, geneIndex, pval = pvalue, fdr = padj)
# res_DE_p_vals[is.na(res_DE_p_vals)] <- 1
# 
# 
# saveRDS(sceb, "Outputs/sceb_mu8_simZinb_220129.rds")
# data_sce <- do.call(cbind, select2(rowData(sceb)@listData, Group_Group2_vs_Group1, LRTPvalue, baseMean)) %>% as.data.frame()
# data_sce <- mutate(data_sce, padj = p.adjust(LRTPvalue, "BH")) %>%
#   mutate(sig = if_else(padj<0.05, "sig", "nosig"))
# names(data_sce)[1] <- "lfc"
# p1 <- ggplot(data_sce, aes(x= log2(baseMean), y = lfc, color =  as.factor(sig))) + geom_point()
# p1
# ggsave("Outputs/plotMA_mu8_simZinb_220129.png", p1)
# 
# select2 <- function(.x, ...) {
#   vars <- rlang::names2(.x)
#   vars <- tidyselect::vars_select(vars, ...)
#   .x[vars]
# }
# 
# #     Table to store the log2FC
# table <- dplyr::select(res_DE, geneIndex, cont = log2FoldChange)
# blas_set_num_threads(10)
# 
# 
# alpha = 0.05
# strata = seq(0, 1, by = 0.2)
# DErslt <- list(table = table, cont = res_DE_p_vals)
# 
# Pow_only_DE1 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE1_only)
# Pow_only_DE2 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE2_only)
# Pow_DE1_DE2 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE1_DE2)
# Pow_DE1_or_DE2 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE1_or_DE2)
# data_sce <- simData$data_sce[simData$ix.DE1_DE2, ]
# ix.DE1_DE2_co <- filter(data_sce, (lfc > 0 & pi.df > 0) | (lfc < 0 & pi.df < 0))$id
# ix.DE1_DE2_op <- filter(data_sce, (lfc > 0 & pi.df < 0) | (lfc < 0 & pi.df > 0))$id
# Pow_DE1_DE2_co <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = ix.DE1_DE2_co)
# Pow_DE1_DE2_op <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = ix.DE1_DE2_op)
# 
# DEid <- ix.DE1_DE2_op
# ngenes <- nrow(simData$sce)
# # Genero dos vectores con 0, de longuitud = numero de genes
# Zg = Zg2 = rep(0, ngenes)
# # A los genes DA (de verdad), les pongo un 1
# Zg[DEid] <- 1
# # Busco el ID de los genes que poseen un pi.df (proporcion de cambio a genes activos) mayor a un umbral (delta)
# # ix <- which(abs(pi.df) > delta) # Comenté esta linea porque con la asignacio de un dpi.g definado, ya no tiene sentido
# # A los genes DA (de verdad), pero que superan el umbral de pi.df les pongo un 1
# Zg2[DEid] <- 1 #Zg2[DEid[ix]] <- 1 #Antes de la asignacion de un dpi.g esta linea era lo comentado
# # Extraigo la matrix de datos simulados
# sce <- simData$sce
# # Extraigo el número de células totales
# ntotal <- ncol(sce)
# # Extraigo la matriz original (de counts) de datos
# Y <- round(assays(sce)[[1]])
# # Para cada gen, calculo la proporción de células en las cuales está inactivo
# rate0 <- rowMeans(Y == 0)
# # De los genes tanto DE1 como DE2, identifico cuales no pertenecen al subset de DEs que quiero analizar ahora
# # Si analizo todos los genes DEs ahora, ix.not = empty
# ix.not <- simData$ix.DE1_or_DE2[!(simData$ix.DE1_or_DE2 %in% DEid)]
# # De los genes analizados, saco los genes que son DEs pero no pertenecen al grupo de DEs analizados ahora
# ix.yes <- DErslt$cont$geneIndex[!(DErslt$cont$geneIndex %in% ix.not)]
# # De los genes totales, me quedo con los finalmente analizados
# rate0 = rate0[ix.yes]
# # Identifico los genes que están inactivos entre 1% y 99% de las celulas (excluyo los puntos extremos, por ej: genes inactivos en todas las célulaso (100%) en ninguna (0%))
# #ix.keep <- intersect(which(rate0 < 0.99), which(rate0 >  0.01)) #esto me va a dar el índice respecto de rate0
# # Divido la proporcion de inactivos de los genes anteriores, en estratos
# xgr <- cut(rate0, strata)#xgr <- cut(rate0[ix.keep], strata)
# #ix.keep <- DErslt$cont$geneIndex[ix.keep] # esto me va a dar el índice respecto de la matriz original
# # De los genes totales, me quedo con los finalmente analizados
# # Tanto para el vector que tiene marcados los genes DA
# Zg <- Zg[ix.yes] 
# # Como el que tiene los gene DA que superan el umbral
# Zg2 = Zg2[ix.yes]
# # Extraigo el FDR del analisis de DE de los datos simulados
# # Esto lo hago con filter porque DErstl es un data frame con un nrows mucho mas chico que Zg2 o rate0, osea, es mas chico el nrows de sce
# # entonces si hiciese DErslt$cont$fdr[ix.yes], estaria agarrando los indices de una tabla de 1000 filas en vez de 10000 filas
# fdrvec <- dplyr::filter(DErslt$cont, geneIndex %in% ix.yes)$fdr
# pvalvec <- dplyr::filter(DErslt$cont, geneIndex %in% ix.yes)$pval
# 
# which(names(sort(sd_vst, decreasing = T)) %in% paste0("g", ix.DE1_DE2_op))
# # [1]  4425  5969  6924  7449  8500
# # [6]  8999 10049 11930 13028 13046
# # [11] 13176 13745

#all_runs <- readRDS("Outputs/pow_scDESeq_220126_size_2000_mu_8_dpig_0.1_round_1.rds")
