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

setwd("/home/usuario/Documentos/Cami/PowerAnalysis/")  
source("Scripts/runPOWSC_genefilter_furrr_dpig_mu_comb.R")

#### Different ways of defining data that with a specific average of UMIs/Cell (library size)

# # First selecting a subset of the data -> Problems: maybe other parameters of the data varies, not only the library size
# sce <- readRDS("Data/sce_Dentate10x_r.rds")
# sce <- sce[, colSums(counts(sce))>2000]
# colSums(counts(sce)) %>% summary()
# #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #2002    2272    2535    2618    2862    5152 
# est_Paras <- Est2Phase(sce, low.prob = 0.99)
# #Simulate data and check if the distribution of library size is simmilar to the original data
# simData = Simulate2SCE_explicada_Zex2(mu = 2, sdP = 0.5, dpi.g = 0.3, n= 1000, estParas1 = est_Paras, estParas2 = est_Paras)
# colSums(counts(simData$sce)) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 1869    2656    3220    3404    3974   22473 
# # Again but with aother minimun threshold
# sce <- readRDS("Data/sce_Dentate10x_r.rds")
# sce <- sce[, colSums(counts(sce))>1000]
# colSums(counts(sce)) %>% summary()
# #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #1091    1931    2320    2338    2713    5152
# # All the data
# sce <- readRDS("Data/sce_Dentate10x_r.rds")
# colSums(counts(sce)) %>% summary()
# #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #1091    1931    2320    2338    2713    5152
# summary(rowMeans(counts(sce) == 0)) # Median = 0.96
# #Entionces no voy  ausar el data set del giro dentado porque quiero probar con datasets que tengan 1000, 3000 y 5000 reads por celula en promedio
# # y en el de dentado, 3/4 partes del data set tienen hasta menos de 3000 UMI/Cell
# # Voy a usat un data set que tenga o llegue a un mayor numero de UMI/Cell
# # With de SNr data
# sce_SNr <- readRDS("Data/sce_SNr.rds")
# sce_SNr <- sce_SNr[, colSums(counts(sce_SNr)) > 0]
# colSums(counts(sce_SNr)) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 474     894    1338    1812    2151   42878 
# dim(sce_SNr) #26434 10049
# rowMeans(counts(sce_SNr) == 0) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.0000  0.9597  0.9967  0.9607  0.9998  1.0000
# sce_SNr <- as.matrix(counts(sce_SNr))
# est_Paras <- Est2Phase(sce_SNr, low.prob = 0.99)
# #Error in asMethod(object) :Cholmod error 'out of memory' at file ../Core/cholmod_memory.c, line 146
# # Again but with lower cell lines
# sce_SNr <- readRDS("Data/sce_SNr.rds")
# sce_SNr <- sce_SNr[, colSums(counts(sce_SNr)) > 0]
# sce_SNr <- sce_SNr[, sample(1:10049, 3000)]
# colSums(counts(sce_SNr)) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 478     898    1338    1822    2131   38628 
# rowMeans(counts(sce_SNr) == 0) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.0000  0.9597  0.9967  0.9607  1.0000  1.0000
# sce_SNr <- as.matrix(counts(sce_SNr))
# est_Paras <- Est2Phase(sce_SNr, low.prob = 0.99)
# simData = Simulate2SCE_explicada_Zex2(mu = 2, sdP = 0.5, dpi.g = 0.3, n= 1000, estParas1 = est_Paras, estParas2 = est_Paras)
# colSums(counts(simData$sce)) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 571    1381    1818    2189    2574   22442 
# # Again but with aother minimun threshold
# sce_SNr <- readRDS("Data/sce_SNr.rds")
# sce_SNr <- sce_SNr[, colSums(counts(sce_SNr)) > 1000]
# colSums(counts(sce_SNr)) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2001    2375    2912    3552    3912   42878
# sce_SNr <- as.matrix(counts(sce_SNr))
# est_Paras <- Est2Phase(sce_SNr, low.prob = 0.99)
# simData = Simulate2SCE_explicada_Zex2(mu = 2, sdP = 0.5, dpi.g = 0.3, n= 1000, estParas1 = est_Paras, estParas2 = est_Paras)
# colSums(counts(simData$sce)) %>% summary()
# # Again but with aother minimun threshold
# sce_SNr <- readRDS("Data/sce_SNr.rds")
# sce_SNr <- sce_SNr[, colSums(counts(sce_SNr)) > 3000]
# colSums(counts(sce_SNr)) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 3003    3424    4028    4853    5228   42878 
# sce_SNr <- as.matrix(counts(sce_SNr))
# est_Paras <- Est2Phase(sce_SNr, low.prob = 0.99)
# simData = Simulate2SCE_explicada_Zex2(mu = 2, sdP = 0.5, dpi.g = 0.3, n= 1000, estParas1 = est_Paras, estParas2 = est_Paras)
# colSums(counts(simData$sce)) %>% summary()
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 2778    4151    4854    5437    5959   33252
# 
# ### Trying to select cell until reach the mean of UMI/Cell desired
# sample_by_mean <- function(thr, matrich){
#   data <- matrich[, sample(1:dim(matrich)[2], 100)]
#   prom <- colSums(matrich) %>% mean()
#   while (prom < thr){
#     others <- colnames(matrich)[!(colnames(matrich) %in% colnames(data))]
#     if (others>100){
#       data <- cbind(data, matrich[, sample(others, 100)])
#       prom <- colSums(matrich) %>% mean()
#     } else {
#       prom <- thr
#     }
#   }
#   return(data)
# }
# sce_1000 <- sample_by_mean(1000, sce_SNr)
# sce_3000 <- sample_by_mean(3000, sce_SNr) #Error
# sce_5000 <- sample_by_mean(5000, sce_SNr) #Error
# 
# # Modifying just only the library size -> all other arameters remains the same
# sce <- readRDS("Data/sce_Dentate10x_r.rds")
# sce <- sce[, colSums(counts(sce))>2000]
# est_Paras <- Est2Phase(sce, low.prob = 0.99)
# # Library size equivalences
# # 5 = 1600, 1000 = 3000, 2000 = 5000
# new_sf <- 2000/1e6 #2000 here equals 5000 uin library size, don know why
# est_Paras$sf <- new_sf
# simData = Simulate2SCE_explicada_Zex2(mu = 2, sdP = 0.5, dpi.g = 0.3, n= 1000, estParas1 = est_Paras, estParas2 = est_Paras)
# colSums(counts(simData$sce)) %>% summary()
# # This work pretty well

runPOWSC_genefilter_furrr <- function (sim_size, mu, sdP, dpi.g, per_DE = 0.05, 
                                       est_Paras, DE_Method = c("MAST", "SC2P", "DESeq2", "ScDESeq2",  "SeuratDESeq2", "DESeq2zw"), Cell_Type = c("PW",                                                                                                     "Multi"), multi_Prob = NULL, alpha = 0.1, disc_delta = 0.1, 
                                       cont_delta = 0.5, lbsz) 
{
  DE_Method = match.arg(DE_Method)
  Cell_Type = match.arg(Cell_Type)
  pow_rslt = NULL
  if (Cell_Type == "PW") {
    print(paste0("Running simulation with ", sim_size, " cells and mu = ", mu, ", sd = ", sdP, " and dpi.g = ", dpi.g, " and lbsz = ", lbsz))
    est_Paras$sf <- lbsz/1e6
    simData = Simulate2SCE_explicada_Zex2(mu = mu, sdP = sdP, dpi.g = dpi.g, n=sim_size, estParas1 = est_Paras, estParas2 = est_Paras)
    lbsz_dist <- colSums(counts(simData$sce))
    rate0cell_dist <- colMeans(counts(simData$sce) == 0)
    de = runDE_explicada(simData$sce, DE_Method = DE_Method)
    pow1_rslt = Pow_Disc_different_DEgenes(de, simData = simData)
    pow2_rslt = Pow_Cont_different_DEgenes(de, simData = simData)
    # if(length(simData$ix.DE2)!= 0){
    #   pow2_rslt = Pow_Cont_explicada_genefilter(de, simData = simData)
    # } else {
    #   pow2_rslt = NULL
    # }
    tmp_rslt = list(pow1 = pow1_rslt, pow2 = pow2_rslt, 
                    ngenes = simData$ngene, data_sce = simData$data_sce, ix.DE1_DE2 = simData$ix.DE1_DE2, 
                    ix.DE1_only = simData$ix.DE1_only, ix.DE2_only = simData$ix.DE2_only, ix.de = de$cont$geneIndex,
                    lbsz_dist = lbsz_dist, rate0cell_dist = rate0cell_dist) # pow2 = pow2_rslt
    pow_rslt[[toString(dpi.g)]][[toString(mu)]] = tmp_rslt
    print("llegue hasta aca")
  }
  class(pow_rslt) = "POWSC"
  return(pow_rslt)
}


sce <- readRDS("Data/sce_Dentate10x_r.rds")
sce <- sce[, colSums(counts(sce))>2000]
est_Paras <- Est2Phase(sce, low.prob = 0.99)

args =  commandArgs(trailingOnly=TRUE)
contador = args[1]
print(args)
class(as.numeric(args[3]))
print(as.numeric(args[3]))
sim_size = eval(parse(text=args[2])) %>% as.list()
sim_mu = as.numeric(args[3])
sim_sd = as.numeric(args[4])
sim_dpi.g = as.numeric(args[5])
lib_size = as.numeric(args[6])

start <- Sys.time()
pow_rslt_t_DESeq <- future_map(sim_size, runPOWSC_genefilter_furrr, sim_mu, sim_sd, sim_dpi.g,
                               est_Paras = est_Paras,
                               per_DE=0.05, # Porcentaje de genes diferencialmente expresados
                               DE_Method = "DESeq2zw",
                               Cell_Type = "PW", lbsz = lib_size)
names(pow_rslt_t_DESeq) <- c(as.character(unlist(sim_size)))
end <- Sys.time()
print(paste0(round(end - start, 0), " ", units(end - start)))
size = paste(unlist(sim_size), collapse = "_")
saveRDS(pow_rslt_t_DESeq, paste0("Outputs/pow_scDESeq_220225_size_", size, "_mu_", sim_mu , "_dpig_", sim_dpi.g, "_round_", contador, "_lbsz_", lib_size,".rds"))
rm(pow_rslt_t_DESeq)

