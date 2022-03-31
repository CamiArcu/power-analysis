#### Analizing the results ####
library(POWSC) # Requires R > =  4.1
library(tidyverse)
library(purrr)
library(viridis)
library(gridExtra)
library(pROC)
setwd("/home/usuario/Documentos/Cami/PowerAnalysis/")

# Here I charge the file names of the respective run
all_runs <- list.files("Outputs") %>%
  grep(pattern = "pow_scDESeq_220202", value = T)

# Obtain details of the run as: number of Cells, abs.lfc (mu), dpig, Round, etc.
Details <- str_split_fixed(all_runs, pattern = "_", 12)[, c(5:12)]
for (i in 1:nrow(Details)){
  if(Details[i, 1] != "250"){
    Details[i, 4:8] <- Details[i, 3:7]
  } else{
    Details[i, 1] <- "250_500"
  }
}
Details[,8] <-  str_split_fixed(Details[,8], pattern = fixed("."), 2)[,1]
Details[, c(1,4,6,8)] %>% `colnames<-`(c("Cells", "abs.lfc", "dpi.g", "Round")) %>% data.frame() -> Details

# Here I retrive Data per strata.
# Could be: power, alpha, CD, TD, FD of the strata of the rate0(pow1) and mu(pow2)
getData<- function(POWSCobj, name, type, var, subs){
  dpi.g <- names(POWSCobj[[1]]) %>% as.numeric()
  abs.lfc <- names(POWSCobj[[1]][[1]]) %>% as.numeric()
  analize <- function(){
    if(type == "I"){
      POWmatfdr <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow1[[subs]][[1]][[var]])) %>%
        `colnames<-` (names(lapply(POWSCobj, function(x) x[[1]][[1]]$pow1[[subs]][[1]]$FDR)[[1]]))
      POWmatpval <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow1[[subs]][[2]][[var]])) %>%
        `colnames<-` (names(lapply(POWSCobj, function(x) x[[1]][[1]]$pow1[[subs]][[2]]$FDR)[[1]]))
    } else if(type == "II"){
      POWmatfdr <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow2[[subs]][[1]][[var]])) %>%
        `colnames<-` (names(lapply(POWSCobj, function(x) x[[1]][[1]]$pow2[[subs]][[1]]$FDR)[[1]]))
      POWmatpval <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow2[[subs]][[2]][[var]])) %>%
        `colnames<-` (names(lapply(POWSCobj, function(x) x[[1]][[1]]$pow2[[subs]][[2]]$FDR)[[1]]))
    } else{
      return(print(paste0(type, " is not a type. Indicate 'I' or 'II'")))
    }
    POWmatfdr <- as.data.frame(POWmatfdr) %>%
      pivot_longer(everything(), names_to = "Strata", values_to = var) %>%
      mutate(nCells = rep(as.integer(names(POWSCobj)), each = ncol(POWmatfdr)),
             dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name), subs = subs)
    POWmatpval <- as.data.frame(POWmatpval) %>%
      pivot_longer(everything(), names_to = "Strata", values_to = var) %>%
      mutate(nCells = rep(as.integer(names(POWSCobj)), each = ncol(POWmatpval)),
             dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name), subs = subs)
    Data <- left_join(POWmatfdr, POWmatpval, by = c("dpi.g", "abs.lfc", "Round", "subs", "nCells", "Strata"))
    names(Data)[c(2, 8)] <- c(paste0(var,".fdr"), paste0(var,".pval"))
    return(Data)
  }
  if (subs == "DE1_DE2_op"){
    rslt <- do.call(rbind, lapply(POWSCobj, function(x){
      if (class(x[[1]][[1]][["pow1"]][[subs]]) == "character"){
        msg <- x[[1]][[1]][["pow1"]][[subs]]
        Data <- data.frame(Strata = NA, var1 = msg, nCells = rep(as.integer(names(POWSCobj)), each = 1),
               dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name), subs = subs, var2 = msg)
        names(Data)[c(2, 8)] <- c(paste0(var,".fdr"), paste0(var,".pval"))
        return(Data)
      } else{
        return(analize())
      }
    }))
    return(unique(rslt))
  } else{
    return(analize())
  }
}

# Here I get the marginal variables:
# It could be: power, alpha and auc
getData_mar<- function(POWSCobj, name, type, var, subs){
  dpi.g <- names(POWSCobj[[1]]) %>% as.numeric()
  abs.lfc <- names(POWSCobj[[1]][[1]]) %>% as.numeric()
  analize <- function(){
    if(type == "I"){
      POWfdr <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow1[[subs]][[1]][[var]]))
      POWpval <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow1[[subs]][[2]][[var]]))
    } else if(type == "II"){
      POWfdr <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow2[[subs]][[1]][[var]]))
      POWpval <- do.call(rbind, lapply(POWSCobj, function(x) x[[1]][[1]]$pow2[[subs]][[2]][[var]]))
      
    } else{
      return(print(paste0(type, " is not a type. Indicate 'I' or 'II'")))
    }
    Data<-data.frame(var1 = POWfdr, var2 = POWpval,
                     nCells = rep(as.integer(names(POWSCobj)), each = 1),
                     dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name), subs = subs)
    names(Data)[1:2]<- c(paste0(var, ".fdr"), paste0(var, ".pval"))
    return(Data)
  }
  if (subs == "DE1_DE2_op"){
    rslt <- do.call(rbind, lapply(POWSCobj, function(x){
      if (class(x[[1]][[1]][["pow1"]][[subs]]) == "character"){
        msg <- x[[1]][[1]][["pow1"]][[subs]]
        Data<-data.frame(var1 = msg, var2 = msg,
                         nCells = rep(as.integer(names(POWSCobj)), each = 1),
                         dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name), subs = subs)
        names(Data)[1:2]<- c(paste0(var, ".fdr"), paste0(var, ".pval"))
        return(Data)
      } else{
        return(analize())
      }
    }))
    return(unique(rslt))
  } else{
    return(analize())
  }
}

#Color of the cells
cells_colors <- RColorBrewer::brewer.pal(5, "Blues")[2:5]

# This function plot the variables per strata
plot_function<- function(A_Tibble, tit, tmpxlab, breaks, var, strata){
  if (strata == "mu"){
    A_Tibble$Strata <- as.factor(A_Tibble$Strata)
    levels(A_Tibble$Strata) <- c(levels(A_Tibble$Strata)[-grep("Inf", levels(A_Tibble$Strata))], grep("Inf", levels(A_Tibble$Strata), value = T))
  }
   ggplot(A_Tibble, aes(x = Strata, y = get(var), size = dpi.g,
                       color = nCells)) +
    geom_point(alpha = 0.5, position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = cells_colors) +
    tit +
    scale_y_continuous(breaks = breaks, limits = c(0, 1 + 0.1),
                       labels = format(breaks, nsmall = 1)) +
    labs(y = var, x = tmpxlab, color = "# Total Cells\n(both groups summed)", size = "Delta pi.g" ) +
    my_theme
}

plot_function2<- function(A_Tibble, tit, tmpxlab, breaks, var, strata){
  if (strata == "mu"){
    A_Tibble$Strata <- as.factor(A_Tibble$Strata)
    levels(A_Tibble$Strata) <- c(levels(A_Tibble$Strata)[-grep("Inf", levels(A_Tibble$Strata))], grep("Inf", levels(A_Tibble$Strata), value = T))
  }
  ggplot(A_Tibble, aes(x = Strata, y = get(var), size = TD.fdr,
                       color = nCells)) +
    geom_point(alpha = 0.5, position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = cells_colors) +
    tit +
    scale_y_continuous(breaks = breaks, limits = c(0, 1 + 0.1),
                       labels = format(breaks, nsmall = 1)) +
    labs(y = var, x = tmpxlab, color = "# Total Cells\n(both groups summed)", size = "#Genes\nAssigned" ) +
    my_theme
}

# The theme of most plots
my_theme <-     theme_classic() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))

# Function to remove grid
remove_grid <- theme(axis.line = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank())

# Get the legend of a plot, used for strata plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Get the data.frame of the mu, pi.g of celltype1, celltype2, abs.lfc, dpig, Round, nCell
get_data_sce<- function(POWSCobj, Round){
  dpi.g <- names(POWSCobj[[1]]) %>% as.numeric()
  abs.lfc <- names(POWSCobj[[1]][[1]]) %>% as.numeric()
  Data <- do.call(rbind, lapply(POWSCobj, function(x){
    data_sce <- x[[1]][[1]]$data_sce
    data_sce <- filter(data_sce, id %in% x[[1]][[1]]$ix.de)
    mutate(data_sce, dpi.g = dpi.g, abs.lfc = abs.lfc, Round = Round)}))
  return(Data)
}

nCells_dpi.g_levels <- c("250_0.1", "250_0.2", "250_0.3", "250_0.4",
                         "500_0.1", "500_0.2", "500_0.3", "500_0.4",
                         "1000_0.1", "1000_0.2", "1000_0.3", "1000_0.4",
                         "2000_0.1", "2000_0.2", "2000_0.3", "2000_0.4")
# Plot the distribution of the aucs
p1 <- function(AUCtypeI_mar_, subs){
  p <- mutate(AUCtypeI_mar_, nCells_dpi.g = paste0(nCells, "_", dpi.g)) %>%
    ggplot(aes(x = factor(nCells_dpi.g, levels = nCells_dpi.g_levels), y = auc.fdr)) +
    geom_violin(aes(fill = as.factor(dpi.g), color = as.factor(nCells))) + geom_point() +
    labs(title = paste0("Genes - ", subs), x = "# Total Cells (both groups summed) _ Delta pi.g", y = "auc", fill = "Delta pi.g", color = "# Total Cells\n(both groups summed)") +
    my_theme +
    scale_colour_manual(values = cells_colors) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1)) +
    facet_wrap(vars(abs.lfc), ncol = 2, labeller = label_both)
  ggsave(paste0("Outputs/POWSC_auc_distribution_",  subs, "_", date, ".png"), p, width = 16,  height = 8)
}

# get the data.frame of the Zg2, xgr and fdr.
# Carefull, it is only prepared to one cell type.
get_ids<- function(POWSCobj, name, type, subs){
  dpi.g <- names(POWSCobj[[1]]) %>% as.numeric()
  abs.lfc <- names(POWSCobj[[1]][[1]]) %>% as.numeric()
  if(type == "I"){
    idsfdr <- POWSCobj[[1]][[1]][[1]]$pow1[[subs]][[1]][["ids"]]
    idspval <- POWSCobj[[1]][[1]][[1]]$pow1[[subs]][[2]][["ids"]]
  } else if(type == "II"){
    idsfdr <- POWSCobj[[1]][[1]][[1]]$pow2[[subs]][[1]][["ids"]]
    idspval <- POWSCobj[[1]][[1]][[1]]$pow2[[subs]][[2]][["ids"]]
  } else{
    return(print(paste0(type, " is not a type. Indicate 'I' or 'II'")))
  }
  idsfdr<-mutate(idsfdr, nCells = "2000",
                 dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name), subs = subs)
  idspval<-mutate(idspval, nCells = "2000",
                  dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name), subs = subs)
  Data <- left_join(idsfdr, idspval, by = c("dpi.g", "abs.lfc", "Round", "subs", "nCells", "xgr", "Zg2")) 
  names(Data)[c(3,9)]<- c("fdr", "pval")
  return(Data)
}

# Load the runs
all_runs<-lapply(paste0("Outputs/", all_runs), readRDS)
# all_runs<-lapply(all_runs, function(x) {
#   names(x[[1]][[1]][[1]])[2] <- "pow2"
#   x
# })

# Obtain the power of strata
# Type1:
#POWtypeI_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "all_DE1"))
POWtypeI_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "DE1_DE2"))
POWtypeI_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "DE1_or_DE2"))
POWtypeI_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "only_DE1"))
POWtypeI_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "only_DE2"))
POWtypeI_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "DE1_DE2_op"))
POWtypeI_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "DE1_DE2_co"))
# Obtain the alpha of strata
#ALPHAtypeI_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "alpha", subs = "all_DE1"))
ALPHAtypeI_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "alpha", subs = "DE1_DE2"))
ALPHAtypeI_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "alpha", subs = "DE1_or_DE2"))
ALPHAtypeI_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "alpha", subs = "only_DE1"))
ALPHAtypeI_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "alpha", subs = "only_DE2"))
ALPHAtypeI_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "alpha", subs = "DE1_DE2_op"))
ALPHAtypeI_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "alpha", subs = "DE1_DE2_co"))
#Type 2:
#POWtypeII_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "power", subs = "all_DE1"))
POWtypeII_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "power", subs = "DE1_DE2"))
POWtypeII_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "power", subs = "DE1_or_DE2"))
POWtypeII_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "power", subs = "only_DE1"))
POWtypeII_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "power", subs = "only_DE2"))
POWtypeII_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "power", subs = "DE1_DE2_op"))
POWtypeII_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "power", subs = "DE1_DE2_co"))
#ALPHAtypeII_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "alpha", subs = "all_DE1"))
ALPHAtypeII_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "alpha", subs = "DE1_DE2"))
ALPHAtypeII_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "alpha", subs = "DE1_or_DE2"))
ALPHAtypeII_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "alpha", subs = "only_DE1"))
ALPHAtypeII_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "alpha", subs = "only_DE2"))
ALPHAtypeII_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "alpha", subs = "DE1_DE2_op"))
ALPHAtypeII_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "alpha", subs = "DE1_DE2_co"))

# Obtain the CD (number of genes assigned to be DE, writted as TD) per strata
# Type1:
#POWtypeI_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "all_DE1"))
CDtypeI_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "TD", subs = "DE1_DE2"))
CDtypeI_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "TD", subs = "DE1_or_DE2"))
CDtypeI_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "TD", subs = "only_DE1"))
CDtypeI_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "TD", subs = "only_DE2"))
CDtypeI_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "TD", subs = "DE1_DE2_op"))
CDtypeI_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "TD", subs = "DE1_DE2_co"))
# Type2:
#POWtypeI_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "I", var = "power", subs = "all_DE1"))
CDtypeII_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "TD", subs = "DE1_DE2"))
CDtypeII_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "TD", subs = "DE1_or_DE2"))
CDtypeII_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "TD", subs = "only_DE1"))
CDtypeII_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "TD", subs = "only_DE2"))
CDtypeII_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "TD", subs = "DE1_DE2_op"))
CDtypeII_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData, type = "II", var = "TD", subs = "DE1_DE2_co"))



# Obtain the marginal power
#POWtypeI_mar_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "power.marginal", subs = "all_DE1"))
POWtypeI_mar_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "power.marginal", subs = "DE1_DE2"))
POWtypeI_mar_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "power.marginal", subs = "DE1_or_DE2"))
POWtypeI_mar_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "power.marginal", subs = "only_DE1"))
POWtypeI_mar_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "power.marginal", subs = "only_DE2"))
POWtypeI_mar_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "power.marginal", subs = "DE1_DE2_op"))
POWtypeI_mar_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "power.marginal", subs = "DE1_DE2_co"))
# Obtain the marginal alpha
#ALPHAtypeI_mar_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "alpha.marginal", subs = "all_DE1"))
ALPHAtypeI_mar_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "alpha.marginal", subs = "DE1_DE2"))
ALPHAtypeI_mar_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "alpha.marginal", subs = "DE1_or_DE2"))
ALPHAtypeI_mar_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "alpha.marginal", subs = "only_DE1"))
ALPHAtypeI_mar_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "alpha.marginal", subs = "only_DE2"))
ALPHAtypeI_mar_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "alpha.marginal", subs = "DE1_DE2_op"))
ALPHAtypeI_mar_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "alpha.marginal", subs = "DE1_DE2_co"))
# Obtain the marginal auc
#AUCtypeI_mar_all_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "auc", subs = "all_DE1"))
AUCtypeI_mar_DE1_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "auc", subs = "DE1_DE2"))
AUCtypeI_mar_DE1_or_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "auc", subs = "DE1_or_DE2"))
AUCtypeI_mar_only_DE1 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "auc", subs = "only_DE1"))
AUCtypeI_mar_only_DE2 <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "auc", subs = "only_DE2"))
AUCtypeI_mar_DE1_DE2_op <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "auc", subs = "DE1_DE2_op"))
AUCtypeI_mar_DE1_DE2_co <- do.call(rbind, map2(all_runs, Details$Round, getData_mar, type = "I", var = "auc", subs = "DE1_DE2_co"))

# Get the data.frame of the Zg2, xgr and fdr.
# Using the strata of rate0
ids_only_DE2_mu8_p1 <- all_runs[[4]][[1]][[1]][[1]]$pow1[["only_DE2"]][[1]][["ids"]]
# get the number of genes analyzed per strata
table(ids_only_DE2_mu8_p1$xgr)
# Using the strata of lfc (with runs previous of 220201_2) or mu (with runs after 220201_3, inclusive)
ids_only_DE2_mu8_p2 <- all_runs[[4]][[1]][[1]][[1]]$pow2[["only_DE2"]][[1]][["ids"]]
# get the number of genes analyzed per strata
table(ids_only_DE2_mu8_p2$xgr)

# Same but get the number of genes aasigned to be differentially expressed of only_DE2
ids_only_DE2_mu8_p1_2 <- all_runs[[1]][[1]][[1]][[1]]$pow1[["only_DE2"]][[1]][["ids"]]
table(ids_only_DE2_mu8_p1_2$xgr[ids_only_DE2_mu8_p1_2$Zg2 == 1])
ids_only_DE2_mu8_p2_2 <- all_runs[[1]][[1]][[1]][[1]]$pow2[["only_DE2"]][[1]][["ids"]]
table(ids_only_DE2_mu8_p2_2$xgr[ids_only_DE2_mu8_p1_2$Zg2 == 1])

# Get the CD and realize that is the TD
all_runs[[1]][[1]][[1]][[1]]$pow1[["only_DE2"]][[1]][["CD"]]

# Analyze the data_sce and try to find if the lfc of the data_sce
# is similar to the lfcass, finds that no, can not get why.
# Same thing happens with pidf and pi.dfass.
# I could understand why it could be different from the strata calculated, because
# the mu and the rate0 simply calculated from the matrix of counts (as is in strata)
# is a little different from the mu and pi.g calculated from Est_Paras (because there are
# some adjustment of the structure of the data, etc.)
lfc <- all_runs[[1]][[1]][[1]][[1]]$data_sce$lfc
pidf <- all_runs[[1]][[1]][[1]][[1]]$data_sce$pi.df
all(lfc[lfc!=0] == all_runs[[1]][[1]][[1]][[1]]$lfcass)
all(pidf[pidf!=0] == all_runs[[1]][[1]][[1]][[1]]$pi.dfass)

# Get the date of the run so I can save properly the plots
date <- "220202"

# Get the data.frame of the lfc, pi.g of celltype1, celltype2,
# and the plot the histogram of the analyzed genes per lfc.
# Finds that the mostly analyzed genes have zero lfc
# and there are still genes with lfc near 1 and 8, which is fine.
# It is logical that the most genes are with 0 lfc, because only
# a 10% of the genes are differetially expressed.
# There is sligthly more genes with positive lfc.
ret_data_sce <- do.call(rbind, map2(all_runs, Details$Round, get_data_sce))
ret_data_sce$abs.lfc <- as.factor(ret_data_sce$abs.lfc)
p <- ggplot(ret_data_sce, aes(x = lfc, fill = abs.lfc)) +
  geom_histogram() + 
  facet_wrap(vars(abs.lfc))
ggsave(paste0("Outputs/Histogram_analyzed_genes_per_lfc_value_", date, ".png"), p)
# Just to take into acount, the lfc is not exactly the lfc of the matrix
# because the mu of the lfc is calculated without zeros, logaritmized and divided by the logaritmized library size

# Get the genes per strata  using the dataframe data_sce, but remember that the lfc and rate0
# of the data_sce is different from de lfcass
ret_data_sce$xgr <- cut(abs(ret_data_sce$lfc), c(0, 2, 4, 8, 16, Inf))
ret_data_sce$rate0 <- cut((1-rowMeans(cbind(ret_data_sce$pi.g1, ret_data_sce$pi.g2))), seq(0, 1, by = 0.2))
count_strata <- filter(ret_data_sce, abs(lfc)>0 & pi.df == 0) %>% group_by(abs.lfc, Round) %>% 
  summarize(counts = table(xgr), Strata = names(table(xgr))) %>% group_by(Strata, abs.lfc) %>%
  summarize(counts = mean(counts))
count_strata_0 <- filter(ret_data_sce, abs(lfc)>0 & pi.df == 0) %>% group_by(abs.lfc, Round) %>% 
  summarize(counts = table(rate0), Strata = names(table(rate0)))

# Do the same but for run 220201_3 or 220201_2
ret_data_sce_2 <- do.call(rbind, map2(all_runs, "1", get_data_sce))
ret_data_sce_2$xgr <- cut(abs(ret_data_sce_2$abs.lfc), seq(0, 1, by = 0.2))
ret_data_sce_2$rate0 <- cut((1-rowMeans(cbind(ret_data_sce_2$pi.g1, ret_data_sce_2$pi.g2))), seq(0, 1, by = 0.2))
count_strata_2 <- filter(ret_data_sce_2, abs(lfc)>0.5 & pi.df == 0) %>% group_by(abs.lfc, Round) %>% 
  summarize(counts = table(xgr), Strata = names(table(xgr))) %>% group_by(Strata, abs.lfc) %>%
  summarize(counts = mean(counts))
count_strata_0_2 <- filter(ret_data_sce_2, abs(lfc)>0 & pi.df == 0 & id %in% all_runs[[1]][[1]][[1]][[1]]$ix.DE2_only) %>% group_by(abs.lfc, Round) %>% 
  summarize(counts = table(rate0), Strata = names(table(rate0)))

# Get the proportion of genes analized that are only_DE2, only_DE1, DE2, DE1.
get_prop<- function(POWSCobj, name){
  dpi.g <- names(POWSCobj[[1]]) %>% as.numeric()
  abs.lfc <- names(POWSCobj[[1]][[1]]) %>% as.numeric()
  prop_DE1_DE2 <- do.call(rbind, lapply(POWSCobj, function(x){
    ix.onlyDE2 <- x[[1]][[1]]$ix.DE2_only[x[[1]][[1]]$ix.DE2_only %in% x[[1]][[1]]$ix.de]
    ix.DE2 <- union(x[[1]][[1]]$ix.DE2_only, x[[1]][[1]]$ix.DE1_DE2)
    ix.DE2 <- ix.DE2[ix.DE2 %in% x[[1]][[1]]$ix.de]
    ix.onlyDE1 <- x[[1]][[1]]$ix.DE1_only[x[[1]][[1]]$ix.DE1_only %in% x[[1]][[1]]$ix.de]
    ix.DE1 <- union(x[[1]][[1]]$ix.DE1_only, x[[1]][[1]]$ix.DE1_DE2)
    ix.DE1 <- ix.DE1[ix.DE1 %in% x[[1]][[1]]$ix.de]
    ix.DE1_or_DE2 <- union(ix.DE2, ix.DE1)
    ix.DE1_DE2 <- x[[1]][[1]]$ix.DE1_DE2[x[[1]][[1]]$ix.DE1_DE2 %in% x[[1]][[1]]$ix.de]
    prop_ix.DE2only <- length(ix.onlyDE2)/length(ix.DE1_or_DE2)
    prop_ix.DE1only <- length(ix.onlyDE1)/length(ix.DE1_or_DE2)
    prop_ix.DE1_DE2 <- length(ix.DE1_DE2)/length(ix.DE1_or_DE2)
    return(list(prop_ix.DE1only, prop_ix.DE2only, prop_ix.DE1_DE2))
  }))
  Data<- data.frame(prop_ix.DE1only = prop_DE1_DE2[[1]],
                    prop_ix.DE2only = prop_DE1_DE2[[2]],
                    prop_ix.DE1_DE2 = prop_DE1_DE2[[3]],
                    nCells = rep(as.integer(names(POWSCobj)), each = 1),
                    dpi.g = dpi.g, abs.lfc = abs.lfc, Round = as.numeric(name))
  return(Data)
}
prop_general <- do.call(rbind, map2(all_runs, Details$Round, get_prop))
prop_general_mean <- group_by(prop_general, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(prop_ix.DE1only = mean(prop_ix.DE1only), prop_ix.DE2only = mean(prop_ix.DE2only), prop_ix.DE1_DE2 = mean(prop_ix.DE1_DE2))
# Save prop_general_mean:
write_csv(prop_general_mean, paste0("Outputs/prop_general_mean_", date, ".csv"))
prop_general_mean <- pivot_longer(prop_general_mean, cols =  starts_with("prop"), names_to = "subs", values_to = "prop") %>%
  mutate(subs = gsub("prop_ix.", "", subs))
prop_plot <- ggplot(prop_general_mean, aes(x = subs, y = prop, fill = as.factor(nCells))) + 
  geom_col(position = position_dodge()) +
  facet_grid(rows = vars(abs.lfc), cols = vars(dpi.g)) +
  scale_fill_manual (values = cells_colors) +
  remove_grid +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title = paste0("Proportion of genes analyzed"), x = "Subset of genes", y = "Proportion over all analized", fill = "#Cells")
ggsave(paste0("Outputs/prop_general_mean_", date, ".png"), prop_plot,  width = 16,  height = 10)

#  scale_fill_manual(values = c("250" = "purple",
#                                "500" = "orange",
#                                "1000" = "steelblue",
#                               "2000" = "steelblue")) 
#

# Plot the distribution of aucs
# subs = "allDE1"
# p1(AUCtypeI_mar_all_DE1, subs)
subs = "onlyDE1"
p1(AUCtypeI_mar_only_DE1, subs)
subs = "onlyDE2"
p1(AUCtypeI_mar_only_DE2, subs)
subs = "DE1&DE2"
p1(AUCtypeI_mar_DE1_DE2, subs)
subs = "DE1orDE2"
p1(AUCtypeI_mar_DE1_or_DE2, subs)
subs = "DE1&DE2op"
p1(AUCtypeI_mar_DE1_DE2_op, subs)
subs = "DE1&DE2co"
p1(AUCtypeI_mar_DE1_DE2_co, subs)

# Resume Rounds:
# TypeI:
# POWtypeI_all_DE1 <- group_by(POWtypeI_all_DE1, Strata, nCells, dpi.g, abs.lfc = abs.lfc) %>% summarize(power = mean(power, na.rm =T))
# ALPHAtypeI_all_DE1 <- group_by(ALPHAtypeI_all_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha = mean(alpha, na.rm =T))
POWtypeI_only_DE1 <- group_by(POWtypeI_only_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeI_only_DE1 <- group_by(ALPHAtypeI_only_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeI_only_DE1 <- group_by(CDtypeI_only_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeI_only_DE2 <- group_by(POWtypeI_only_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>%  summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeI_only_DE2 <- group_by(ALPHAtypeI_only_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeI_only_DE2 <- group_by(CDtypeI_only_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeI_DE1_DE2 <- group_by(POWtypeI_DE1_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeI_DE1_DE2 <- group_by(ALPHAtypeI_DE1_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeI_DE1_DE2 <- group_by(CDtypeI_DE1_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeI_DE1_or_DE2 <- group_by(POWtypeI_DE1_or_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeI_DE1_or_DE2 <- group_by(ALPHAtypeI_DE1_or_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeI_DE1_or_DE2 <- group_by(CDtypeI_DE1_or_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeI_DE1_DE2_op <- group_by(POWtypeI_DE1_DE2_op, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeI_DE1_DE2_op <- group_by(ALPHAtypeI_DE1_DE2_op, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeI_DE1_DE2_op <- group_by(CDtypeI_DE1_DE2_op, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeI_DE1_DE2_co <- group_by(POWtypeI_DE1_DE2_co, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeI_DE1_DE2_co <- group_by(ALPHAtypeI_DE1_DE2_co, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeI_DE1_DE2_co <- group_by(CDtypeI_DE1_DE2_co, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
# TypeII:
# POWtypeII_all_DE1 <- group_by(POWtypeII_all_DE1, Strata, nCells, dpi.g, abs.lfc = abs.lfc) %>% summarize(power = mean(power, na.rm =T))
# ALPHAtypeII_all_DE1 <- group_by(ALPHAtypeII_all_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha = mean(alpha, na.rm =T))
POWtypeII_only_DE1 <- group_by(POWtypeII_only_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeII_only_DE1 <- group_by(ALPHAtypeII_only_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeII_only_DE1 <- group_by(CDtypeII_only_DE1, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeII_only_DE2 <- group_by(POWtypeII_only_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>%  summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeII_only_DE2 <- group_by(ALPHAtypeII_only_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeII_only_DE2 <- group_by(CDtypeII_only_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeII_DE1_DE2 <- group_by(POWtypeII_DE1_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeII_DE1_DE2 <- group_by(ALPHAtypeII_DE1_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeII_DE1_DE2 <- group_by(CDtypeII_DE1_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeII_DE1_or_DE2 <- group_by(POWtypeII_DE1_or_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeII_DE1_or_DE2 <- group_by(ALPHAtypeII_DE1_or_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeII_DE1_or_DE2 <- group_by(CDtypeII_DE1_or_DE2, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeII_DE1_DE2_op <- group_by(POWtypeII_DE1_DE2_op, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeII_DE1_DE2_op <- group_by(ALPHAtypeII_DE1_DE2_op, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeII_DE1_DE2_op <- group_by(CDtypeII_DE1_DE2_op, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))
POWtypeII_DE1_DE2_co <- group_by(POWtypeII_DE1_DE2_co, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.fdr = mean(power.fdr, na.rm =T), power.pval = mean(power.pval, na.rm =T))
ALPHAtypeII_DE1_DE2_co <- group_by(ALPHAtypeII_DE1_DE2_co, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.fdr = mean(alpha.fdr, na.rm =T), alpha.pval = mean(alpha.pval, na.rm =T))
CDtypeII_DE1_DE2_co <- group_by(CDtypeII_DE1_DE2_co, Strata, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(TD.fdr = mean(TD.fdr, na.rm =T), TD.pval = mean(TD.pval, na.rm =T))

# Obtain the marginal power, alpha and auc:
# POWtypeI_mar_all_DE1 <- group_by(POWtypeI_mar_all_DE1, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.marginal = mean(power.marginal, na.rm =T))
# ALPHAtypeI_mar_all_DE1 <- group_by(ALPHAtypeI_mar_all_DE1, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.marginal = mean(alpha.marginal, na.rm =T))
# AUCtypeI_mar_all_DE1 <- group_by(AUCtypeI_mar_all_DE1, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(auc = mean(auc, na.rm =T))
POWtypeI_mar_only_DE1 <- group_by(POWtypeI_mar_only_DE1, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.marginal.fdr = mean(power.marginal.fdr, na.rm =T), power.marginal.pval = mean(power.marginal.pval, na.rm =T))
ALPHAtypeI_mar_only_DE1 <- group_by(ALPHAtypeI_mar_only_DE1, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(alpha.marginal.fdr = mean(alpha.marginal.fdr, na.rm =T), alpha.marginal.pval = mean(alpha.marginal.pval, na.rm =T))
AUCtypeI_mar_only_DE1 <- group_by(AUCtypeI_mar_only_DE1, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(auc.fdr = mean(auc.fdr, na.rm =T), auc.pval = mean(auc.pval, na.rm =T))
POWtypeI_mar_only_DE2 <- group_by(POWtypeI_mar_only_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.marginal.fdr = mean(power.marginal.fdr, na.rm =T), power.marginal.pval = mean(power.marginal.pval, na.rm =T))
ALPHAtypeI_mar_only_DE2 <- group_by(ALPHAtypeI_mar_only_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>%summarize(alpha.marginal.fdr = mean(alpha.marginal.fdr, na.rm =T), alpha.marginal.pval = mean(alpha.marginal.pval, na.rm =T))
AUCtypeI_mar_only_DE2 <- group_by(AUCtypeI_mar_only_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(auc.fdr = mean(auc.fdr, na.rm =T), auc.pval = mean(auc.pval, na.rm =T))
POWtypeI_mar_DE1_DE2 <- group_by(POWtypeI_mar_DE1_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.marginal.fdr = mean(power.marginal.fdr, na.rm =T), power.marginal.pval = mean(power.marginal.pval, na.rm =T))
ALPHAtypeI_mar_DE1_DE2 <- group_by(ALPHAtypeI_mar_DE1_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>%  summarize(alpha.marginal.fdr = mean(alpha.marginal.fdr, na.rm =T), alpha.marginal.pval = mean(alpha.marginal.pval, na.rm =T))
AUCtypeI_mar_DE1_DE2 <- group_by(AUCtypeI_mar_DE1_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(auc.fdr = mean(auc.fdr, na.rm =T), auc.pval = mean(auc.pval, na.rm =T))
POWtypeI_mar_DE1_or_DE2 <- group_by(POWtypeI_mar_DE1_or_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.marginal.fdr = mean(power.marginal.fdr, na.rm =T), power.marginal.pval = mean(power.marginal.pval, na.rm =T))
ALPHAtypeI_mar_DE1_or_DE2 <- group_by(ALPHAtypeI_mar_DE1_or_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>%  summarize(alpha.marginal.fdr = mean(alpha.marginal.fdr, na.rm =T), alpha.marginal.pval = mean(alpha.marginal.pval, na.rm =T))
AUCtypeI_mar_DE1_or_DE2 <- group_by(AUCtypeI_mar_DE1_or_DE2, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(auc.fdr = mean(auc.fdr, na.rm =T), auc.pval = mean(auc.pval, na.rm =T))
POWtypeI_mar_DE1_DE2_op <- group_by(POWtypeI_mar_DE1_DE2_op, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.marginal.fdr = mean(power.marginal.fdr, na.rm =T), power.marginal.pval = mean(power.marginal.pval, na.rm =T))
ALPHAtypeI_mar_DE1_DE2_op <- group_by(ALPHAtypeI_mar_DE1_DE2_op, nCells, dpi.g,  abs.lfc = abs.lfc) %>%  summarize(alpha.marginal.fdr = mean(alpha.marginal.fdr, na.rm =T), alpha.marginal.pval = mean(alpha.marginal.pval, na.rm =T))
AUCtypeI_mar_DE1_DE2_op <- group_by(AUCtypeI_mar_DE1_DE2_op, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(auc.fdr = mean(auc.fdr, na.rm =T), auc.pval = mean(auc.pval, na.rm =T))
POWtypeI_mar_DE1_DE2_co <- group_by(POWtypeI_mar_DE1_DE2_co, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(power.marginal.fdr = mean(power.marginal.fdr, na.rm =T), power.marginal.pval = mean(power.marginal.pval, na.rm =T))
ALPHAtypeI_mar_DE1_DE2_co <- group_by(ALPHAtypeI_mar_DE1_DE2_co, nCells, dpi.g,  abs.lfc = abs.lfc) %>%  summarize(alpha.marginal.fdr = mean(alpha.marginal.fdr, na.rm =T), alpha.marginal.pval = mean(alpha.marginal.pval, na.rm =T))
AUCtypeI_mar_DE1_DE2_co <- group_by(AUCtypeI_mar_DE1_DE2_co, nCells, dpi.g,  abs.lfc = abs.lfc) %>% summarize(auc.fdr = mean(auc.fdr, na.rm =T), auc.pval = mean(auc.pval, na.rm =T))

# Ploting variables:
p2 <- function(POWtypeI_mar, subs){
  p <- ggplot(POWtypeI_mar, aes(x = as.factor(nCells), y = as.factor(dpi.g), fill = power.marginal.fdr)) + geom_tile() +
  scale_fill_continuous(type = "viridis", direction = -1) +
  labs(title = paste0("Genes - ", subs), x = "# Total Cells (both groups summed)", y = "Delta pi.g") +
  my_theme +
  remove_grid +
  facet_wrap(vars(abs.lfc), ncol = 2, labeller = label_both)
ggsave(paste0("Outputs/POWSC_power_", subs, "_", date, ".png"), p, width = 8, height = 8)
}

p3 <- function(AUCtypeI_mar, subs){
  p <- ggplot(AUCtypeI_mar, aes(x = as.factor(nCells), y = as.factor(dpi.g), fill = auc.fdr)) + geom_tile() +
    scale_fill_continuous(type = "viridis", direction = -1) +
    labs(title = paste0("Genes - ", subs), x = "# Total Cells (both groups summed)", y = "Delta pi.g") +
    my_theme +
    remove_grid +
    facet_wrap(vars(abs.lfc), ncol = 2, labeller = label_both)
  ggsave(paste0("Outputs/POWSC_auc_", subs, "_", date, ".png"), p, width = 8, height = 8)
}

# Ploting variables per Strata
p4 <- function(POWtypeI, ALPHAtypeI, subs, strata, tmpxlab_I, type){
  tit_I = ggtitle(paste0("Genes ", subs , " - ", type, " - ", date))
  breaks_I = round(seq(0, 1, length = 6), 1)
  p1 <- POWtypeI %>%
    mutate(nCells = as.factor(nCells), muLFC = as.factor(abs.lfc), dpig = as.factor(dpi.g)) %>%
    plot_function(tit_I, tmpxlab_I, breaks_I, "power.fdr", strata) +
    facet_wrap(vars(abs.lfc), ncol = 2, labeller = label_both)
  p2 <- ALPHAtypeI %>%
    mutate(nCells = as.factor(nCells), muLFC = as.factor(abs.lfc), dpig = as.factor(dpi.g)) %>%
    plot_function(tit_I, tmpxlab_I, breaks_I, "alpha.fdr", strata) +
    facet_wrap(vars(abs.lfc), ncol = 2, labeller = label_both)
  pLegendStrataCells <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none", plot.title = element_blank())
  pStrataCells <- grid.arrange(p1, pLegendStrataCells, p2, layout_matrix = rbind(c(1,2), c(3,2)), nrow = 2, widths = c(4.4, 1), heights = c(2.5, 2.5))
  ggsave(paste0("Outputs/POWSC_StrataCells_", subs, "_", type, "_", date, ".png"), pStrataCells, width = 12, height = 10)
}

# Ploting CD variables per Strata
p5 <- function(POWtypeI, CDtypeI, subs, strata, tmpxlab_I, type){
  tit_I = ggtitle(paste0("Genes ", subs , " - ", type, " - ", date))
  breaks_I = round(seq(0, 1, length = 6), 1)
  typeI <- left_join(POWtypeI, CDtypeI)
  typeI <- typeI %>% filter(dpi.g == 0.4) %>%
    mutate(nCells = as.factor(nCells), muLFC = as.factor(abs.lfc)) %>%
    plot_function2(tit_I, tmpxlab_I, breaks_I, "power.fdr", strata) +
    facet_wrap(vars(abs.lfc), ncol = 2, labeller = label_both)
  ggsave(paste0("Outputs/POWSC_StrataCells_CD_", subs, "_", type, "_", date, ".png"), typeI, width = 12, height = 10)
}


# Ploting variables, marginal and per strata
# subs = "allDE1"
# p2(POWtypeI_mar_all_DE1, subs)
# p3(AUCtypeI_mar_all_DE1, subs)
# p4(POWtypeI_all_DE1, ALPHAtypeI_all_DE1, subs)
tmpxlab_I = "strata of zero ratios"
subs = "onlyDE1"
p2(POWtypeI_mar_only_DE1, subs)
p3(AUCtypeI_mar_only_DE1, subs)
p4(POWtypeI_only_DE1, ALPHAtypeI_only_DE1, subs, strata = "pi", tmpxlab_I, "0ratio")
p5(POWtypeI_only_DE1, CDtypeI_only_DE1, subs, strata = "pi", tmpxlab_I, "0ratio")
subs = "onlyDE2"
p2(POWtypeI_mar_only_DE2, subs)
p3(AUCtypeI_mar_only_DE2, subs)
p4(POWtypeI_only_DE2, ALPHAtypeI_only_DE2, subs, strata = "pi", tmpxlab_I, "0ratio")
p5(POWtypeI_only_DE2, CDtypeI_only_DE2, subs, strata = "pi", tmpxlab_I, "0ratio")
subs = "DE1&DE2"
p2(POWtypeI_mar_DE1_DE2, subs)
p3(AUCtypeI_mar_DE1_DE2, subs)
p4(POWtypeI_DE1_DE2, ALPHAtypeI_DE1_DE2, subs, strata = "pi", tmpxlab_I, "0ratio")
p5(POWtypeI_DE1_DE2, CDtypeI_DE1_DE2, subs, strata = "pi", tmpxlab_I, "0ratio")
subs = "DE1orDE2"
p2(POWtypeI_mar_DE1_or_DE2, subs)
p3(AUCtypeI_mar_DE1_or_DE2, subs)
p4(POWtypeI_DE1_or_DE2, ALPHAtypeI_DE1_or_DE2, subs, strata = "pi", tmpxlab_I, "0ratio")
p5(POWtypeI_DE1_or_DE2, CDtypeI_DE1_or_DE2, subs, strata = "pi", tmpxlab_I, "0ratio")
subs = "DE1&DE2op"
p2(POWtypeI_mar_DE1_DE2_op, subs)
p3(AUCtypeI_mar_DE1_DE2_op, subs)
p4(POWtypeI_DE1_DE2_op, ALPHAtypeI_DE1_DE2_op, subs, strata = "pi", tmpxlab_I, "0ratio")
p5(POWtypeI_DE1_DE2_op, CDtypeI_DE1_DE2_op, subs, strata = "pi", tmpxlab_I, "0ratio")
subs = "DE1&DE2co"
p2(POWtypeI_mar_DE1_DE2_co, subs)
p3(AUCtypeI_mar_DE1_DE2_co, subs)
p4(POWtypeI_DE1_DE2_co, ALPHAtypeI_DE1_DE2_co, subs, strata = "pi", tmpxlab_I, "0ratio")
p5(POWtypeI_DE1_DE2_co, CDtypeI_DE1_DE2_co, subs, strata = "pi", tmpxlab_I, "0ratio")
#Pow2
tmpxlab_II = "strata of the mean expression"
subs = "onlyDE1_mean"
p4(POWtypeII_only_DE1, ALPHAtypeII_only_DE1, subs, strata = "mu", tmpxlab_II, "MeanExp")
p5(POWtypeII_only_DE1, CDtypeII_only_DE1, subs, strata = "mu", tmpxlab_II, "MeanExp")
subs = "onlyDE2_mean"
p4(POWtypeII_only_DE2, ALPHAtypeII_only_DE2, subs, strata = "mu", tmpxlab_II, "MeanExp")
p5(POWtypeII_only_DE2, CDtypeII_only_DE2, subs, strata = "mu", tmpxlab_II, "MeanExp")
subs = "DE1&DE2_mean"
p4(POWtypeII_DE1_DE2, ALPHAtypeII_DE1_DE2, subs, strata = "mu", tmpxlab_II, "MeanExp")
p5(POWtypeII_DE1_DE2, CDtypeII_DE1_DE2, subs, strata = "mu", tmpxlab_II, "MeanExp")
subs = "DE1orDE2_mean"
p4(POWtypeII_DE1_or_DE2, ALPHAtypeII_DE1_or_DE2, subs, strata = "mu", tmpxlab_II, "MeanExp")
p5(POWtypeII_DE1_or_DE2, CDtypeII_DE1_or_DE2, subs, strata = "mu", tmpxlab_II, "MeanExp")
subs = "DE1&DE2op_mean"
p4(POWtypeII_DE1_DE2_op, ALPHAtypeII_DE1_DE2_op, subs, strata = "mu", tmpxlab_II, "MeanExp")
p5(POWtypeII_DE1_DE2_op, CDtypeII_DE1_DE2_op, subs, strata = "mu", tmpxlab_II, "MeanExp")
subs = "DE1&DE2co_mean"
p4(POWtypeII_DE1_DE2_co, ALPHAtypeII_DE1_DE2_co, subs, strata = "mu", tmpxlab_II, "MeanExp")
p5(POWtypeII_DE1_DE2_co, CDtypeII_DE1_DE2_co, subs, strata = "mu", tmpxlab_II, "MeanExp")

# Plot ROC  
pow_rslt_t_DESeq <- readRDS("Outputs/pow_scDESeq_220112_size_250_mu_1_dpig_0.1_round_1.rds")
sens.ci <- pow_rslt_t_DESeq[[1]][[1]][["pow1"]][["roc"]] %>%
  ci.se(specificities = c(seq(0, 0.9, by = 0.1), 0.93, 0.97, 0.95, 0.97, 0.98, 0.99, 1))
plot.roc(pow_rslt_t_DESeq[[1]][[1]][["pow1"]][["roc"]], print.thres = T,
         print.auc = T, auc.polygon = T, ci = T, ci.type = "bars", xlim = c(1, 0))
plot(sens.ci, type = "shape", col = "lightgray")
plot(sens.ci, type = "bars")

plotroc <- ggroc(pow_rslt_t_DESeq[[1]][[1]][["pow1"]][["roc"]]) +
  scale_x_reverse(breaks = seq(1, 0, by = -0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
roc_thr <- coords(pow_rslt_t_DESeq[[1]][[1]][["pow1"]][["roc"]], x = "best")
roc_0.1 <- coords(pow_rslt_t_DESeq[[1]][[1]][["pow1"]][["roc"]], x = (1 - 0.1), input = "threshold")
roc_0.05 <- coords(pow_rslt_t_DESeq[[1]][[1]][["pow1"]][["roc"]], x = (1 - 0.05), input = "threshold")
plotroc +
  geom_vline(xintercept = roc_thr$specificity, linetype = "dashed", colour = "gray85") +
  geom_point(aes(x = roc_thr$specificity, y = roc_thr$sensitivity), colour = "gray85") +
  geom_text(aes(x = roc_thr$specificity, y = roc_thr$sensitivity, label = round(1 - roc_thr$threshold, 2)), hjust = -0.1, vjust = -1.2, colour = "gray85") +
  geom_vline(xintercept = roc_0.1$specificity, linetype = "dashed", colour = "gray70") +
  geom_point(aes(x = roc_0.1$specificity, y = roc_0.1$sensitivity), colour = "gray70") +
  geom_text(aes(x = roc_0.1$specificity, y = roc_0.1$sensitivity, label = round(1 - roc_0.1$threshold, 2)), hjust = -0.1, vjust = -1.2, colour = "gray70") +
  geom_vline(xintercept = roc_0.05$specificity, linetype = "dashed", colour = "gray47") +
  geom_point(aes(x = roc_0.05$specificity, y = roc_0.05$sensitivity), colour = "gray47") +
  geom_text(aes(x = roc_0.05$specificity, y = roc_0.05$sensitivity, label = round(1 - roc_0.05$threshold, 2)), hjust = -0.1, vjust = -1.2, colour = "gray47") +
  my_theme +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
