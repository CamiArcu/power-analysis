runPOWSC_genefilter_furrr <- function (sim_size, mu, sdP, dpi.g, per_DE = 0.05, 
                                       est_Paras, DE_Method = c("MAST", "SC2P", "DESeq2", "ScDESeq2",  "SeuratDESeq2", "DESeq2zw"), Cell_Type = c("PW",                                                                                                     "Multi"), multi_Prob = NULL, alpha = 0.1, disc_delta = 0.1, 
                                       cont_delta = 0.5) 
{
  DE_Method = match.arg(DE_Method)
  Cell_Type = match.arg(Cell_Type)
  pow_rslt = NULL
  if (Cell_Type == "PW") {
    print(paste0("Running simulation with ", sim_size, " cells and mu = ", mu, ", sd = ", sdP, " and dpi.g = ", dpi.g))
    simData = Simulate2SCE_explicada_Zex2(mu = mu, sdP = sdP, dpi.g = dpi.g, n=sim_size, estParas1 = est_Paras, estParas2 = est_Paras)
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
                    ix.DE1_only = simData$ix.DE1_only, ix.DE2_only = simData$ix.DE2_only, ix.de = de$cont$geneIndex) # pow2 = pow2_rslt
    pow_rslt[[toString(dpi.g)]][[toString(mu)]] = tmp_rslt
    print("llegue hasta aca")
  }
  class(pow_rslt) = "POWSC"
  return(pow_rslt)
}

Simulate2SCE_explicada_Zex2 <- function (n = 100, perDE = 0.05, estParas1, estParas2, mu, sdP, dpi.g) 
{
  start <- Sys.time()
  # Divido el dataset en 2
  n1 = n2 = round(n/2)
  # Me quedo con el el parametro pi.g estimado para el grupo 1, que es el numero de celulas con transcripcion activa 
  pi.g1 = estParas1$pi.g
  # Me quedo con el el parametro pi.g estimado para el grupo 2, que es el numero de celulas con transcripcion activa 
  pi.g2 = estParas2$pi.g
  # # es igual a p0 de la ZIP, la proporcion estimada de genes con counts = 0, punto de masa cero.
  p01 = estParas1$p0
  p02 = estParas2$p0
  # lambda es el parametro de la distribucion ZIP, para transcripcion inactiva. Estimo este parametro para ambos grupos
  lambda1 = estParas1$lambda
  lambda2 = estParas2$lambda
  # Me quedo con la media estimada para cada grupo que corresponde a los genes de la transcripcion activa
  mu1 = estParas1$mu
  mu2 = estParas2$mu
  # Me quedo con el desvio estadar estimado para cada grupo que corresponde a los genes de la transcripcion activa
  sigma1 = estParas1$sd
  sigma2 = estParas2$sd
  # Creo que es la media de expresion por celula de los genes que superan la media estimada
  sf1 = estParas1$sf
  sf2 = estParas2$sf
  if (length(pi.g1) == length(pi.g2)) {
    ngene = length(pi.g1)
  }
  else {
    stop("Please Input Proper Parameter Estimation Objects")
  }
  # Me fijo cual grupo tiene mas medias mayores que 3 (Creo que se determina por default que una media menor o igual a 3, corresponde a un gen de baja expresion)
  # De este grupo, me quedo con el numero de genes cuyas medias son mayor a 3, osea, que estan expresandose activamente.
  n0 = max(sum(mu1 > 3), sum(mu2 > 3))
  # Determino el numero de genes DE en cada grupo
  nDE1 = nDE2 = ngene * perDE
  # Ordeno los genes de acuerdo a su media de expresion de mayor a menor valor y me quedo con los genes expresandose activamente
  ix1.highGenes = order(mu1, decreasing = TRUE)[seq_len(n0)]
  ix2.highGenes = order(mu2, decreasing = TRUE)[seq_len(n0)]
  # El cambio que introduzco para que los ix.DE1 no sean ix.DE2:
  ix.high = union(ix1.highGenes, ix2.highGenes)
  # De los genes expresados activamente, tomo al azar nDE1 genes del total para ser genes DE del grupo 1
  # que van a ser los genes que hagan fase transition
  ix.DE1 = sample(ix.high, nDE1)
  # Hago lo mismo para el grupo 2, que van a ser los genes que hagan transcription change
  ix.DE2 = sample(ix.high, nDE2) # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  # Me quedo con los genes DE tabto del grupo 1 como del 2
  ix.DEs = union(ix.DE1, ix.DE2)
  # Tomo al azar expresiones medias que superan la media estimada, tanto para las celulas del grupo 1
  sf1 = sample(sf1, n1, replace = TRUE)
  # Como para las celulas del grupo 2
  sf2 = sample(sf2, n2, replace = TRUE)
  # Tomo al azar proporciones estimadas de genes con counts = 0 para cada grupo (se pueden repetir)
  p0.1 = sample(p01, n1, replace = TRUE)
  p0.2 = sample(p02, n2, replace = TRUE)
  # Tomo al azar lambdas para cada gen para cada grupo (se pueden repetir)
  lambda1 = sample(lambda1, n1, replace = TRUE)
  lambda2 = sample(lambda2, n2, replace = TRUE)
  # Tomo la proporcion de celulas con transcripcion activa del grupo 2 de los genes que haran phase transition.
  tmp = pi.g2[ix.DE1]
  # Genero los desvios dpi.g
  # Para evitar que los valores de pi.g2 superen a 1 o sean menores que 0, defino primero que genes
  # tienen incialmente un pi.g que se encuentren entre un rango de 0 al 1-dpi.g y de dpi.g a 1.
  # Ojo, para que esto funcione el dpi.g tiene que estar entre (0 y 0.5) (sin incluir los extremos).
  sum_tmp <- ((tmp >= 0) & (tmp <= (1-dpi.g))) #Los genes que entren en este rango se les va a poder sumar dpi.g
  sub_tmp <- ((tmp >= dpi.g) & (tmp <=1)) #Los genes que entren en este rango se les va a poder restar dpi.g
  # sin embargo estos rangos se superponen, entonces a los genes que entren en ambos rangos se les puede tanto sumar como restar dpi.g
  # Entonces lo primero que hago es sumarle dpi.g a los genes que entren solo en el primer rango (osea, que vayan de 0 a dpi.g)
  tmp[sum_tmp & (sum_tmp != sub_tmp)] <- tmp[sum_tmp & (sum_tmp != sub_tmp)] + dpi.g
  # Luego le resto dpi.g a los genes que entren solo en el segundo rango (osea, que vayan de 1-dpi.g a 1)
  tmp[sub_tmp & (sum_tmp != sub_tmp)] <- tmp[sub_tmp & (sum_tmp != sub_tmp)] - dpi.g
  # Luego extraigo el indice de los genes que  entran en ambos rangos (osea que van de dpi.g a 1-dpig)
  ix.subsum_tmp <- (1:nDE1)[(sum_tmp == sub_tmp)]
  # Me fijo cuantos genes son
  n.subsum_tmp <- length(ix.subsum_tmp)
  # A la mitad de los genes que entran en ambos rangos les sumo dpi.g
  ix.subsum_tmp_sum<- sample(ix.subsum_tmp, (n.subsum_tmp%/%2))
  tmp[ix.subsum_tmp_sum] <- tmp[ix.subsum_tmp_sum] + dpi.g
  # A la otra mitad de los genes que entran en ambos rangos les resto dpi.g
  ix.subsum_tmp_sub <- ix.subsum_tmp[!(ix.subsum_tmp %in% ix.subsum_tmp_sum)]
  tmp[ix.subsum_tmp_sub] <- tmp[ix.subsum_tmp_sub] - dpi.g
  # Para la proporcion de celulas con transc activa mayor al 50%, entonces, modifico esta proporcion, restandole
  # valores uniformes, aleatoreos, entre 0.1 y 0.3
  # Agrego estas modificaciones, creo que la gracia de este paso es homogeneizar las proporciones celulares con transcripcion activa para los genes de fase 1
  pi.g2[ix.DE1] = tmp
  # Las diferencias entre las proporciones del grupo 1 y el grupo 2.
  #   pi.dfass = pi.g2[ix.DE1] - pi.g1[ix.DE1]
  # Reemplazo tmp por una mezcla de valores provenientes de 2 distribuciones normales que representaran los log2FC
  # una con media -1 y otra con media 1
  #tmp = c(rgamma(1000, shape = 4, rate = 4/(-mu)), rgamma(1000, shape = 4, rate = 4/(mu)))
  # tmp = c(rnorm(1000, mean = -mu, sd = sdP), rnorm(1000, mean = mu, sd = sdP))
  # entonces selecciono logFC de esta mezcla de valores
  magnitude = sample(c(-1,1), length(ix.DE2), replace = T)
  # Genero las media del grupo 2 como la media del gripo 1, con este determinado cambio
  mu2[ix.DE2] = mu1[ix.DE2] + mu * magnitude
  y1 = GenerateCountMatrix_Zex2(pi.g1, p0.1, lambda1, mu1, sigma1, sf1)
  y2 = GenerateCountMatrix_Zex2(pi.g2, p0.2, lambda2, mu2, sigma2, sf2)
  y = cbind(y1, y2)
  rownames(y) = paste0("g", seq_len(nrow(y)))
  celltypes = rep(paste0("celltype", c(1, 2)), c(n1, n2))
  sce = SingleCellExperiment(assays = list(counts = y), rowData = data.frame(geneNames = rownames(y), stringsAsFactors = FALSE), colData = data.frame(cellTypes = celltypes, stringsAsFactors = FALSE))
  # ALGUNAS PRUEBAS PARA VER QUE EL MODELO COINCIDE CON LAS ASIGNACIONES
  # L1 = colSums(y1[ix.DE2,] * (y1[ix.DE2,]>0))/1e+06
  # mu.g1 = log2(rowSums(y1[ix.DE2,] * (y1[ix.DE2,]>0))/rowSums(sweep(y1[ix.DE2,]>0, 2, L1, FUN = "*")))
  # L2 = colSums(y2[ix.DE2,] * (y2[ix.DE2,]>0))/1e+06
  # mu.g2 = log2(rowSums(y2[ix.DE2,] * (y2[ix.DE2,]>0))/rowSums(sweep(y2[ix.DE2,]>0, 2, L2, FUN = "*")))
  # lfc = mu.g2 - mu.g1
  # DE2_and_changeMu  <-    (lfc >= lfcass-0.2 | lfc <= lfcass+0.2) # No voy a tener una lfc exactamente igual al asignado porque depende tambien si hay valores muy extremos (que superan al maximo de enteros entonces estos se llevan al limite maximo) 
  # Además de que podrian haber otras variaciones sutiles. ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  # ix.DE2real <- ix.DE2[DE2_and_changeMu] # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  # pi.df = rowMeans(y2[ix.DE1,] > 2) - rowMeans(y1[ix.DE1,] > 2) # Pongo 2 porque es el numero de counts promedio de una transcripcion inactiva que no es cero. Ver est_Param$lambda
  # DE1_and_changePig  <- ix.DE1 %in% which(abs(pi.df) >= dpi.g-0.1 | abs(pi.df) <= dpi.g+0.1 ) # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  # ix.DE1real <- ix.DE1[DE1_and_change       Pig] # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  # is.DE1_DE2_real = ix.DE1[ix.DE1 %in% ix.DE2real] # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  # LO MODELADO COINCIDE CON LAS ASIGNACIONES!! QUE FELICIDADDDDDD
  ix.DE1_or_DE2 = union(ix.DE1, ix.DE2) # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  ix.DE1_DE2 = ix.DE1[ix.DE1 %in% ix.DE2] # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  ix.DE1_only = ix.DE1[!(ix.DE1 %in% ix.DE2)] # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  ix.DE2_only = ix.DE2[!(ix.DE2 %in% ix.DE1)] # ESTA ES LA LINEA QUE ESTOY CAMBIANDO RESPECTO A runPOWSC_genefilter_furrr.R
  data_sce <- data.frame(id = seq_len(nrow(y)), lfc = mu2-mu1, pi.df = pi.g2 - pi.g1,
                         mu1 = mu1, mu2 = mu2, pi.g1 = pi.g1, pi.g2 = pi.g2)
  print("ya pase por Simulate2SCE")
  end <- Sys.time() 
  print(paste0("Time of Simulate2SCE_explicada_Zex2: ", round(end - start, 0), " ", units(end - start)))
  list(ix.DE1 = ix.DE1, ix.DE2 = ix.DE2, ix.DE1_or_DE2 = ix.DE1_or_DE2, data_sce = data_sce,
       sce = sce, ngenes = ngene, 
       ix.DE1_DE2 = ix.DE1_DE2, ix.DE1_only= ix.DE1_only, ix.DE2_only = ix.DE2_only) 
  # ix.DE2real = ix.DE2real, ix.DE2ass = ix.DE2, ix.DE1_DE2_real = ix.DE1_DE2_real, pi.df = pi.df
}

runDE_explicada<- function (sce, DE_Method = c("MAST", "SC2P", "DESeq2", "ScDESeq2",  "SeuratDESeq2", "DESeq2zw")) 
{
  DE_Method <- match.arg(DE_Method)
  if (DE_Method == "MAST") {
    DE_rslt <- runMAST(sce)
  }
  else if(DE_Method == "SC2P") {
    DE_rslt <- runSC2P(sce)
  }
  else if(DE_Method == "DESeq2"){
    DE_rslt <- runDESeq2(sce)
  }
  else if(DE_Method == "ScDESeq2"){
    DE_rslt <- runScDESeq2(sce)
  }
  else if(DE_Method == "SeuratDESeq2"){
    DE_rslt <- runSeuratDESeq2(sce)
  }
  else if(DE_Method == "DESeq2zw"){
    DE_rslt <- runDESeq2zinbwave(sce)
  }
  else{
    return(print(paste0(DE_Method, " is not an available method. Try 'MAST', 'SC2P' or 'DESeq2'")))
  }
  return(DE_rslt)
}

Pow_Disc_different_DEgenes <- function(DErslt, simData, alpha = 0.05, strata = seq(0, 1, by = 0.2)){
  print("pow_diffDEgenes")
  Pow_only_DE1 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE1_only)
  Pow_only_DE2 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE2_only)
  Pow_DE1_DE2 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE1_DE2)
  Pow_DE1_or_DE2 <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = simData$ix.DE1_or_DE2)
  data_sce <- simData$data_sce[simData$ix.DE1_DE2, ]
  ix.DE1_DE2_co <- filter(data_sce, (lfc > 0 & pi.df > 0) | (lfc < 0 & pi.df < 0))$id
  ix.DE1_DE2_op <- filter(data_sce, (lfc > 0 & pi.df < 0) | (lfc < 0 & pi.df > 0))$id
  Pow_DE1_DE2_co <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = ix.DE1_DE2_co)
  if (length(DErslt$cont$geneIndex[DErslt$cont$geneIndex %in% ix.DE1_DE2_op]) == 0){
    Pow_DE1_DE2_op <- "This subset of genes was not considered HVG, therefore, not analyzed"
  } else{
    Pow_DE1_DE2_op <- Pow_Disc_explicada_genefilter(DErslt, simData, alpha, strata, DEid = ix.DE1_DE2_op)
  }
    return(list(only_DE1 = Pow_only_DE1, only_DE2 = Pow_only_DE2, DE1_DE2 = Pow_DE1_DE2, DE1_DE2_co = Pow_DE1_DE2_co, DE1_DE2_op = Pow_DE1_DE2_op, DE1_or_DE2 = Pow_DE1_or_DE2))
}

Pow_Cont_different_DEgenes <- function(DErslt, simData, alpha = 0.05, delta = 0.5, strata = c(0, 2, 4, 8, 16, Inf)){
  print("pow_diffDEgenes")
  Pow_only_DE1 <- Pow_Cont_explicada_genefilter(DErslt, simData, alpha, delta, strata, DEid = simData$ix.DE1_only)
  Pow_only_DE2 <- Pow_Cont_explicada_genefilter(DErslt, simData, alpha, delta, strata, DEid = simData$ix.DE2_only)
  Pow_DE1_DE2 <- Pow_Cont_explicada_genefilter(DErslt, simData, alpha, delta, strata, DEid = simData$ix.DE1_DE2)
  Pow_DE1_or_DE2 <- Pow_Cont_explicada_genefilter(DErslt, simData, alpha, delta, strata, DEid = simData$ix.DE1_or_DE2)
  data_sce <- simData$data_sce[simData$ix.DE1_DE2, ]
  ix.DE1_DE2_co <- filter(data_sce, (lfc > 0 & pi.df > 0) | (lfc < 0 & pi.df < 0))$id
  ix.DE1_DE2_op <- filter(data_sce, (lfc > 0 & pi.df < 0) | (lfc < 0 & pi.df > 0))$id
  Pow_DE1_DE2_co <- Pow_Cont_explicada_genefilter(DErslt, simData, alpha, delta, strata, DEid = ix.DE1_DE2_co)
  if (length(DErslt$cont$geneIndex[DErslt$cont$geneIndex %in% ix.DE1_DE2_op]) == 0){
    Pow_DE1_DE2_op <- "This subset of genes was not considered HVG, therefore, not analyzed"
  } else{
    Pow_DE1_DE2_op <- Pow_Cont_explicada_genefilter(DErslt, simData, alpha, delta, strata, DEid = ix.DE1_DE2_op)
  }
  return(list(only_DE1 = Pow_only_DE1, only_DE2 = Pow_only_DE2, DE1_DE2 = Pow_DE1_DE2, DE1_DE2_co = Pow_DE1_DE2_co, DE1_DE2_op = Pow_DE1_DE2_op, DE1_or_DE2 = Pow_DE1_or_DE2))
}

Pow_Disc_explicada_genefilter<- function (DErslt, simData, alpha, strata, DEid) 
{
  print("pow_disc_new")
  # Extraigo el numero de genes totales
  ngenes <- nrow(simData$sce)
  # Genero dos vectores con 0, de longuitud = numero de genes
  Zg = Zg2 = rep(0, ngenes)
  # A los genes DA (de verdad), les pongo un 1
  Zg[DEid] <- 1
  # Busco el ID de los genes que poseen un pi.df (proporcion de cambio a genes activos) mayor a un umbral (delta)
  # ix <- which(abs(pi.df) > delta) # Comenté esta linea porque con la asignacio de un dpi.g definado, ya no tiene sentido
  # A los genes DA (de verdad), pero que superan el umbral de pi.df les pongo un 1
  Zg2[DEid] <- 1 #Zg2[DEid[ix]] <- 1 #Antes de la asignacion de un dpi.g esta linea era lo comentado
  # Extraigo la matrix de datos simulados
  sce <- simData$sce
  # Extraigo el número de células totales
  ntotal <- ncol(sce)
  # Extraigo la matriz original (de counts) de datos
  Y <- round(assays(sce)[[1]])
  # Para cada gen, calculo la proporción de células en las cuales está inactivo
  rate0 <- rowMeans(Y == 0)
  # De los genes tanto DE1 como DE2, identifico cuales no pertenecen al subset de DEs que quiero analizar ahora
  # Si analizo todos los genes DEs ahora, ix.not = empty
  ix.not <- simData$ix.DE1_or_DE2[!(simData$ix.DE1_or_DE2 %in% DEid)]
  # De los genes analizados, saco los genes que son DEs pero no pertenecen al grupo de DEs analizados ahora
  ix.yes <- DErslt$cont$geneIndex[!(DErslt$cont$geneIndex %in% ix.not)]
  # De los genes totales, me quedo con los finalmente analizados
  rate0 = rate0[ix.yes]
  # Identifico los genes que están inactivos entre 1% y 99% de las celulas (excluyo los puntos extremos, por ej: genes inactivos en todas las célulaso (100%) en ninguna (0%))
  #ix.keep <- intersect(which(rate0 < 0.99), which(rate0 >  0.01)) #esto me va a dar el índice respecto de rate0
  # Divido la proporcion de inactivos de los genes anteriores, en estratos
  xgr <- cut(rate0, strata)#xgr <- cut(rate0[ix.keep], strata)
  #ix.keep <- DErslt$cont$geneIndex[ix.keep] # esto me va a dar el índice respecto de la matriz original
  # De los genes totales, me quedo con los finalmente analizados
  # Tanto para el vector que tiene marcados los genes DA
  Zg <- Zg[ix.yes] 
  # Como el que tiene los gene DA que superan el umbral
  Zg2 = Zg2[ix.yes]
  # Extraigo el FDR del analisis de DE de los datos simulados
  # Esto lo hago con filter porque DErstl es un data frame con un nrows mucho mas chico que Zg2 o rate0, osea, es mas chico el nrows de sce
  # entonces si hiciese DErslt$cont$fdr[ix.yes], estaria agarrando los indices de una tabla de 1000 filas en vez de 10000 filas
  fdrvec <- dplyr::filter(DErslt$cont, geneIndex %in% ix.yes)$fdr
  pvalvec <- dplyr::filter(DErslt$cont, geneIndex %in% ix.yes)$pval
  # Calculo el poder usando el FDR, los vectores que idican los genes DE, los estratos de medias y el alpha (valor de corte)
  pow1 = PowerEst(fdrvec, alpha, Zg, Zg2, xgr = xgr)
  pow2 = PowerEst(pvalvec, alpha, Zg, Zg2, xgr = xgr)
  return(list(pow1, pow2))
}

#c(0, 10, 2^(seq_len(4)) * 10, Inf)
Pow_Cont_explicada_genefilter<- function (DErslt, simData, alpha, delta, strata,  DEid) 
{
  print("pow_cont_new")
  # EXtraigo el log fold-change de la matriz de los datos simulados
  lfc = simData$lfcass
  # Extraigo el numero de genes totales
  ngenes = nrow(simData$sce)
  # Extraigo el índice de los genes asignados para que esten diferencialmente expresados o no se modifiquen
  # DEid = simData$ix.DE2
  # Genero dos vectores con 0, de longuitud = numero de genes
  Zg = Zg2 = rep(0, ngenes)
  # A los genes asignadas para que .. bla bla bla, les pongo un 1
  Zg[DEid] = 1
  # Busco el ID de los genes que poseen un lfc mayor a un umbral (delta) (osea, los genes a los que efectivamente se les cambio la media)
  # ix = which(abs(lfc) > delta) Comento esta linea porque como controlo el valor del lfc, no tiene tanta relevancia esto
  # A los genes DE (de verdad), les pongo un 1 en Zg2
  Zg2[DEid] = 1 # Zg2[DEid[ix]] = 1
  # Extraigo la matrix de datos simulados
  sce = simData$sce
  Y = round(assays(sce)[[1]])
  # Calculo la media de cada gen
  X.bar = rowMeans(Y) #X.bar = abs(lfc)
  # De los genes tanto DE1 como DE2, identifico cuales no pertenecen al subset de DEs que quiero analizar ahora
  # Si analizo todos los genes DEs ahora, ix.not = empty
  ix.not <- simData$ix.DE1_or_DE2[!(simData$ix.DE1_or_DE2 %in% DEid)]
  # De los genes analizados, saco los genes que son DEs pero no pertenecen al grupo de DEs analizados ahora
  ix.yes <- DErslt$cont$geneIndex[!(DErslt$cont$geneIndex %in% ix.not)]
  # De los genes totales, me quedo con los finalmente analizados
  X.bar = X.bar[ix.yes]
  # Extraigo el ID de los genes que tienen una media > 0, osea,que no son 0 para todas las celulas
  # ix.keep = which(X.bar > 0)
  # Divido la medias de expresion de los genes que no son 0 para todas las celulas, en estratos
  xgr = cut(X.bar, strata) # xgr = cut(X.bar[ix.keep], strata)
  # De los genes que no son cero, me quedo con los finalmente analizados
  # ix.keep = DErslt$cont$geneIndex[ix.keep]
  # Me quedo con los genes que no son 0 para todas las celulas y que fueron analizados
  # Tanto para el vector que tiene marcados los genes asignados para que... bla bla bla
  Zg = Zg[ix.yes]
  # Como el que tiene los gene DE (de verdad)
  Zg2 = Zg2[ix.yes]
  # Extraigo el FDR del analisis de DE de los datos simulados
  # Esto lo hago con filter porque DErstl es un data frame con un nrows mucho mas chico que Zg2 o X.bar, osea, es mas chico el nrows de sce
  # entonces si hiciese DErslt$cont$fdr[ix.yes], estaria agarrando los indices de una tabla de 1000 filas en vez de 10000 filas
  fdrvec <- dplyr::filter(DErslt$cont, geneIndex %in% ix.yes)$fdr
  pvalvec <- dplyr::filter(DErslt$cont, geneIndex %in% ix.yes)$pval
  # Calculo el poder usando el FDR, los vectores que idican los genes DE, los estratos de medias y el alpha (valor de corte)
  pow1 = PowerEst(fdrvec, alpha, Zg, Zg2, xgr = xgr)
  pow2 = PowerEst(pvalvec, alpha, Zg, Zg2, xgr = xgr)
  return(list(pow1, pow2))
}


PowerEst <- function (fdr, alpha, Zg, Zg2, xgr) 
{
  # Me quedo con el ID de los genes detectados DE, osea, con FDR menor o igual a alpha
  ix.D = fdr <= alpha
  # Sumo los genes detectados DE
  N = sum(ix.D)
  # Sumos los genes detectados DE por estratos
  # Si no hay genes para un determinado estrato, devuelve NA, no 0
  N.stratified = tapply(ix.D, xgr, sum)
  # Identifico los genes que son realmente DE y superan el umbral de lfc
  id.TP = Zg2 == 1
  # Indentifico, de los genes que son realmente DE y superan el umbral, cuales son detectados DE, osea su FDR es menor o igual a alpha
  # Luego sumo estos genes por estratos
  TD = tapply(fdr[id.TP] <= alpha, xgr[id.TP], sum)
  TD[is.na(TD)] = 0
  # Identifico los genes que realmente no son DE.
  id.TN = Zg == 0
  # Me fijo cuales de estos fueron detectados como DE, osea que su FDR es menor o igual a alpha
  FD = tapply(fdr[id.TN] <= alpha, xgr[id.TN], sum)
  FD[is.na(FD)] = 0
  # Calculo el alpha o el False Posituve Rate (o 1-Especificidad):
  # Divido, por estratos, el numero de genes detectados como DE que realmente no lo son (FP) sobre los genes que realmente no son DE
  alpha = as.vector(FD/table(xgr[id.TN]))
  # Lo mismo, pero todo junto, sin estratos
  alpha.marginal = sum(FD)/sum(id.TN)
  # Divido, por estratos, el numero de genes detectados como DE que realmente lo son (TP) sobre los genes que realmente son DE y superan el umbral
  power = as.vector(TD/table(xgr[id.TP]))
  # Lo mismo, pero todo junto, sin estratos  
  power.marginal = sum(TD, na.rm = TRUE)/sum(id.TP)
  # Divido, por estratos, el numero de FP por el numero de genes DE detectados
  FDR = FD/N.stratified
  # Lo mismo, pero todo junto, sin estratos
  FDR.marginal = sum(FD, na.rm = TRUE)/N
  # Cuento por estratos, el numero de genes asignados como DE
  CT = table(xgr[id.TP])
  # AUC
  AUC <- pROC::auc(response = Zg2, predictor = (1-fdr))
  ROC <- pROC::roc(response = Zg2, predictor = (1-fdr), smoothed = F,
                   ci = T, ci.alpha = 0.9, stratified = F, plot = T,
                   auc.polygon = T, max.auc.polygon = T, grid = T, 
                   print.auc = T, show.thres = T)
  list(CD = TD, FD = FD, TD = CT, alpha.nominal = alpha, alpha = alpha, 
       alpha.marginal = alpha.marginal, power = power, power.marginal = power.marginal, 
       FDR = FDR, FDR.marginal = FDR.marginal, auc = AUC, roc = ROC, ids = data.frame("xgr" = xgr, "Zg2" = Zg2, "fdr" = fdr))
}

rzip = function(n, p0, lambda) {
  y = rep(0, n)
  ix = rbinom(n, size=1, prob=p0)
  y[!ix]  = rpois(sum(!ix), lambda)
  y
}

### generate LNP random variable
## The real data is not exactly LNP, the right tail is much shorter
## Using LNP to simulate there are too many large numbers.
## How to shave the right tail??
rLNP = function(n, mu, sigma, sf) {
  theta = 2^rnorm(n, mean=mu, sd=sigma)
  ## should shave the right tail of theta a bit to avoid excessive large number??
  
  y = suppressWarnings(rpois(n, theta*sf))
  if (sum(is.na(y)) > 0){
    m = max(y[!is.na(y)])
    y[is.na(y)] = m
  }
  return(y)
}

## main function for generating the count matrix
GenerateCountMatrix_Zex2 = function(pi.g, p0, lambda, mu, sigma, sf){
  
  stopifnot(identical(length(p0), length(lambda), length(sf)),
            identical(length(pi.g), length(mu), length(sigma)))
  
  N = length(p0)
  G = length(mu)
  
  ## simulate Z, indicator for being in foreground.
  Z = matrix(rbinom(N*G, size=1, prob=pi.g), ncol=N)
  
  ## loop over cells to generate counts
  Y = matrix(0, nrow=G, ncol=N)
  for(i in seq_len(N)) {
    ix <- Z[,i] == 0
    Y[ix, i]  = rzip(sum(ix), p0[i], lambda[i])
    Y[!ix, i]  = rLNP(sum(!ix), mu[!ix], sigma[!ix], sf[i])
  }
  return(Y)
}

runDESeq2zinbwave <- function (sceb) 
{
  start <- Sys.time()
  blas_set_num_threads(5)
  # Estimates DE genes
  # Seteo a la variable de los 'grupos' o 'cellTypes' para que esté en formato factor
  colData(sceb)$cellTypes <- factor(x = colData(sceb)$cellTypes)
  # Seteo los rownames de colData y los colnames de assay. Si todo está bien, deberían ser los mismos
  rownames(colData(sceb)) <- 1:length(colData(sceb)[,1])
  colnames(assay(sceb)) <- 1:length(assay(sceb)[1,])
  sceb <-sceb[rowSums(assay(sceb))>0,]
  sceb <- zinbwave::zinbwave(sceb, K=0, observationalWeights=TRUE, epsilon=1e12, BPPARAM = BiocParallel::MulticoreParam(5))
  print("ya pase por zinwave")
  sceb <- scran::convertTo(sceb, type = "DESeq2")
  sceb <- scran::computeSumFactors(sceb)
  design(sceb) <- formula(~ cellTypes)
  sd_vst <- DESeq2::varianceStabilizingTransformation(sceb, fitType = "local", blind = TRUE) %>% 
    assay() %>% apply(MARGIN = 1, sd) %>% sort(decreasing = T)
  #VariableFeatures <- names(which(sd_vst>1))
  VariableFeatures <- names(sort(sd_vst, decreasing = T)[1:2000])
  sceb <- DESeq2::estimateDispersions(object = sceb, fitType = "local")
  sceb <- sceb[VariableFeatures, ]
  sceb <- DESeq2::nbinomLRT(object = sceb, reduced= ~1, minmu = 1e-6)
  sceb <-  DESeq2::replaceOutliers(sceb, minReplicates = Inf)
  res_DE <- DESeq2::results(object = sceb, contrast = c("cellTypes", "celltype1", "celltype2"), alpha = 0.05, independentFiltering=FALSE)
  print("ya pase por results de DESeq2")
  res_DE <- data.frame(res_DE) %>% mutate(geneIndex = str_split_fixed(rownames(res_DE), "g", 2)[,2]) %>% 
    mutate(geneIndex = as.integer(geneIndex))
  res_DE_p_vals <- dplyr::select(res_DE, geneIndex, pval = pvalue, fdr = padj) 
  res_DE_p_vals[is.na(res_DE_p_vals)] <- 1
  # Table to store the log2FC
  table <- dplyr::select(res_DE, geneIndex, cont = log2FoldChange)
  blas_set_num_threads(10)
  end <- Sys.time() 
  print(paste0("Time of runDESeq2zinbwave: ", round(end - start, 0), " ", units(end - start)))
  return(list(table = table, cont = res_DE_p_vals))
}

Est2Phase = function(sce, low.prob=0.99){
  # check the input whether it is an SingleCellExperiment Object or a matrix
  if (is(sce, "SingleCellExperiment")){
    Y = round(assays(sce)[[1]])
  }else{Y = round(sce)}
  ## ## initial estimate of prob(X\in bg)
  Cell0=colMeans(Y==0) # each cell has this percentage 0
  par1= apply(Y,2,function(yy) {
    yyy=yy[yy<=1]
    ### consider (0s, 3, 4)
    if (length(table(yyy)) <= 2){
      # print(which(colSums(Y==yy) == nrow(Y)))
      yyy = yy[yy<=2]
    }
    RobustPoi0(yyy)
  })
  ncell = ncol(sce)
  pi0.hat=Cell0/(par1[1,]+(1-par1[1,])*dpois(0,par1[2,]))
  if (any((pi0.hat > 1))) {warning("Zero proportion is greater than estimation.")}
  pi0.hat <- pmin(pi0.hat, 1)
  prob0=pi0.hat*par1[1,]+ pi0.hat*(1-par1[1,])*dpois(0,par1[2,]) ## ZIP prob at 0
  ## First round
  ## get the 1-low.prob quantile of ZIP
  x0=qpois(pmax(1-(1-low.prob)/(1-par1[1,]),0),par1[2,])
  Z= sweep(Y,2,x0)>0 # indicate if a gene is > bg
  L=colSums(Y*Z)/1e6 # so far it is like simple total..
  
  mu.g1=log2(rowSums(Z*Y)/rowSums(sweep(Z,2,L,FUN="*")))
  mu.g1[is.na(mu.g1)]=0 ## if allZ is 0, it gets NA,
  ### but we should shrink mu.g1 as well since some mu.g1 is estimated by only a few observations
  ## leave it here for now.
  n.g1=rowSums(Z)
  y1=log2(sweep(Y,2,L,FUN="/")+1) #like TPM**
  s.g1=sqrt(rowSums(Z*sweep(y1,1,mu.g1)^2)/(n.g1-1)) ## TPM type of SD
  mu.g2 = shrink.mu(mu.g1,s.g1,n.g1)
  ## get sd.g
  res.g1=log2(sweep(Y,2,L,FUN="/")+1)-mu.g1
  ## mad of those res.g1 that are associated with Z==1
  tmp=array(0,dim=c(dim(res.g1),2))
  tmp[,,1]=res.g1;tmp[,,2]=Z
  sd.g1=apply(tmp,1,function(xx) my.mad(xx[xx[,2]==1,1])) #si quiero que asea una BN, calcular la dispersion
  sd.g1[is.na(sd.g1)]=0## if all bg, there's no info about fg sd
  ## add a shrinkage for sd.g1
  sd.prior=squeezeVar(sd.g1^2,n.g1-1)
  sd.g2=sqrt(sd.prior$var.post)
  #####  gene specific bg. Z_gi
  den.fg = den.bg = NA*Y
  for(i in seq_len(ncell)){
    den.bg[,i]=dZinf.pois(Y[,i], par1[1,i], par1[2,i]) #Qiero que sea ZI?
    den.fg[,i]=dLNP2(x=Y[,i], mu=mu.g1, sigma=sd.g2, l=L[i]) #Cambiar esto por una binomial negativa
  }
  Z.fg=sweep(den.fg,2,1-pi0.hat,FUN="*")
  Z.bg=sweep(den.bg,2,pi0.hat,FUN="*")
  post.Z=Z.fg/(Z.fg+Z.bg)
  post.Z[is.na(post.Z)] <- 1
  
  # ### if I shrink mu.g
  # den.fg2 = NA*Y
  # for (i in seq_len(ncell)){
  #   den.fg2[,i]= dLNP2(x=Y[,i], mu=mu.g2, sigma=sd.g2, l=L[i])
  # }
  # Z.fg2=sweep(den.fg2,2,1-pi0.hat,FUN="*")
  # post.Z2=Z.fg2/(Z.fg2+Z.bg)
  # post.Z2[is.na(post.Z2)] <- 1
  
  pi.g = rowMeans(post.Z)  # pi.g = rowMeans(post.Z2)
  est_rslt = list(exprs = Y, pi.g = pi.g, p0 = par1[1,], lambda = par1[2,],
                  mu = mu.g1, sd = sd.g2, sf = L, bg = Z) #mu = mu.g2
  return(est_rslt)
}

RobustPoi0 = function(x){
  tt=table(x)
  n0=sum(x==0)
  if (n0 == length(x)){
    c("p0"=1,"mu"=0)
  }else{
    if (names(tt)[1]=="0") {
      xx=as.numeric(names(tt))[-1]
      tt=as.vector(tt)[-1]
    }else{
      xx=as.numeric(names(tt))
      tt=as.vector(tt)
    }
    tt=log(tt)+log(gamma(1+xx))
    ##fit the regression without N0
    beta=lm(tt~xx,weight=1/exp(xx))$coef
    mu.hat=exp(beta[2])
    p0=(n0-exp(beta[1]))/(exp(beta[1]+mu.hat)+n0-exp(beta[1]))
    if (any(p0<0)) {warning("Proportion of zero inflation is negative")}
    p0 <- pmax(p0, 0)
    c("p0"=p0,"mu"=mu.hat)
  }
}


## Shrink the mu for sample sizeas are so small
shrink.mu=function(y,s,n){
  mu.g=rep(NA,length(y))
  k=which(n>1)
  if (length(k)<length(n)) {fill=TRUE} else {fill=FALSE}
  s=s[k];y=y[k];n=n[k]
  
  mu0=weighted.mean(y,w=n)
  
  s2.total=sum(s^2*(n-1))+sum(n*(y-mu0)^2)
  s2.total=s2.total/sum(n)
  
  s2.1=sum(s^2*(n-1))/sum(n)
  s2.0=s2.total-s2.1
  ### shrink mu
  mu.sub=  (y*n/s2.1+mu0/s2.0)/(n/s2.1+1/s2.0)
  mu.g[k]=mu.sub
  if (fill) mu.g[-k]=mu0
  
  mu.g
}

my.mad=function (x, center = median(x), constant = 1.4826, na.rm = FALSE){
  if (na.rm)
    x <- x[!is.na(x)]
  res=x-center
  constant * median(res[res>0])
}

dZinf.pois=function(x, p0, mu){
  (x==0)*(p0+(1-p0)*dpois(x,mu))+(x>0)*(1-p0)*dpois(x,mu)
}

dLNP2 <- function(x, mu, sigma, l=1){
  x.min <- pmax(0, x-0.5)
  pnorm(log2((x+0.5)/l), mu, sigma) - pnorm(log2(x.min/l), mu, sigma)
}




### for the use of SC2P
eset2Phase <- function(eset, low.prob=0.99){  ## takes eSet as input
  Y <- round(exprs(eset))
  #################################################
  ## ## initial estimate of prob(X\in bg)
  ##################################################
  Cell0=colMeans(Y==0) # each cell has this percentage 0
  par1=apply(Y,2,function(yy) {
    yy=yy[yy<=15]
    RobustPoi0(yy)}
  )
  pi0.hat=Cell0/(par1[1,]+(1-par1[1,])*dpois(0,par1[2,]))
  if (any((pi0.hat > 1))) {warning("Zero proportion is greater than estimation.")}
  pi0.hat <- pmin(pi0.hat, 1)
  prob0=pi0.hat*par1[1,]+ pi0.hat*(1-par1[1,])*dpois(0,par1[2,]) ## ZIP prob at 0
  ############################################
  ## First round
  ###########################################
  ## get the 1-low.prob quantile of ZIP
  x0=qpois(pmax(1-(1-low.prob)/(1-par1[1,]),0),par1[2,])
  Z= sweep(Y,2,x0)>0 # indicate if a gene is > bg
  L=colSums(Y*Z)/1e6 # so far it is like simple total..
  
  mu.g1=log2(rowSums(Z*Y)/rowSums(sweep(Z,2,L,FUN="*")))
  mu.g1[is.na(mu.g1)]=0 ## if allZ is 0, it gets NA,
  ### but we should shrink mu.g1 as well since some mu.g1 is estimated by only a few observations
  ## leave it here for now.
  n.g1=rowSums(Z)
  y1=log2(sweep(Y,2,L,FUN="/")+1) #like CPM**
  s.g1=sqrt(rowSums(Z*sweep(y1,1,mu.g1)^2)/(n.g1-1)) ## CPM type of SD
  mu.g2 = shrink.mu(mu.g1,s.g1,n.g1)
  ###############################################
  ## get sd.g
  ############################################
  res.g1=log2(sweep(Y,2,L,FUN="/")+1)-mu.g1
  ## mad of those res.g1 that are associated with Z==1
  tmp=array(0,dim=c(dim(res.g1),2))
  tmp[,,1]=res.g1;tmp[,,2]=Z
  sd.g1=apply(tmp,1,function(xx) my.mad(xx[xx[,2]==1,1]))
  sd.g1[is.na(sd.g1)]=0## if all bg, there's no info about fg sd
  ## add a shrinkage for sd.g1
  sd.prior=squeezeVar(sd.g1^2,n.g1-1)
  sd.g2=sqrt(sd.prior$var.post)
  ####################################### ########
  #####  gene specific bg. Z_gi
  #######################
  den.fg = den.bg = NA*Y
  ncell = ncol(Y)
  for(i in seq_len(ncell)){
    den.bg[,i]=dZinf.pois(Y[,i], par1[1,i], par1[2,i])
    den.fg[,i]=dLNP2(x=Y[,i], mu=mu.g1, sigma=sd.g2, l=L[i])
  }
  Z.fg=sweep(den.fg,2,1-pi0.hat,FUN="*")
  Z.bg=sweep(den.bg,2,pi0.hat,FUN="*")
  post.Z=Z.fg/(Z.fg+Z.bg)
  post.Z[is.na(post.Z)] <- 1
  
  ### if I shrink mu.g
  den.fg2 = NA*Y
  for (i in seq_len(ncell)){
    den.fg2[,i]= dLNP2(x=Y[,i], mu=mu.g2, sigma=sd.g2, l=L[i])
  }
  Z.fg2=sweep(den.fg2,2,1-pi0.hat,FUN="*")
  post.Z2=Z.fg2/(Z.fg2+Z.bg)
  post.Z2[is.na(post.Z2)] <- 1
  ##################################################
  ## compute offsets
  ##################################################
  Offset = Y*0
  Ylim=range(log2(1+Y)-mu.g1);Xlim=range(mu.g1)
  
  for(i in seq_len(ncell)){
    tmp.y=log2(1+Y[,i])-mu.g2
    subset= post.Z2[,i] > .99
    tmp.Z2 = post.Z2[, i]
    
    ## lm1 <- loess(tmp.y~mu.g1,
    ##              weights=post.Z2[,i]*mu.g2,subset=subset,degree=1,span=.3)
    if (sum(subset) < 2)
      next
    else
      lm1 <- tryCatch(expr = loess(tmp.y ~ mu.g1, weights = tmp.Z2 * mu.g2,
                                   subset = subset, degree = 1, span = span),
                      error = function(e) loess(tmp.y ~ mu.g1, weights = tmp.Z2 * mu.g2,
                                                subset = subset, degree = 1, span = 0.8))
    
    Offset[subset,i]=lm1$fitted
    ## par(mfrow=c(1,2))
    ## plot(mu.g1, log2(1+Y[,i])-mu.g1, pch=16,cex=.6,ylab="",
    ##      col=rgb(1-post.Z2[,i],0,post.Z2[,i],alpha=rowMeans(post.Z2))
    ##     ,ylim=Ylim,xlim=Xlim,main=i)
    ## points(lm1$x,lm1$fitted,col=5)
    
    ## plot(mu.g1, log2(1+Y[,i])-Offset[,i]-mu.g1, pch=16,cex=.6,ylab="",
    ##      col=rgb(1-post.Z2[,i],0,post.Z2[,i],alpha=rowMeans(post.Z2))
    ##     ,ylim=Ylim,xlim=Xlim,main=i)
    ## tmp.y2 <- tmp.y-Offset[,i]
    ## lm2 <- loess(tmp.y2 ~ mu.g1,
    ##              weights=post.Z2[,i]*mu.g2, subset=subset,degree=1,span=.3)
    ## points(lm2$x,lm2$fitted,col=5)
  }
  ##################################################
  ## assemble the estimators into sc2pSet object
  ##################################################
  ## add mu and sd to feature data
  fdata <- fData(eset)
  fdata2 <- as.data.frame(cbind(fdata, mu.g2, sd.g2))
  colnames(fdata2) <- c(colnames(fdata), "mean", "sd")
  fvar <- rbind(fvarMetadata(eset), "mean"="shrinkage estimated foreground mean",
                "sd"="shrinkage estimated foreground standard deviation")
  featureData <- new("AnnotatedDataFrame", data=fdata2,
                     varMetadata=fvar)
  ## add lambda and p0 to phenoData
  pdata <- pData(eset)
  pdata2 <- as.data.frame(cbind(pdata, par1[1,], par1[2,], L))
  colnames(pdata2) <- c(colnames(pdata), "p0", "lambda", "L")
  pvar <-rbind(varMetadata(eset), "p0"="proportion of zero inflation",
               "lambda"="mean of background poisson",
               "L"="foreground library size")
  phenoData <- new("AnnotatedDataFrame", data=pdata2, varMetadata=pvar)
  
  out <- list(exprs=Y, Z=post.Z2, Offset=Offset,
              phenoData=phenoData,
              featureData=featureData,
              experimentData=experimentData(eset),
              annotation=annotation(eset))
  out
}

