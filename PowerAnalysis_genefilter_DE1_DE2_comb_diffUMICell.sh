#!/bin/bash

# Este script tiene el proposito de correrel script PowerAnalysis_genefilter.R reiteradas veces para que pueda generar replicas tecnicas del putput del script y promediar los resultados.
# La repeticion no la hago en R directamente porque creo que me colapsa la memoria (ver error que aparece mas abajo) ya que por mas que remueva el objeto pow_rslt_t_DESeq, no libera el espacio correspondiente. Correrlo mediante un script de bash, obliga a cerrar la sesion de R, eliminar los procesos que R esta corriendo y restaurar el espacio en memoria para la repeticion
# Library size (l) equivalences:
# 5 = 1600, 1000 = 3000, 2000 = 5000 (in mean(UMIs/Cell))

# Processing
for i in  1 2 3 4
do
	j="c(250,500)"
	mu=2
	sd=0.5
	k=0.3
	for l in 5 1000 2000
	do
		"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter_dpig_mu_comb_differentUMIcell.R $i $j $mu $sd $k $l;
		pkill R
	done
	j="c(1000)"
	for l in 5 1000 2000
	do
		"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter_dpig_mu_comb_differentUMIcell.R $i $j $mu $sd $k $l;
		pkill R
	done
	j="c(2000)"
	for l in 5 1000 2000
	do
		"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter_dpig_mu_comb_differentUMIcell.R $i $j $mu $sd $k $l;
		pkill R
	done
done


