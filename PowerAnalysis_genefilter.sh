#!/bin/bash

# Este script tiene el proposito de correrel script PowerAnalysis_genefilter.R reiteradas veces para que pueda generar replicas tecnicas del putput del script y promediar los resultados.
# La repeticion no la hago en R directamente porque creo que me colapsa la memoria (ver error que aparece mas abajo) ya que por mas que remueva el objeto pow_rslt_t_DESeq, no libera el espacio correspondiente. Correrlo mediante un script de bash, obliga a cerrar la sesion de R, eliminar los procesos que R esta corriendo y restaurar el espacio en memoria para la repeticion

# Processing
for i in  1 2 3 4
do
	j="c(250,500)"
	mu=(1) # 2 4 8)
	sd=(0.5) # 1 2 4)
	for k in 0.1 0.2 0.3 0.4
	do
		for m in 0 #1 2 3
		do
			"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter.R $i $j ${mu[m]} ${sd[m]} $k;
			pkill R
		done
	done
	j="c(1000)"
	for k in 0.1 0.2 0.3 0.4
	do
		for m in 0 #1 2 3
		do
			"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter.R $i $j ${mu[m]} ${sd[m]} $k;
			pkill R
		done
	done
	j="c(2000)"
	for k in 0.1 0.2 0.3 0.4
	do
		for m in 0 #1 2 3
		do
			"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter.R $i $j ${mu[m]} ${sd[m]} $k;
			pkill R
		done
	done
#	j="c(1000)"
#	for k in 0 1 2 3
#	do
#		"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter.R $i $j ${mu[k]} ${sd[k]};
#		pkill R
#	done
#	j="c(2000)"
#	for k in 0 1 2 3
#	do
#		"/opt/R/4.1.1/bin"/Rscript Scripts/PowerAnalysis_genefilter.R $i $j ${mu[k]} ${sd[k]};
#		pkill R
#	done
done


