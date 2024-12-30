---
title: "ChIP-seq workflow using R"
author: "Adil Hannaoui Anaaoui"
date: "December 12th, 2024"
---

Antes de proceder con el análisis diferencial de los datos, R nos permite, a través del paquete ChIPQC, evaluar la calidad de los archivos ChIP-seq antes de trabajar con ellos.

```
## Load libraries
library(ChIPQC)

## Load sample data
samples <- read.csv('meta/metadata.csv')
View(samples)
```

Para poder cargar correctamente y analizar las muestras con las que vamos a trabajar, vamos a necesitar un archivo en formato csv, con los metadatos de las muestras, reflejando el tipo de muestra de cada una, su tratamiento, el tejido, y demás información necesaria. Esto es importante ya que permite realizar distintos controles de calidad para conocer la fortaleza de nuestros datos. Por ejemplo, para realizar correctamente un PCA, es necesario clasificar correctamente las muestras entre IP e input.

Entre las distintas características que se deben recoger en el csv, tenemos los siguientes parámetros:

* **SampleID**: Identifier string for sample
* **Tissue, Factor, Condition**: Identifier strings for up to three different factors (You will need to have all columns listed. If you don't have infomation, then set values to NA)
* **Replicate**: Replicate number of sample
* **bamReads**: file path for BAM file containing aligned reads for ChIP sample
* **ControlID**: an identifier string for the control sample
* **bamControl**: file path for bam file containing aligned reads for control sample
* **Peaks**: path for file containing peaks for sample
* **PeakCaller**: Identifier string for peak caller used. Possible values include “raw”, “bed”, “narrow”, “macs”

```
## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="hg19") 
```

Una vez creado el archivo de metadatos, podemos proceder con la creación del objeto. Es importante recalcar que ChIPQC, hasta la fecha, solo tiene compatibilidad con humanos y ratones. Si estuviéramos trabajando con cerevisiae, no se podría realizar el control de calidad de esta manera.

```
## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: input and IP", reportFolder="ChIPQCreport")
```

El siguiente paso consiste en crear el report de ChIPQC, para cuando queramos ejecutar el objeto

El resultado del control de calidad tendrá una estructura similar a la siguiente:

<img src="../img/QCsummary.png">

- Read characteristics
- Enrichment of reads in peaks
- Peak signal strength
- Peak profiles

