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
samples <- read.csv('meta/samplesheet_chr12.csv')
View(samples)
```

