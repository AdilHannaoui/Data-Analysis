# Análisis de Secuencias Biológicas: RNA-seq, ChIP-seq y RIP-seq

Este repositorio está diseñado para proporcionar flujos de trabajo y herramientas relacionados con el análisis de datos biológicos de secuencias, tanto de ARN como de ADN. Se enfoca en tres técnicas ampliamente utilizadas en biología molecular y bioinformática:

- **RNA-seq**: Para el análisis de expresión génica y cuantificación de ARN.  
- **ChIP-seq**: Para estudiar interacciones proteína-ADN y regiones de unión en el genoma.  
- **RIP-seq**: Para identificar interacciones proteína-ARN y explorar regiones de unión específicas.  

Las similitudes metodológicas entre **ChIP-seq** y **RIP-seq** se aprovechan para optimizar los flujos de trabajo, mientras que las particularidades de cada técnica se tratan de manera específica para garantizar resultados precisos.

---

## Contenido del Repositorio

- **Pipelines Automatizados**  
  Scripts y herramientas para llevar a cabo análisis de principio a fin, desde la calidad de las secuencias hasta la obtención de resultados interpretables.

- **Análisis Diferencial**  
  Métodos para identificar diferencias significativas en los datos, incluyendo picos diferenciales (ChIP-seq, RIP-seq) y genes diferencialmente expresados (RNA-seq).  

- **Visualización**  
  Herramientas para la creación de gráficos y visualizaciones, como heatmaps, gráficos de enriquecimiento GO, y anotaciones genómicas.

- **Documentación Detallada**  
  Instrucciones paso a paso y ejemplos prácticos para la ejecución de los análisis.

---

## Tecnologías y Herramientas

### RNA-seq
- HISAT2 para el alineamiento y cuantificación.
- featureCounts para la asignación de lecturas.
- DESeq2, Limma-Voom o edgeR para análisis diferencial.

### ChIP-seq y RIP-seq
- MACS2 para la identificación de picos.  
- Bedtools, Samtools y deepTools para el manejo y visualización de datos.  
- DiffBind y clusterProfiler para análisis diferencial y enriquecimiento funcional.  
