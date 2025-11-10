# Instalar librerías si aún no están instaladas
install.packages(c("readxl", "pheatmap", "dplyr", "tibble"))

# Cargar librerías
library(readxl)
library(pheatmap)
library(dplyr)
library(tibble)

# Ruta exacta al archivo

# Carga el archivo
genes_heatmap <- read_excel("D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Nuevo analisis Jose/GCKIIIvsWTG/Genes_heatmap_GCKIIvsWTG.xlsx")
View(genes_heatmap)

# Leer el archivo Excel (asumiendo que los datos están en la primera hoja)

print(colnames(genes_heatmap))

# Convertir la primera columna (genes) en nombres de fila
# Asegúrate de que la primera columna realmente contiene los nombres de los genes, ademas quito los valores duplicados para que no de error
genes_heatmap$gene_name <- make.unique(as.character(genes_heatmap$gene_name))

datos <- genes_heatmap %>% column_to_rownames(var = "gene_name")
print(colnames(datos))
datos[] <- lapply(datos, function(x) as.numeric(as.character(x)))

# Transformar los valores de expresión a log2(FPKM + 1)
datos_log <- log2(datos + 1)

# Verifica que solo quedan columnas numéricas (fpkm)
print(str(datos))

# Limpiar nombres de columnas: quitar "_fpkm"
colnames(datos) <- gsub("_fpkm$", "", colnames(datos))

# Crear anotación: extraer grupo y subgrupo
grupos <- sapply(strsplit(colnames(datos), "_"), `[`, 1)

anotacion <- data.frame(
  Grupo = factor(grupos))

rownames(anotacion) <- colnames(datos)

# Log2(FPKM + 1)
datos_log <- log2(datos + 1)


# 9. Definir colores para los grupos
grupo_colores <- c(
  "WTG" = "#6a3d9a",   # morado
  "KOS" = "#1b9e77",   # verde azulado
  "KOM" = "#1f78b4",   # azul claro
  "DKO" = "#e7298a"    # rosa fuerte
)

ann_colors <- list(Grupo = grupo_colores)


salida_png <- file.path(
  "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Nuevo analisis Jose/GCKIIIvsWTG",
  "Heatmap_final2_GCKIIIvsWTG.png"
)

salida_pdf <- file.path(
  "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Nuevo analisis Jose/GCKIIIvsWTG",
  "Heatmap_final2_GCKIIIvsWTG.pdf"
)

# Verificar que la ruta se ve correcta
print(salida)

# Heatmap estándar publicable
pheatmap(
  datos_log,
  annotation_col = anotacion,
  annotation_colors = ann_colors,
  scale = "row",                        # escala por gen
  clustering_distance_rows = "correlation",  # Pearson
  clustering_distance_cols = "correlation",  # Pearson
  clustering_method = "average",        # método estándar
  show_rownames = FALSE,                # genes ocultos por claridad
  fontsize_col = 10,
  main = "Heatmap GCKIIIvsWTG",
  filename = salida_pdf, width = 10, height = 10
)

pheatmap(
  datos_log,
  annotation_col = anotacion,
  annotation_colors = ann_colors,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  show_rownames = FALSE,
  fontsize_col = 10,
  main = "Heatmap GCKIIIvsWTG",
  filename = salida_png, width = 10, height = 10
)
