# Cargar librerías
library(readxl)
library(pheatmap)
library(dplyr)
library(tibble)

# Cargar datos
genes_heatmap <- read_excel("D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis Jose 8/GCKIII/Genes Heatmap GCKIII Final.xlsx")
print(colnames(genes_heatmap))

# Salida: carpeta nueva dentro de TESIS (CREA SI NO EXISTE)
carpeta_tesis <- "D:/Carpetas/Tesis"
nombre_carpeta_salida <- "Figuras_Heatmap"

carpeta_salida <- file.path(carpeta_tesis, nombre_carpeta_salida)

# Crear carpeta si no existe ✅ CORREGIDO
if (!dir.exists(carpeta_salida)) {
  dir.create(carpeta_salida, recursive = TRUE)
  cat("✅ Carpeta creada:", carpeta_salida, "\n")
} else {
  cat("✅ Carpeta ya existe:", carpeta_salida, "\n")
}

# Preparar datos
genes_heatmap$gene_name <- make.unique(as.character(genes_heatmap$gene_name))
datos <- genes_heatmap %>% column_to_rownames(var = "gene_name")

# Asegurar numérico
datos[] <- lapply(datos, function(x) as.numeric(as.character(x)))
str(datos)

# Limpiar nombres de columnas
colnames(datos) <- gsub("_fpkm$", "", colnames(datos))

# Log2(FPKM + 1)
datos_log <- log2(datos + 1)

# Anotación por grupo ✅ AÑADIDO KOP
grupos <- sapply(strsplit(colnames(datos), "_"), `[`, 1)
anotacion <- data.frame(Grupo = factor(grupos))
rownames(anotacion) <- colnames(datos)

grupo_colores <- c(
  "WTG" = "#1f78b4",
  "DKO" = "#ff7f00",
  "KOM" = "#e7298a",
  "KOS" = "#6a3d9a"  # 
)
ann_colors <- list(Grupo = grupo_colores)

# Rutas finales de salida
salida_png <- file.path(carpeta_salida, "Heatmap_FinalFiltrado_GCKIII_2.png")
salida_pdf <- file.path(carpeta_salida, "Heatmap_FinalFiltrado_GCKIII_2.pdf")

cat("📁 Guardará:\n", salida_png, "\n", salida_pdf, "\n")

# Heatmap PDF ✅ TÍTULO CORREGIDO

crear_heatmap <- function(filename, width=12, height=11) {
  pheatmap(
    datos_log,
    annotation_col = anotacion,
    annotation_colors = ann_colors,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "average",
    
    # ETIQUETAS SUBGRUPOS: NEGRITA Y GRANDES
    fontsize_col = 14,
    
    # LEYENDA GRUPOS: NEGRITA Y GRANDE
    fontsize = 12,
    annotation_legend_type = "heatmap",
    
    # TÍTULO SEPARADO Y GRANDE
    main = "Heatmap GCKIII",
    fontsize_main = 18,
    
    # MARGENES Y LAYOUT
    show_rownames = FALSE,
    gaps_row = 10,
    treeheight_row = 50,
    treeheight_col = 20,
    
    filename = filename, 
    width = width, 
    height = height
  )
}

# Generar PDF (tesis/publicación)
crear_heatmap(salida_pdf, width=12, height=11)

# Generar PNG (presentaciones)
crear_heatmap(salida_png, width=12, height=11)

cat("✅ ¡Heatmaps guardados en la carpeta de tesis!\n")
