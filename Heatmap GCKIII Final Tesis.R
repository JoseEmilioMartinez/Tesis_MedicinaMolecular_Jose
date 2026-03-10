# Cargar librerías
library(readxl)
library(pheatmap)
library(dplyr)
library(tibble)
library(stringr)

# Cargar datos
genes_heatmap <- read_excel("D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis Jose 8/GCKIII/Genes Heatmap GCKIII Final.xlsx")
print(colnames(genes_heatmap))

# Ruta específica que pediste
carpeta_salida <- "D:\\Carpetas\\Experimentos\\Celulas Endoteliales\\RNA\\RNAseq\\Resultados\\Graficos Final Tesis"

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

# Guardar grupos ORIGINALES antes recode
grupos_originales <- sapply(strsplit(colnames(datos), "_"), `[`, 1)

# RECODE PARA COLORES
grupos_recodificados <- recode(grupos_originales,
                               "WTG" = "STK24 WT STK25 WT",
                               "DKO" = "STK24 KO STK25ECKO",
                               "KOM" = "STK24 KO",
                               "KOS" = "STK25ECKO")
                               

anotacion <- data.frame(Grupo = factor(grupos_recodificados))
rownames(anotacion) <- colnames(datos)
print(table(anotacion$Grupo))

  
# grupo_colores

grupo_colores <- c(
  "STK24 WT STK25 WT" = "#1f78b4",
  "STK24 KO STK25ECKO" = "#ff7f00",
  "STK24 KO" = "#e7298a",
  "STK25ECKO" = "#6a3d9a" 
)

ann_colors <- list(Grupo = grupo_colores)

# Rutas finales de salida
salida_png <- file.path(carpeta_salida, "Heatmap_FinalFiltrado_GCKIII_3.png")
salida_pdf <- file.path(carpeta_salida, "Heatmap_FinalFiltrado_GCKIII_3.pdf")

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
    main = "Heatmap STK24/25 vs STK24 vs STK25",
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

