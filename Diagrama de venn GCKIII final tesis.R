#========================================================
# VENN GCKIII con NOMBRES COMPLETOS STK24/STK25 WT
#========================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(ggvenn)

ruta_venn <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis Jose 8/GCKIII/Diagrama de Venn Final GCKIII.xlsx"  # Ajusta ruta

venn_excel <- read_excel(ruta_venn, sheet = 1)
colnames(venn_excel) <- trimws(colnames(venn_excel))
print("Columnas:")
print(colnames(venn_excel))

# MAPEO NOMBRES COMPLETOS (exactos de tu GSEA Excel)
nombres_completos <- list(
  "STK24 KO STK25ECKO\nvs STK24 WT STK25 WT" = c("exclusively in DKOvsWTG", 
                                                 "common DKOvsWTG and KOMvsWTG", 
                                                 "DKOvsWTG and KOSvsWTG", 
                                                 "DKOvsWTG, KOMvsWTG and KOSvsWTG"),
  
  "STK24 KO\nvs STK24 WT STK25 WT" = c("exclusively in KOMvsWTG", 
                                       "common DKOvsWTG and KOMvsWTG", 
                                       "KOMvsWTG and KOSvsWTG", 
                                       "DKOvsWTG, KOMvsWTG and KOSvsWTG"),
  
  "STK25ECKO\nvs STK24 WT STK25 WT" = c("KOSvsWTG", 
                                        "KOMvsWTG and KOSvsWTG", 
                                        "DKOvsWTG and KOSvsWTG", 
                                        "DKOvsWTG, KOMvsWTG and KOSvsWTG")
)

# Construir listas
venn_list <- list()
for (nombre_largo in names(nombres_completos)) {
  cols <- nombres_completos[[nombre_largo]]
  genes <- unlist(lapply(cols[cols %in% colnames(venn_excel)], 
                         function(col) venn_excel[[col]]))
  venn_list[[nombre_largo]] <- unique(na.omit(as.character(genes)))
}

print("Tamaños:")
print(sapply(venn_list, length))

# COLORES COHERENTES CON GSEA (DKO=rojo/naranja, KOM=azul, KOS=verde?)
p_venn <- ggvenn(
  venn_list,
  show_elements = FALSE,
  show_percentage = TRUE,
  digits = 1,
  fill_color = c(
    "STK24 KO STK25ECKO\nvs STK24 WT STK25 WT" = "#ff7f00",  
    "STK24 KO\nvs STK24 WT STK25 WT" = "#e7298a",           
    "STK25ECKO\nvs STK24 WT STK25 WT" = "#6a3d9a"           
  ),
  stroke_color = "black",
  stroke_size = 1.5,
  set_name_size = 4,      # Tamaño leyenda nombres largos
  text_size = 4
) +
  labs(title = "Diagrama de Venn STK24/25 vs STK24 vs STK25") +
    theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
    legend.text = element_text(size = 10, face = "bold"),
    plot.margin = margin(25, 25, 25, 25)
  )

print(p_venn)

# Guardar
salida_dir <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Graficos Final Tesis"
dir.create(salida_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(salida_dir, "Venn_GCKIII_STK24_STK25_completo_2.png"), p_venn, 
       width = 10, height = 8, dpi = 400, bg = "white")
ggsave(file.path(salida_dir, "Venn_GCKIII_STK24_STK25_completo_2.pdf"), p_venn, 
       width = 10, height = 8, bg = "white")

cat("✅ Venn con nombres STK24 WT STK25 WT guardado!\n")
