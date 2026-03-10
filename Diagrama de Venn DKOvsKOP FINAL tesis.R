#========================================================
# VENN GCKIII con NOMBRES COMPLETOS STK24/STK25 WT
#========================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(ggvenn)

ruta_venn2 <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis Jose 8/DKOvsKOP/Diagrama de Venn Final DKOvsKOP.xlsx"

venn2_excel <- read_excel(ruta_venn2, sheet = 1)
colnames(venn2_excel) <- trimws(colnames(venn2_excel))
print("Columnas detectadas:")
print(colnames(venn2_excel))

# NOMBRES COMPLETOS EXACTOS (de tu GSEA)
nombres_completos2 <- list(
  "STK24 KO STK25ECKO\nvs STK24 WT STK25 WT" = c("DKOvsWTG", "DKOvsWTG and KOPvsWTP"),
  "CCM3ECKO\nvs CCM3 WT" = c("KOPvsWTP", "DKOvsWTG and KOPvsWTP")
)

# Construir listas automáticamente
venn_list2 <- list()
for (nombre_largo in names(nombres_completos2)) {
  cols <- nombres_completos2[[nombre_largo]]
  cols_exist <- cols[cols %in% colnames(venn2_excel)]
  if (length(cols_exist) > 0) {
    genes <- unlist(lapply(cols_exist, function(col) venn2_excel[[col]]))
    venn_list2[[nombre_largo]] <- unique(na.omit(as.character(genes)))
  }
}

print("Tamaños conjuntos:")
print(sapply(venn_list2, length))

# Venn 2 conjuntos (DKO rojo, CCM3 azul - coherente con GSEA)
p_venn2 <- ggvenn(
  venn_list2,
  show_elements = FALSE,
  show_percentage = TRUE,
  digits = 1,
  fill_color = c(
    "STK24 KO STK25ECKO\nvs STK24 WT STK25 WT" = "#ff7f00",  # Rojo DKO
    "CCM3ECKO\nvs CCM3 WT" = "#33a02c"                      # Azul CCM3ECKO
  ),
  stroke_color = "black",
  stroke_size = 1.5,
  set_name_size = 4,      # Tamaño leyenda nombres largos
  text_size = 4
) +
  labs(title = "Diagrama de Venn CCM3 vs STK24/25") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
    legend.text = element_text(size = 10, face = "bold"),
    plot.margin = margin(25, 25, 25, 25)
  )

print(p_venn2)

# Guardar en carpeta unificada GSEA/Venn
salida_dir <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Graficos Final Tesis"
dir.create(salida_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(salida_dir, "Venn_DKOvsWTG_CCM3ECKOvsWT_2.png"), p_venn2, 
       width = 9, height = 7, dpi = 400, bg = "white")
ggsave(file.path(salida_dir, "Venn_DKOvsWTG_CCM3ECKOvsWT_2.pdf"), p_venn2, 
       width = 9, height = 7, bg = "white")

cat("✅ Venn DKO vs CCM3ECKO guardado en:", salida_dir, "\n")
cat("Genes comunes DKO∩CCM3ECKO:", length(Reduce(intersect, venn_list2)), "\n")
