#========================================================
# VENN DKOvsWTG vs KOPvsWTP (3 REGIONES)
#========================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(ggvenn)

# 1) CARGAR EXCEL ---------------------------------------
ruta_venn2 <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis Jose 8/DKOvsKOP/Diagrama de Venn Final DKOvsKOP.xlsx"
# Ajusta el nombre si es otro

venn2_excel <- read_excel(ruta_venn2, sheet = 1)
colnames(venn2_excel) <- trimws(colnames(venn2_excel))
print(colnames(venn2_excel))
# Esperado:
# "DKOvsWTG"  "DKOvsWTG and KOPvsWTP"  "KOPvsWTP"

# 2) CONSTRUIR LOS TRES CONJUNTOS -----------------------

genes_DKO2 <- c(
  venn2_excel$`DKOvsWTG`,                   # solo DKOvsWTG
  venn2_excel$`DKOvsWTG and KOPvsWTP`      # comunes DKOvsWTG + KOPvsWTP
) |>
  na.omit() |>
  as.character()

genes_KOP2 <- c(
  venn2_excel$`KOPvsWTP`,                  # solo KOPvsWTP
  venn2_excel$`DKOvsWTG and KOPvsWTP`      # comunes DKOvsWTG + KOPvsWTP
) |>
  na.omit() |>
  as.character()

# Si solo quieres 2 conjuntos (DKO vs KOP), el tercero sería vacío,
# pero como pides 3 columnas, defino un tercero por si acaso:
# Por ejemplo, genes solo comunes (zona intermedia) como conjunto aparte:
genes_COMUN2 <- venn2_excel$`DKOvsWTG and KOPvsWTP` |>
  na.omit() |>
  as.character()

# Lista para venn
venn_list2 <- list(
  "DKOvsWTG"     = unique(genes_DKO2),
  "KOPvsWTP"     = unique(genes_KOP2)
)

# 3) DIAGRAMA DE VENN CON ggvenn ------------------------[web:120]

p_venn2 <- ggvenn(
  venn_list2,
  show_elements   = FALSE,                 # solo números
  show_percentage = TRUE,                  # n + %
  digits          = 1,                     # 1 decimal en %
  fill_color      = c(
    "DKOvsWTG"      = "#ff7f00",           # DKO naranja
    "KOPvsWTP"      = "#33a02c",           # KOP verde
    "DKO+KOP comun" = "#b2df8a"           # verde claro para el tercero
  ),
  stroke_color    = "black",
  stroke_size     = 1,
  set_name_size   = 3.8,
  text_size       = 3.8
) +
  ggtitle("Diagrama de Venn WTG vs WTP vs DKO vs KOP") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    text       = element_text(face = "bold")   # todo en negrita
  )

print(p_venn2)

# 4) GUARDAR --------------------------------------------

salida_dir <- "D:/Carpetas/Tesis/Figuras_Heatmap"
if (!dir.exists(salida_dir)) dir.create(salida_dir, recursive = TRUE)

ggsave(file.path(salida_dir, "Venn_DKOvsWTG_KOPvsWTP_ggvenn_2.png"),
       p_venn2, width = 8, height = 6, dpi = 400, bg = "white")
ggsave(file.path(salida_dir, "Venn_DKOvsWTG_KOPvsWTP_ggvenn_2.pdf"),
       p_venn2, width = 8, height = 6, bg = "white")

cat("✅ Venn DKOvsWTG / KOPvsWTP guardado en:", salida_dir, "\n")
