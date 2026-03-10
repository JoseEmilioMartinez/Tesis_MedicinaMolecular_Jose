#========================================================
# DIAGRAMA DE VENN GCKIII (DKOvsWTG, KOMvsWTG, KOSvsWTG)
# USANDO ggvenn (más control que ggVennDiagram)
#========================================================

# 1) CARGAR LIBRERÍAS
# Ejecuta esta línea una sola vez si no tienes ggvenn instalado:
# install.packages("ggvenn")
install.packages("ggvenn")
library(readxl)
library(dplyr)
library(ggplot2)
library(ggvenn)

# 2) CARGAR EXCEL CON LAS LISTAS
ruta_venn <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis Jose 8/GCKIII/Diagrama de Venn Final GCKIII.xlsx"

venn_excel <- read_excel(ruta_venn, sheet = 1)
colnames(venn_excel) <- trimws(colnames(venn_excel))
print(colnames(venn_excel))

# Espero columnas:
# "exclusively in DKOvsWTG"
# "common DKOvsWTG and KOMvsWTG"
# "exclusively in KOMvsWTG"
# "KOMvsWTG and KOSvsWTG"
# "KOSvsWTG"
# "DKOvsWTG and KOSvsWTG"
# "DKOvsWTG, KOMvsWTG and KOSvsWTG"

# 3) CONSTRUIR LOS TRES CONJUNTOS

genes_DKO <- c(
  venn_excel$`exclusively in DKOvsWTG`,
  venn_excel$`common DKOvsWTG and KOMvsWTG`,
  venn_excel$`DKOvsWTG and KOSvsWTG`,
  venn_excel$`DKOvsWTG, KOMvsWTG and KOSvsWTG`
) |> na.omit() |> as.character()

genes_KOM <- c(
  venn_excel$`exclusively in KOMvsWTG`,
  venn_excel$`common DKOvsWTG and KOMvsWTG`,
  venn_excel$`KOMvsWTG and KOSvsWTG`,
  venn_excel$`DKOvsWTG, KOMvsWTG and KOSvsWTG`
) |> na.omit() |> as.character()

genes_KOS <- c(
  venn_excel$`KOSvsWTG`,
  venn_excel$`KOMvsWTG and KOSvsWTG`,
  venn_excel$`DKOvsWTG and KOSvsWTG`,
  venn_excel$`DKOvsWTG, KOMvsWTG and KOSvsWTG`
) |> na.omit() |> as.character()

# Lista para venn (nombres tal y como quieres que salgan)
venn_list <- list(
  "DKOvsWTG" = unique(genes_DKO),
  "KOMvsWTG" = unique(genes_KOM),
  "KOSvsWTG" = unique(genes_KOS)
)

# 4) DIAGRAMA DE VENN CON ggvenn (MUY CONFIGURABLE) [web:120][web:124]

p_venn <- ggvenn(
  venn_list,
  show_elements   = FALSE,                 # no mostrar nombres de genes
  show_percentage = TRUE,                  # mostrar también %
  digits          = 1,                     # 1 decimal en el porcentaje
  fill_color      = c(                    # color de cada círculo
    "DKOvsWTG" = "#ff7f00",  # DKO
    "KOMvsWTG" = "#e7298a",  # KOM
    "KOSvsWTG" = "#6a3d9a"   # KOS
  ),  # colores de cada círculo
  stroke_color    = "black",               # bordes negros
  stroke_size     = 1,                   # grosor del borde
  set_name_size   = 3.5,                     # tamaño de DKOvsWTG, KOMvsWTG, KOSvsWTG
  text_size       = 3.5                    # tamaño de los números n / %
) +
  ggtitle("Diagrama de Venn GCKIII") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    text       = element_text(face = "bold")
  )

print(p_venn)

# 5) GUARDAR EN LA CARPETA DE TESIS

salida_dir <- "D:/Carpetas/Tesis/Figuras_Heatmap"
if (!dir.exists(salida_dir)) dir.create(salida_dir, recursive = TRUE)

ggsave(file.path(salida_dir, "Venn_GCKIII_DKO_KOM_KOS_ggvenn_3.png"),
       p_venn, width = 8, height = 6, dpi = 400, bg = "white")
ggsave(file.path(salida_dir, "Venn_GCKIII_DKO_KOM_KOS_ggvenn_3.pdf"),
       p_venn, width = 8, height = 6, bg = "white")

cat("✅ Venn ggvenn guardado en:", salida_dir, "\n")
