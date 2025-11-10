# Cargar librerías
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidyr)
library(zoo)

# Leer el archivo Excel
ruta <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis nuevo datos filtrados/GSEA filtrado heatmap/Genesets GSEA.xlsx"

# Leer la primera hoja con nombres de columna
df <- read_excel(ruta, sheet = 1, col_names = TRUE)

# 🔹 LIMPIEZA DE DATOS
# Quitar filas completamente vacías
df <- df[rowSums(is.na(df) | df == "") < ncol(df), ]

# Dar nombres de columnas esperados (ajusta si es necesario)
colnames(df) <- c("Comparacion", "GS_DETAILS", "Size", "ES", "NES", "NOM_pval", "FDR_qval")

# Rellenar la columna "Comparacion" hacia abajo
df$Comparacion <- zoo::na.locf(df$Comparacion, na.rm = FALSE)

# Quitar filas sin datos relevantes
df <- df %>% filter(!is.na(GS_DETAILS) & GS_DETAILS != "")

# Asegurar que las columnas numéricas sean numéricas
df <- df %>%
  mutate(
    Size = as.numeric(Size),
    ES = as.numeric(ES),
    NES = as.numeric(NES),
    NOM_pval = as.numeric(NOM_pval),
    FDR_qval = as.numeric(FDR_qval)
  )
# --- 🔹 FILTRADO ---
# Quedarse solo con FDR < 0.05
df_filt <- df %>% filter(FDR_qval < 0.05)

# Seleccionar las 5 vías más significativas por comparación (mayor |NES|)
df_top <- df_filt %>%
  group_by(Comparacion) %>%
  slice_max(order_by = abs(NES), n = 5, with_ties = FALSE) %>%
  ungroup()

# Ordenar gene sets para mostrar los más extremos arriba
df_top <- df_top %>%
  mutate(GS_DETAILS = fct_reorder(GS_DETAILS, NES))

# --- 🎨 DOT PLOT ---
g_dot <- ggplot(df_top, aes(x = Comparacion, 
                            y = GS_DETAILS,
                            size = Size, 
                            color = NES)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size(range = c(2, 8), name = "Tamaño") +
  labs(
    title = "GSEA — Principales vías por comparación (FDR < 0.05)",
    x = "Comparación",
    y = "GENE SET"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# --- 💾 EXPORTAR ---
output_dir <- "D:/Carpetas/Tesis/Graficos_GSEA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(paste0(output_dir, "/DotPlot_GSEA.png"), plot = g_dot, width = 9, height = 6, dpi = 300)
ggsave(paste0(output_dir, "/DotPlot_GSEA.pdf"), plot = g_dot, width = 9, height = 6)
