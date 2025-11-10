# Cargar librerías
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidyr)
library(zoo)
library(pheatmap)
library(RColorBrewer)

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

# --- Seleccionar top vías significativas ---
df_top <- df %>%
  filter(FDR_qval < 0.05) %>%
  group_by(Comparacion) %>%
  slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
  ungroup()
df_comp$GS_DETAILS <- str_wrap(df_comp$GS_DETAILS, width = 35)

# --- Gráfico combinado ---
for (comp in unique(df_top$Comparacion)) {
  
  df_comp <- df_top %>% filter(Comparacion == comp)
  
  g_comb <- ggplot(df_comp, aes(x = NES, y = fct_reorder(GS_DETAILS, NES), fill = NES)) +
    geom_col(width = 0.6) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "NES"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title = paste0(comp, " — RUTAS AUMENTADAS Y DISMINUIDAS"),
      x = "NES",
      y = "GENE SET"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 10)
    )
  
 # Crear carpeta de salida
  output_dir <- "D:/Carpetas/Tesis/Graficos_GSEA"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Guardar en PNG y PDF
  ggsave(
    filename = paste0(output_dir, "/", comp, "_PanelCombinado.png"),
    plot = g_comb, width = 8, height = 5, dpi = 300
  )
  ggsave(
    filename = paste0(output_dir, "/", comp, "_PanelCombinado.pdf"),
    plot = g_comb, width = 8, height = 5
  )
}
