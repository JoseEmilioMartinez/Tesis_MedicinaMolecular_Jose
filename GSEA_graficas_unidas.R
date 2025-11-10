#  BARPLOT GSEA: Top 5 vías aumentadas y disminuidas por comparación
# =======================================================

# Cargar librerías
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidyr)
library(zoo)
library(gridExtra) 
library(patchwork)


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

# Crear columna indicando dirección
df <- df %>%
  mutate(Direccion = ifelse(NES > 0, "Aumentado", "Disminuido"))

# Crear carpeta de salida
output_dir <- "D:/Carpetas/Tesis/Graficos_GSEA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Crear listas para almacenar los gráficos
plots_aumentadas <- list()
plots_disminuidas <- list()

# --- Bucle por comparación ---
for (comp in unique(df$Comparacion)) {
  
  df_comp <- df %>% filter(Comparacion == comp)
  
  # Top 5 aumentadas y top 5 disminuidas
  df_pos <- df_comp %>% filter(NES > 0) %>% slice_max(order_by = NES, n = 5)
  df_neg <- df_comp %>% filter(NES < 0) %>% slice_max(order_by = abs(NES), n = 5)
  
  # Ajuste de nombres largos (máximo 40 caracteres por línea)
  df_pos$GS_DETAILS <- str_wrap(df_pos$GS_DETAILS, width = 40)
  df_neg$GS_DETAILS <- str_wrap(df_neg$GS_DETAILS, width = 40)
  
  # --- Gráfico AUMENTADAS ---
  if(nrow(df_pos) > 0) {
    g_pos <- ggplot(df_pos, aes(x = NES, y = fct_reorder(GS_DETAILS, NES), fill = NES)) +
      geom_col(width = 0.5) +
      scale_fill_gradient(low = "indianred", high = "darkred") +
      coord_cartesian(xlim = c(0, max(df$NES, na.rm = TRUE))) +  # misma escala
      labs(
        title = paste0(comp, " AUMENTADAS"),
        x = "NES",
        y = "GENE SET"
      ) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )
    plots_aumentadas[[comp]] <- g_pos
  }
  
  # --- Gráfico DISMINUIDAS ---
  if(nrow(df_neg) > 0) {
    g_neg <- ggplot(df_neg, aes(x = abs(NES), y = fct_reorder(GS_DETAILS, abs(NES)), fill = abs(NES))) +
      geom_col(width = 0.5) +
      scale_fill_gradient(low = "lightblue", high = "darkblue") +
      coord_cartesian(xlim = c(0, max(df$NES, na.rm = TRUE))) +  # misma escala
      labs(
        title = paste0(comp, " DISMINUIDAS"),
        x = "NES",
        y = "GENE SET"
      ) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )
    plots_disminuidas[[comp]] <- g_neg
  }
}

# --- Combinar y exportar ---
# Combina todos los gráficos en un solo archivo PDF/PNG
if (length(plots_aumentadas) > 0) {
  g_final_aum <- wrap_plots(plots_aumentadas, ncol = 1)
  ggsave(paste0(output_dir, "/Comparaciones_Aumentadas.png"),
         plot = g_final_aum, width = 8, height = 2.5 * length(plots_aumentadas), dpi = 300)
  ggsave(paste0(output_dir, "/Comparaciones_Aumentadas.pdf"),
         plot = g_final_aum, width = 8, height = 2.5 * length(plots_aumentadas))
}


if (length(plots_disminuidas) > 0) {
  g_final_dis <- wrap_plots(plots_disminuidas, ncol = 1)
  ggsave(paste0(output_dir, "/Comparaciones_Disminuidas.png"),
         plot = g_final_dis, width = 8, height = 2.5 * length(plots_disminuidas), dpi = 300)
  ggsave(paste0(output_dir, "/Comparaciones_Disminuidas.pdf"),
         plot = g_final_dis, width = 8, height = 2.5 * length(plots_disminuidas))
}
