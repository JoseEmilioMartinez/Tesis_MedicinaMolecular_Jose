# Cargar librerías
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidyr)
library(zoo)


# Leer el archivo Excel
ruta <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Otros analisis/Analisis nuevo datos filtrados/GSEA/All Compare Filtrado/Genesets GSEA.xlsx"

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

# 🎨 FUNCIÓN PARA LIMPIAR NOMBRES (MUY LEGIBLES)
limpiar_nombre <- function(x) {
  x %>% 
    str_replace_all("^GOBP_", "GO: ") %>%  # GOBP_ → GO:
    str_replace_all("_", " ") %>%          # Guiones bajos → espacios
    str_replace_all("GOBP", "") %>%        # Quitar GOBP residual
    str_to_title() %>%                     # Primera letra mayúscula
    str_trim()                             # Limpiar espacios
}

df$GS_LABEL <- limpiar_nombre(df$GS_DETAILS)

# Crear carpeta de salida
output_dir <- "D:/Carpetas/Tesis/Graficos_GSEA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 🎨 BUCLE PRINCIPAL CON LEYENDA DE VALORES
for (comp in unique(df$Comparacion)) {
  
  df_comp <- df %>% filter(Comparacion == comp)
  
  df_pos <- df_comp %>% filter(NES > 0) %>% slice_max(NES, n = 5, with_ties = FALSE)
  df_neg <- df_comp %>% filter(NES < 0) %>% slice_max(abs(NES), n = 5, with_ties = FALSE)
  
  # AUMENTADAS - PROFESIONAL
  if(nrow(df_pos) > 0) {
    g_pos <- ggplot(df_pos, aes(x = NES, y = fct_reorder(GS_LABEL, NES), fill = NES)) +
      geom_col(width = 0.65, alpha = 0.85, color = "white", linewidth = 0.3) +
      scale_fill_gradient(low = "pink", high = "red", name = "NES") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(title = paste(comp, "Aumentadas", sep = "\n"), 
           x = "NES", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 11, hjust = 1, face = "bold", margin = margin(r = 8)),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.margin = margin(15, 55, 15, 20),  # Más espacio para valores
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    ggsave(file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", comp), "_Aumentadas.png")),
           g_pos, width = 11, height = 5.2, dpi = 400, bg = "white")
    ggsave(file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", comp), "_Aumentadas.pdf")),
           g_pos, width = 11, height = 5.2, bg = "white")
  }
  
  # DISMINUIDAS - PROFESIONAL
  if(nrow(df_neg) > 0) {
    g_neg <- ggplot(df_neg, aes(x = abs(NES), y = fct_reorder(GS_LABEL, NES), fill = abs(NES))) +
      geom_col(width = 0.65, alpha = 0.85, color = "white", linewidth = 0.3) +
      scale_fill_gradient(low = "blue", high = "lightblue", name = "|NES|") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(title = paste(comp, "Disminuidas", sep = "\n"), 
           x = "|NES|", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 11, hjust = 1, face = "bold", margin = margin(r = 8)),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.height = unit(0.8, "cm"),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.margin = margin(15, 55, 15, 20),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    ggsave(file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", comp), "_Disminuidas.png")),
           g_neg, width = 11, height = 5.2, dpi = 400, bg = "white")
    ggsave(file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", comp), "_Disminuidas.pdf")),
           g_neg, width = 11, height = 5.2, bg = "white")
  }
}

cat("¡VERSIÓN PROFESIONAL sin ggtext lista en:", output_dir, "\n")

######################################

# Crear columna indicando dirección
df <- df %>%
  mutate(Direccion = ifelse(NES > 0, "Aumentado", "Disminuido"))

# Crear carpeta de salida
output_dir <- "D:/Carpetas/Tesis/Graficos_GSEA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 🎨 Bucle por cada comparación
for (comp in unique(df$Comparacion)) {
  
  df_comp <- df %>% filter(Comparacion == comp)
  df_pos <- df_comp %>% filter(NES > 0) %>% slice_max(order_by = abs(NES), n = 5)
  df_neg <- df_comp %>% filter(NES < 0) %>% slice_max(order_by = abs(NES), n = 5)  # usar NES_abs para misma dirección
  
  # --- Gráfico AUMENTADAS ---
  if(nrow(df_pos) > 0) {
    g_pos <- ggplot(df_pos, aes(x = abs(NES), y = fct_reorder(GS_DETAILS, NES), fill = abs(NES))) +
      geom_col() +
      scale_fill_gradient(low = "pink", high = "red") +
      labs(
        title = paste0(comp, " AUMENTADAS"),
        x = "NES",
        y = "Gene Set"
      ) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold")
      )
    
    # Guardar gráfico
    ggsave(filename = paste0(output_dir, "/", comp, "_Aumentadas.png"),
           plot = g_pos, width = 8, height = 5, dpi = 300)
  }
  
  # --- Gráfico DISMINUIDAS ---
  if(nrow(df_neg) > 0) {
    g_neg <- ggplot(df_neg, aes(x = abs(NES), 
                                y = fct_reorder(GS_DETAILS, NES),  # usar NES original
                                fill = abs(NES))) +
      geom_col() +
      scale_fill_gradient(low = "lightblue", high = "blue") +
      labs(
        title = paste0(comp, " DISMINUIDAS"),
        x = "NES (abs)",
        y = "GENE SET"
      ) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold")
      )
    
    # Guardar gráfico
    ggsave(filename = paste0(output_dir, "/", comp, "_Disminuidas.png"),
           plot = g_neg, width = 8, height = 5, dpi = 300)
  }
}

# Bucle por cada comparación
for (comp in unique(df$Comparacion)) {
  
  df_comp <- df %>% filter(Comparacion == comp)
  
  # Top 5 aumentadas y top 5 disminuidas por magnitud
  df_pos <- df_comp %>% filter(NES > 0) %>% slice_max(order_by = NES, n = 5)  
  df_neg <- df_comp %>% filter(NES < 0) %>% slice_max(order_by = abs(NES), n = 5) # NES negativo más pequeño = más disminuid
  # --- Gráfico AUMENTADAS ---
  if(nrow(df_pos) > 0) {
    g_pos <- ggplot(df_pos, aes(x = NES, y = fct_reorder(GS_DETAILS, NES), fill = NES)) +
      geom_col(width = 0.5) +
      scale_fill_gradient(low = "pink", high = "red") +
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
    
    ggsave(filename = paste0(output_dir, "/", comp, "_Aumentadas.png"),
           plot = g_pos, width = 8, height = 5, dpi = 300)
    # Guardar PDF
    ggsave(filename = paste0(output_dir, "/", comp, "_Aumentadas.pdf"),
           plot = g_pos, width = 8, height = 5)
  }
  
  # --- Gráfico DISMINUIDAS ---
  if(nrow(df_neg) > 0) {
    # Usamos abs(NES) para la longitud de la barra, pero fct_reorder con NES original para ordenar
    g_neg <- ggplot(df_neg, aes(x = abs(NES), y = fct_reorder(GS_DETAILS, NES), fill = NES)) +
      geom_col(width = 0.5) +
      scale_fill_gradient(low = "blue", high = "lightblue") +
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
    
    ggsave(filename = paste0(output_dir, "/", comp, "_Disminuidas.png"),
           plot = g_neg, width = 8, height = 5, dpi = 300)
    # Guardar PDF
    ggsave(filename = paste0(output_dir, "/", comp, "_Disminuidas.pdf"),
           plot = g_neg, width = 8, height = 5)
  }
}