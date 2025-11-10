
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
ruta <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis nuevo datos filtrados/Enrichr.xlsx"

# Leer la hoja (ajusta el nombre si es necesario)
df <- read_excel(ruta, sheet = 1, col_names = TRUE)

# --- LIMPIEZA Y PREPARACIÓN ---
# Rellenar la columna Comparacion hacia abajo si tiene huecos
df$Comparacion <- zoo::na.locf(df$Comparacion, na.rm = FALSE)
df <- df %>% filter(!is.na(Term) & Term != "")

# Asegurar formato numérico correcto (a veces los valores vienen con coma)
df <- df %>%
  mutate(
    p_value = as.numeric(gsub(",", ".", `p-value`)),
    q_value = as.numeric(gsub(",", ".", `q-value`)),
    Direction = tolower(Direction),
    log_qvalue = -log10(q_value)
  )

# Eliminar filas vacías o sin término
df <- df %>% filter(!is.na(Term) & Term != "")

# Crear carpeta de salida
output_dir <- "D:/Carpetas/Tesis/Graficos_Enrichr"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Listas para almacenar los gráficos
plots_up <- list()
plots_down <- list()

# Crear listas para almacenar los gráficos
plots <- list()

# --- BUCLE POR COMPARACIÓN ---
for (comp in unique(df$Comparacion)) {
  
  df_comp <- df %>% filter(Comparacion == comp)
  
  # ---- AUMENTADOS ----
  df_up <- df_comp %>%
    filter(Direction == "up") %>%
    arrange(q_value) %>%
    mutate(Term = str_wrap(Term, width = 60),
           Term = factor(Term, levels = rev(unique(Term))))
  
  if (nrow(df_up) > 0) {
    g_up <- ggplot(df_up, aes(x = log_qvalue, y = Term, fill = log_qvalue)) +
      geom_col(width = 0.6) +
      scale_fill_gradient(low = "#FAD4D4", high = "#B22222") +
      labs(
        title = paste0(comp, " AUMENTADOS"),
        x = "-log10(q-value)",
        y = "GENE SET"
      ) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )
    plots_up[[comp]] <- g_up
  }
  
  # ---- DISMINUIDOS ----
  df_down <- df_comp %>%
    filter(Direction == "down") %>%
    arrange(q_value) %>%
    mutate(Term = str_wrap(Term, width = 60),
           Term = factor(Term, levels = rev(unique(Term))))
  
  if (nrow(df_down) > 0) {
    g_down <- ggplot(df_down, aes(x = log_qvalue, y = Term, fill = log_qvalue)) +
      geom_col(width = 0.6) +
      scale_fill_gradient(low = "#BBDFFF", high = "#08306B") +
      labs(
        title = paste0(comp, " DISMINUIDOS"),
        x = "-log10(q-value)",
        y = "GENE SET"
      ) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )
    plots_down[[comp]] <- g_down
  }
}

# ---- Combinar y guardar ----
if (length(plots_up) > 0) {
  g_final_up <- wrap_plots(plots_up, ncol = 1)
  ggsave(paste0(output_dir, "/Enrichr_Aumentados.png"),
         plot = g_final_up, width = 9, height = 2.5 * length(plots_up), dpi = 300)
  ggsave(paste0(output_dir, "/Enrichr_Aumentados.pdf"),
         plot = g_final_up, width = 9, height = 2.5 * length(plots_up))
}

if (length(plots_down) > 0) {
  g_final_down <- wrap_plots(plots_down, ncol = 1)
  ggsave(paste0(output_dir, "/Enrichr_Disminuidos.png"),
         plot = g_final_down, width = 9, height = 2.5 * length(plots_down), dpi = 300)
  ggsave(paste0(output_dir, "/Enrichr_Disminuidos.pdf"),
         plot = g_final_down, width = 9, height = 2.5 * length(plots_down))
}



