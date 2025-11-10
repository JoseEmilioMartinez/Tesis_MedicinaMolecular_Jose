##instalar paquetes
install.packages(c("ggplot2", "FactoMineR", "factoextra", "ggrepel", "RColorBrewer", "readxl", "dplyr", "tibble"))

##cargar librerias

library(ggplot2)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(tibble)

# 1. Leer el Excel (ajusta sheet si es necesario)
ruta_archivo <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Nuevo analisis Jose/All compare/All_compare_fpkmn.xlsx"
datos_raw <- read_excel(ruta_archivo)

# 2. Seleccionar sólo columnas de expresión (por ejemplo FPKM)

expr_cols <- grep("_fpkm$", colnames(datos_raw), value = TRUE)
datos_expr <- datos_raw %>% select(all_of(expr_cols))

# 3. Convertir comas decimales a puntos y pasar a numérico
datos_expr <- datos_raw %>% 
  select(all_of(expr_cols)) %>% 
  mutate(across(everything(), ~as.numeric(gsub(",", ".", .))))

# 4. Filtrar genes con mayor varianza (11000 más variables)
varianzas <- apply(datos_expr, 1, var, na.rm = TRUE)
genes_mas_variables <- order(varianzas, decreasing = TRUE)[1:11000]
datos_expr_filtrado <- datos_expr[genes_mas_variables, ]

# 5. Transformar log2(x+1)
datos_log <- log2(datos_expr_filtrado + 1)

# 6. Nombres de muestras desde las columnas
muestras <- colnames(datos_log)

# 7. Extraer grupo y réplica del nombre (antes y después del "_")
grupos_split <- strsplit(gsub("_fpkm$", "", muestras), "_")
grupos_df <- data.frame(
  Muestra = muestras,
  Grupo = sapply(grupos_split, `[`, 1),   # WTG, KOS, KOM, DKO, WTP, KOP
  Replica = sapply(grupos_split, `[`, 2)  # 1,2,3,4
)

# 8. PCA (transponemos: filas = muestras, columnas = genes)
res_pca <- prcomp(t(datos_log), scale. = TRUE)
pca_data <- as.data.frame(res_pca$x)
pca_data$Grupo <- grupos_df$Grupo
pca_data$Replica <- grupos_df$Replica

# 9. Filtrar SOLO los grupos de interés
grupos_interes <- c("WTG", "WTP", "KOP", "DKO") 

pca_data <- pca_data %>% filter(Grupo %in% grupos_interes)

# 10. Porcentaje varianza explicada
var_explicada <- summary(res_pca)$importance[2, 1:2] * 100

# 11. Definir paleta de colores personalizada
grupo_colores <- c(
  "WTG" = "#1f78b4",   # azul fuerte
  "WTP" = "#e41a1c",   # rojo intenso
  "KOS" = "#6a3d9a",   # morado
  "KOM" = "#e7298a",   # rosa fuerte
  "DKO" = "#ff7f00",   # naranja fuerte
  "KOP" = "#33a02c"    # verde fuerte
)


# 11. Gráfico PCA
ggplot(pca_data, aes(x = PC1, y = PC2, fill = Grupo, label = Replica)) +
  geom_point(shape = 21, color = "black", size = 5, stroke = 1.2) +
  geom_text_repel(size = 3.5, color = "black", max.overlaps = Inf) +
  scale_fill_manual(values = colores_grupos) +
  scale_color_manual(values = colores_grupos) +
  labs(
    title = "PCA WTGvsWTPvsDKOvsKOP",
    x = paste0("PC1 (", round(var_explicada[1], 1), "% varianza)"),
    y = paste0("PC2 (", round(var_explicada[2], 1), "% varianza)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )



