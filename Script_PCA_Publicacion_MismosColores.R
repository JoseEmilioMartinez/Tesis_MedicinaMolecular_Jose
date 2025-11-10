install.packages(c("ggplot2", "FactoMineR", "factoextra"))


# Cargar librerías
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(RColorBrewer)

# Transponer matriz: ahora filas = muestras, columnas = genes
pca_matrix <- t(datos_log)

# Crear los nombres de muestra desde las columnas del dataset
muestras <- colnames(datos_log)

# Separar grupo y subgrupo a partir del nombre, por ejemplo: "DKO_1"
grupos_split <- strsplit(muestras, "_")
grupos_df <- data.frame(
  Muestra = muestras,
  Grupo = sapply(grupos_split, `[`, 1),
  Subgrupo = sapply(grupos_split, `[`, 2)
)

# Realizar PCA
res_pca <- prcomp(t(datos_log), scale. = TRUE)

# 2. Extraer scores y construir data.frame
pca_data <- as.data.frame(res_pca$x)

# 3. Extraer grupo y subgrupo desde los nombres de las columnas
muestras <- colnames(datos_log)
grupos_split <- strsplit(muestras, "_")

# Crear tabla con anotaciones
grupos_df <- data.frame(
  Muestra = muestras,
  Grupo = sapply(grupos_split, `[`, 1),
  Subgrupo = sapply(grupos_split, `[`, 2)
)

# Añadir info al PCA
pca_data$Grupo <- grupos_df$Grupo
pca_data$Subgrupo <- grupos_df$Subgrupo

# 4. Porcentaje de varianza explicada
var_explicada <- summary(res_pca)$importance[2, 1:2] * 100

# 5. PCA plot elegante y limpio
ggplot(pca_data, aes(x = PC1, y = PC2, fill = Grupo, label = Subgrupo)) +
  geom_point(shape = 21, color = "black", size = 5, stroke = 1.2) +  # puntos con borde negro
  geom_text_repel(size = 3.5, color = "black", max.overlaps = Inf) +  # números de subgrupo
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Análisis PCA WTGvsWTPvsDKOvsKOP",
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


##pca con los mismo colores que Heatmap
# Añadir info al PCA
pca_data$Grupo <- factor(pca_data$Grupo, levels = names(grupo_colores))
pca_data$Subgrupo <- grupos_df$Subgrupo
unique(pca_data$Grupo)
table(pca_data$Grupo, useNA = "ifany")

# Porcentaje de varianza explicada
var_explicada <- summary(res_pca)$importance[2, 1:2] * 100

# PCA plot elegante y limpio con colores definidos
ggplot(pca_data, aes(x = PC1, y = PC2, fill = Grupo, label = Subgrupo)) +
  geom_point(shape = 21, color = "black", size = 5, stroke = 1.2) +  # puntos con borde negro
  geom_text_repel(size = 3.5, color = "black", max.overlaps = Inf) +  # números de subgrupo
  scale_fill_manual(values = grupo_colores) +
  scale_color_manual(values = grupo_colores) +
  labs(
    title = "Análisis PCA GCKIII",
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

