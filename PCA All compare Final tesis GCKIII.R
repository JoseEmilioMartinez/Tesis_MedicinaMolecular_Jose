# ===== PCA PROFESIONAL COMPLETO - TU NUEVA RUTA =====

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(readxl)
library(dplyr)

# 📁 TU NUEVO ARCHIVO
ruta_datos <- "D:\\Carpetas\\Experimentos\\Celulas Endoteliales\\RNA\\RNAseq\\Resultados\\Analisis Jose 8\\All compare filtrado FPKM\\ALL_COMPARE_filtrado.xlsx"

# Carga robusta
suppressWarnings({
  datos_raw <- read_excel(ruta_datos, sheet = 1)
})

# Limpieza
datos_raw <- datos_raw[rowSums(is.na(datos_raw)) < ncol(datos_raw), ]
nombres_genes <- as.character(datos_raw[[1]])
datos_raw <- datos_raw[, -1]

suppressWarnings({
  datos_log <- as.matrix(apply(datos_raw, 2, as.numeric))
})
rownames(datos_log) <- nombres_genes

# Quitar muestras con varianza cero
varianzas <- apply(datos_log, 2, var, na.rm = TRUE)
datos_log <- datos_log[, varianzas > 0]

# PCA
res_pca <- prcomp(t(datos_log), scale. = TRUE)
muestras <- colnames(datos_log)
grupos <- strsplit(muestras, "_")
pca_data <- data.frame(
  PC1 = res_pca$x[,1], PC2 = res_pca$x[,2],
  Muestra = muestras, 
  Grupo = sapply(grupos, `[`, 1),
  Subgrupo = sapply(grupos, `[`, 2)
)

var_exp <- round(summary(res_pca)$importance[2,1:2]*100, 1)
ratio_pca <- var_exp[1] / var_exp[2]

# ===== EXACTO COMO TU IMAGEN PERO MÁS LIMPIO =====

# Tus 4 grupos
pca4 <- pca_data[pca_data$Grupo %in% c("WTG", "DKO", "KOM", "KOS"), ]
pca4$Grupo <- recode(pca4$Grupo,
                     "KOM" = "STK24 KO",
                     "KOS" = "STK25ECKO",
                     "WTG" = "STK24 WT STK25 WT",
                     "DKO" = "STK24 KO STK25ECKO") 

pca4$Grupo <- factor(pca4$Grupo, 
                     levels = c("STK24 WT STK25 WT", "STK24 KO STK25ECKO", "STK24 KO", "STK25ECKO"))

colores <- c("#1f78b4", "#ff7f00", "#e7298a", "#6a3d9a")

# EXPRESSIONES CORREGIDAS (sin guiones complejos)

nombres_super <- list(
  expression(STK24 ~ WT ~ STK25 ~ WT),
  expression(STK24 ~ KO ~ STK25^{ECKO}),
  expression(STK24 ~ KO),
  expression(STK25^{ECKO})
)

print(nombres_super)


ruta_salida <- "D:\\Carpetas\\Experimentos\\Celulas Endoteliales\\RNA\\RNAseq\\Resultados\\Graficos Final Tesis"

p_limpio <- ggplot(pca4, aes(PC1, PC2, fill = Grupo)) +
  geom_point(size = 7, shape = 21, stroke = 0.5, alpha = 1, color = "white") +  # Borde blanco fino
  geom_text_repel(aes(label = Subgrupo), size = 3.5, color = "black", fontface = "bold",
                  point.padding   = unit(1, "lines"),  # separa el número del círculo
                  box.padding     = 0.3,                 # separa etiquetas entre sí
                  max.overlaps    = Inf                  # no oculta etiquetas
  ) +
  scale_fill_manual(values = colores, 
                    labels = nombres_super,
                    name = NULL) +
  labs(title = "PCA STK24/25 vs STK24 vs STK25",
       x = paste0("PC1 (", round(var_exp[1],1), "%)"),
       y = paste0("PC2 (", round(var_exp[2],1), "%)")) +
  coord_fixed(ratio = 2/ratio_pca) +
  theme_classic(base_size = 10) +  # Clásico publicación
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5,  margin = margin(b = 25)),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 10),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 10),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        axis.line = element_blank()
  )

print(p_limpio)
ggsave(file.path(ruta_salida, "PCA_GCKIII_final_2.png"), p_limpio,  width = 9, height = 7, dpi = 400, bg = "white")
 