
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages(c("readxl", "dplyr", "ggplot2", "openxlsx"))

# 1️⃣ Cargar librerías
library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)   # Cambiar a org.Mm.eg.db si tus datos son de ratón
library(enrichplot)
library(ggplot2)

# 2️⃣ Cargar archivo Excel
deg_data <- read_excel("D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis nuevo datos filtrados/KOPvsWTP filtrado.xlsx")

# 2️⃣ Reemplazar comas por puntos y convertir a numérico
numeric_cols <- names(deg_data)[1:(ncol(deg_data)-1)]  # todas menos gene_name

deg_data <- deg_data %>%
  mutate(across(all_of(numeric_cols),
                ~as.numeric(gsub(",", ".", as.character(.)))))


# 1️⃣ Extraer columna de genes
genes_symbols <- deg_data$gene_name

# 2️⃣ Limpiar espacios y duplicados
genes_symbols <- unique(trimws(genes_symbols))

# 3️⃣ Convertir a ENTREZ ID
genes_entrez <- bitr(genes_symbols,
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Mm.eg.db)

# 4️⃣ Ver primeros genes convertidos
head(genes_entrez)

# Enriquecimiento KEGG
kegg_res <- enrichKEGG(gene         = genes_entrez$ENTREZID,
                       organism     = "mmu",  # 'mmu' para ratón, 'hsa' para humano
                       pAdjustMethod= "BH",
                       pvalueCutoff = 0.05)

# Ver los resultados
head(kegg_res)

barplot(kegg_res, showCategory = 10) + ggtitle("KEGG Enrichment - DKOvsWTG")


# Visualización: barplot
p <- barplot(kegg_res, showCategory = 10)   # número de vías a mostrar
                 

# Personalización de colores y título
p <- p + aes(fill = p.adjust) + 
  scale_fill_gradient(low = "red", high = "blue", name = "p.adj") + 
  ggtitle("KEGG Enrichment - KOSvsWTG") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )
# 🔹 Ruta de la carpeta donde se guardarán los gráficos
output_dir <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis nuevo datos filtrados/graficos_KEGG"

# Crear carpeta si no existe
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 🔹 Nombre de archivos
pdf_file <- file.path(output_dir, "KEGG_barplot.pdf")
png_file <- file.path(output_dir, "KEGG_barplot.png")

# 🔹 Guardar PDF
ggsave(pdf_file, plot = p, device = "pdf", width = 10, height = 8)

# 🔹 Guardar PNG
ggsave(png_file, plot = p, device = "png", width = 10, height = 8, dpi = 300)

##Separacion aumentados y disminuidos
# 1️⃣ Realizar análisis KEGG
kegg_res <- enrichKEGG(gene         = genes_entrez$ENTREZID,
                       organism     = "mmu",  # Cambiar a "hsa" si es humano
                       pAdjustMethod= "BH",
                       pvalueCutoff = 0.05)

# 2️⃣ Convertir a data.frame para manipular
kegg_df <- as.data.frame(kegg_res)


# Ver estructura para asegurarnos
str(kegg_df)

# Forzar todas las columnas simples a tipo base
kegg_df <- kegg_df %>%
  mutate(
    Description = as.character(Description),
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust)
  )

# A veces el error persiste por columnas tipo list: las eliminamos
kegg_df <- kegg_df %>% select(where(~!is.list(.)))

# Verificamos que ahora sea un data.frame simple
str(kegg_df)

# 1️⃣ Asegurarte de que p.adjust sea numérico
kegg_df <- kegg_df %>%
  dplyr::mutate(
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust)
  )

# Visualización: barplot
barplot(kegg_res, showCategory = 10) + ggtitle("KEGG Enrichment - DKOvsWTG")



# 3️⃣ Gráfico con gradiente corregido (rojos)
ggplot(top_up, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(
    low = "#FFD1D1",  # rosa claro
    high = "#B30000", # rojo oscuro
    name = "p.adjust"
  ) +
  geom_text(aes(label = Count), hjust = -0.2, size = 4) +
  labs(title = "Top 5 vías KEGG aumentadas", x = "", y = "Número de genes") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  ylim(0, max(top_up$Count) * 1.2)

# 4️⃣ Gráfico con gradiente corregido (azules)
ggplot(top_down, aes(x = reorder(Description, -Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(
    low = "#D1E5FF",  # azul claro
    high = "#003366", # azul oscuro
    name = "p.adjust"
  ) +
  geom_text(aes(label = Count), hjust = -0.2, size = 4) +
  labs(title = "Top 5 vías KEGG disminuidas", x = "", y = "Número de genes") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  ylim(0, max(top_down$Count) * 1.2)

# 7️⃣ Exportar resultados a Excel (opcional)
library(openxlsx)
write.xlsx(as.data.frame(kegg_res), 
           file = "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis nuevo datos filtrados/KEGG_DKOvsWTG.xlsx",
           overwrite = TRUE)