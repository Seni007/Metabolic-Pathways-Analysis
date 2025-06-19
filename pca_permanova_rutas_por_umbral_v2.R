
# Librerías necesarias
library(tidyverse)
library(pheatmap)
library(vegan)

# Establecer carpeta de salida
setwd("/Users/inesmaciasgamero/Desktop/MAster/TFM/Data/output")

# --- PCA + PERMANOVA para log2FC > 0.58 ---

top_rutas_058 <- resultado_final %>%
  filter(`p_Treatment_Time` < 0.05, abs(log2FC_Lactoferrina) > 0.58) %>%
  pull(Ruta)

matriz_058 <- data %>%
  select(`#OTU ID`, Treatment, Time, all_of(top_rutas_058)) %>%
  column_to_rownames("#OTU ID")

metadata_058 <- matriz_058 %>%
  select(Treatment, Time)

matriz_058_valores <- matriz_058 %>%
  select(-Treatment, -Time)

pca_058 <- prcomp(scale(matriz_058_valores), center = TRUE)
pca_df_058 <- as.data.frame(pca_058$x[, 1:2])
pca_df_058$Treatment <- factor(metadata_058$Treatment, levels = c(0, 1), labels = c("Placebo", "Lactoferrina"))
pca_df_058$Time <- factor(metadata_058$Time, levels = c(0, 1), labels = c("Inicio", "Fin"))

p_058 <- ggplot(pca_df_058, aes(x = PC1, y = PC2, color = Treatment, shape = Time)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Treatment), type = "norm", linetype = "dashed") +
  labs(title = "PCA (log₂FC > 0.58)", x = "PC1", y = "PC2") +
  theme_minimal()

ggsave("pca_058.png", plot = p_058, width = 8, height = 6)

dist_058 <- vegdist(scale(matriz_058_valores), method = "bray")
permanova_058 <- adonis2(dist_058 ~ Treatment * Time, data = metadata_058, permutations = 999)
capture.output(permanova_058, file = "permanova_058.txt")

# --- PCA + PERMANOVA para log2FC > 1 ---

top_rutas_1 <- resultado_final %>%
  filter(`p_Treatment_Time` < 0.05, abs(log2FC_Lactoferrina) > 1) %>%
  pull(Ruta)

matriz_1 <- data %>%
  select(`#OTU ID`, Treatment, Time, all_of(top_rutas_1)) %>%
  column_to_rownames("#OTU ID")

metadata_1 <- matriz_1 %>%
  select(Treatment, Time)

matriz_1_valores <- matriz_1 %>%
  select(-Treatment, -Time)

pca_1 <- prcomp(scale(matriz_1_valores), center = TRUE)
pca_df_1 <- as.data.frame(pca_1$x[, 1:2])
pca_df_1$Treatment <- factor(metadata_1$Treatment, levels = c(0, 1), labels = c("Placebo", "Lactoferrina"))
pca_df_1$Time <- factor(metadata_1$Time, levels = c(0, 1), labels = c("Inicio", "Fin"))

p_1 <- ggplot(pca_df_1, aes(x = PC1, y = PC2, color = Treatment, shape = Time)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(x = PC1, y = PC2, group = Treatment, color = Treatment), type = "norm", linetype = "dashed") +
  labs(title = "PCA (log₂FC > 1)", x = "PC1", y = "PC2") +
  theme_minimal()

ggsave("pca_1.png", plot = p_1, width = 8, height = 6)

dist_1 <- vegdist(scale(matriz_1_valores), method = "bray")
permanova_1 <- adonis2(dist_1 ~ Treatment * Time, data = metadata_1, permutations = 999)
capture.output(permanova_1, file = "permanova_1.txt")
