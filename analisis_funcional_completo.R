
# ---------------------
# Análisis funcional + visualización
# ---------------------

# Cargar librerías necesarias
install.packages("readxl")
library(readxl)
install.packages("tidyverse")
library(tidyverse)
install.packages("broom")
library(broom)
install.packages("pheatmap")
library(pheatmap)
install.packages("ggrepel")
library(ggrepel)

# Establecer el directorio de trabajo
setwd("/Users/inesmaciasgamero/Desktop/MAster/TFM/Data")

# Leer el archivo de datos
data <- read_excel("Analisis funcional-29-05-25.xlsx", sheet = "Sheet2")

# Transformar a formato largo
data_long <- data %>%
  pivot_longer(cols = -(1:3), names_to = "Ruta", values_to = "Valor")

# Convertir factores
data_long$Treatment <- factor(data_long$Treatment, levels = c(0, 1), labels = c("Placebo", "Lactoferrina"))
data_long$Time <- factor(data_long$Time, levels = c(0, 1), labels = c("Inicio", "Fin"))

# Modelo ANOVA por ruta
resultados <- data_long %>%
  group_by(Ruta) %>%
  do(tidy(aov(Valor ~ Treatment * Time, data = .))) %>%
  ungroup()

# Filtrar términos relevantes
resultados_filtrados <- resultados %>%
  filter(term %in% c("Treatment", "Time", "Treatment:Time")) %>%
  pivot_wider(names_from = term, values_from = c(statistic, p.value), names_sep = "_") %>%
  rename_with(~ gsub("p.value_", "p_", .x))

# Calcular log2FC entre Fin e Inicio por grupo
resumen_cambio <- data_long %>%
  group_by(Ruta, Treatment, Time) %>%
  summarise(media = mean(Valor), .groups = "drop") %>%
  pivot_wider(names_from = c("Treatment", "Time"), values_from = "media") %>%
  mutate(log2FC_Lactoferrina = log2(`Lactoferrina_Fin` / `Lactoferrina_Inicio`),
         log2FC_Placebo = log2(`Placebo_Fin` / `Placebo_Inicio`))

# Juntar todo
resultado_final <- left_join(resultados_filtrados, resumen_cambio, by = "Ruta")

# Guardar CSV con resultados
write.csv(resultado_final, "resultados_analisis_funcional.csv", row.names = FALSE)

# --- VOLCANO PLOT ---
#log2fc = 0.58
alpha <- 0.05
log2fc_thresh <- 0.58  # cambio mínimo equivalente a fold change de 1.5x

# Asegurar que el nombre de la columna sea correcto
resultado_final <- resultado_final %>%
  rename(p_Treatment_Time = `p_Treatment:Time`)

# Clasificar las rutas
resultado_final <- resultado_final %>%
  mutate(significativo = ifelse(p_Treatment_Time < alpha & abs(log2FC_Lactoferrina) > log2fc_thresh,
                                ifelse(log2FC_Lactoferrina > 0, "Upregulated", "Downregulated"),
                                "No significativo"))

# Eliminar rutas con valores faltantes
volcano_data <- resultado_final %>%
  filter(!is.na(log2FC_Lactoferrina), !is.na(p_Treatment_Time))

# Generar volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FC_Lactoferrina, y = -log10(p_Treatment_Time), color = significativo)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No significativo" = "gray")) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_vline(xintercept = c(-log2fc_thresh, log2fc_thresh), linetype = "dashed") +
  labs(title = "Volcano plot: Rutas funcionales (log₂FC > 0.58 y p < 0.05)",
       x = "log₂(Fold Change) Lactoferrina",
       y = "-log₁₀(p-valor Interacción)") +
  theme_minimal()

# Guardar como imagen
ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6)

# --- VOLCANO PLOT ---
# log2fc 1
alpha <- 0.05
log2fc_thresh <- 1

resultado_final <- resultado_final %>%
  rename(p_Treatment_Time = p_Treatment_Time)

resultado_final <- resultado_final %>%
  mutate(significativo = ifelse(p_Treatment_Time < alpha & abs(log2FC_Lactoferrina) > log2fc_thresh,
                                ifelse(log2FC_Lactoferrina > 0, "Upregulated", "Downregulated"),
                                "No significativo"))

volcano_data <- resultado_final %>%
  filter(!is.na(log2FC_Lactoferrina), !is.na(p_Treatment_Time))

volcano_plot <- ggplot(resultado_final, aes(x = log2FC_Lactoferrina, y = -log10(p_Treatment_Time), color = significativo)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No significativo" = "gray")) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_vline(xintercept = c(-log2fc_thresh, log2fc_thresh), linetype = "dashed") +
  labs(title = "Volcano plot: Rutas funcionales",
       x = "log2(Fold Change) Lactoferrina",
       y = "-log10(p-valor Interacción)") +
  theme_minimal()

ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6)


# --- HEATMAP ---
top_rutas <- resultado_final %>%
  filter(significativo != "No significativo") %>%
  arrange(p_Treatment_Time) %>%
  slice_head(n = 20) %>%
  pull(Ruta)

datos_matrix <- data %>%
  select(`#OTU ID`, Treatment, Time, all_of(top_rutas)) %>%
  column_to_rownames("#OTU ID") %>%
  select(-Treatment, -Time)

# Mostrar en pantalla
heatmap_obj <- pheatmap(scale(datos_matrix), 
                        main = "Heatmap - Top 20 rutas significativas")

# Guardar como imagen PNG
png("heatmap_rutas1.png", width = 1000, height = 800)
pheatmap(scale(datos_matrix), 
         main = "Heatmap - Top 20 rutas significativas")
dev.off()


# --- PCA y PERMANOVA sobre rutas funcionales ---

# Crear matriz de rutas funcionales
ruta_matrix <- data %>%
  column_to_rownames("#OTU ID") %>%
  select(-Treatment, -Time)

# Extraer metadata
metadata <- data %>%
  select(`#OTU ID`, Treatment, Time) %>%
  column_to_rownames("#OTU ID")

# Escalar los datos
ruta_scaled <- scale(ruta_matrix)

# PCA
pca <- prcomp(ruta_scaled, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])

# Añadir metadatos y convertir a factores
pca_df$Treatment <- factor(metadata$Treatment, levels = c(0, 1), labels = c("Placebo", "Lactoferrina"))
pca_df$Time <- factor(metadata$Time, levels = c(0, 1), labels = c("Inicio", "Fin"))

# Gráfico PCA (mostrar en pantalla)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Time)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA de rutas funcionales",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal()

# Mostrar el gráfico en RStudio
print(pca_plot)

# Guardar imagen
ggsave("pca_rutas.png", plot = pca_plot, width = 8, height = 6)

# PERMANOVA (adonis)
install.packages("vegan")
library(vegan)

# Calcular distancia Bray-Curtis
bray_dist <- vegdist(ruta_scaled, method = "bray")

# Añadir tratamiento y tiempo como factores
permanova_result <- adonis2(bray_dist ~ Treatment * Time, data = metadata, permutations = 999)

# Guardar resultados PERMANOVA
capture.output(permanova_result, file = "permanova_resultados.txt")
