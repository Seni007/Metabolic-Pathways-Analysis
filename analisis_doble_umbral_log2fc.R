
# ---------------------
# Análisis funcional con comparación de dos umbrales log2FC (0.58 vs 1)
# ---------------------

# === BLOQUE PARA log2FC > 0.58 ===

alpha <- 0.05
log2fc_thresh_058 <- 0.58

resultado_058 <- resultado_final %>%
  mutate(significativo = ifelse(`p_Treatment_Time` < alpha & abs(log2FC_Lactoferrina) > log2fc_thresh_058,
                                ifelse(log2FC_Lactoferrina > 0, "Upregulated", "Downregulated"),
                                "No significativo"))

# Volcano plot 0.58
volcano_058 <- ggplot(resultado_058, aes(x = log2FC_Lactoferrina, y = -log10(`p_Treatment_Time`), color = significativo)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No significativo" = "gray")) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_vline(xintercept = c(-log2fc_thresh_058, log2fc_thresh_058), linetype = "dashed") +
  labs(title = "Volcano plot (log₂FC > 0.58)", x = "log₂FC Lactoferrina", y = "-log₁₀(p)") +
  theme_minimal()

ggsave("volcano_plot_058.png", plot = volcano_058, width = 8, height = 6)

# Heatmap 0.58
top_rutas_058 <- resultado_058 %>%
  filter(significativo != "No significativo") %>%
  arrange(`p_Treatment_Time`) %>%
  slice_head(n = 20) %>%
  pull(Ruta)

matriz_058 <- data %>%
  select(`#OTU ID`, Treatment, Time, all_of(top_rutas_058)) %>%
  column_to_rownames("#OTU ID") %>%
  select(-Treatment, -Time)

png("heatmap_058.png", width = 1000, height = 800)
pheatmap(scale(matriz_058),
         main = "Heatmap (log₂FC > 0.58) - Top 20 rutas")
dev.off()

# === BLOQUE PARA log2FC > 1 ===

log2fc_thresh_1 <- 1

resultado_1 <- resultado_final %>%
  mutate(significativo = ifelse(`p_Treatment_Time` < alpha & abs(log2FC_Lactoferrina) > log2fc_thresh_1,
                                ifelse(log2FC_Lactoferrina > 0, "Upregulated", "Downregulated"),
                                "No significativo"))

# Volcano plot 1
volcano_1 <- ggplot(resultado_1, aes(x = log2FC_Lactoferrina, y = -log10(`p_Treatment_Time`), color = significativo)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No significativo" = "gray")) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_vline(xintercept = c(-log2fc_thresh_1, log2fc_thresh_1), linetype = "dashed") +
  labs(title = "Volcano plot (log₂FC > 1)", x = "log₂FC Lactoferrina", y = "-log₁₀(p)") +
  theme_minimal()

ggsave("volcano_plot_1.png", plot = volcano_1, width = 8, height = 6)

# Heatmap 1
top_rutas_1 <- resultado_1 %>%
  filter(significativo != "No significativo") %>%
  arrange(`p_Treatment_Time`) %>%
  slice_head(n = 20) %>%
  pull(Ruta)

matriz_1 <- data %>%
  select(`#OTU ID`, Treatment, Time, all_of(top_rutas_1)) %>%
  column_to_rownames("#OTU ID") %>%
  select(-Treatment, -Time)

png("heatmap_1.png", width = 1000, height = 800)
pheatmap(scale(matriz_1),
         main = "Heatmap (log₂FC > 1) - Top 20 rutas")
dev.off()
