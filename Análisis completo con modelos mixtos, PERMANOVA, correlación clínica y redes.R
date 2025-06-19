# ========================
# INSTALAR Y CARGAR LIBRERÍAS
# ========================
install.packages(c(
  "readxl", "tidyverse", "broom", "pheatmap", "ggrepel",
  "vegan", "haven", "lme4", "WGCNA", "Hmisc",
  "igraph", "tidygraph", "ggraph"
))
install.packages(c("readr", "dplyr", "ggplot2"))

# Para usar SPIEC-EASI, descomenta:
# devtools::install_github("zdk123/SpiecEasi")
# install.packages("SpiecEasi")

library(readxl)
library(tidyverse)
library(broom)
library(pheatmap)
library(ggrepel)
library(vegan)
library(haven)
library(lme4)
library(WGCNA)
library(Hmisc)
library(igraph)
library(tidygraph)
library(ggraph)
library(readr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(vegan)
# library(SpiecEasi)

# ========================
# 1. CARGAR DATOS
# ========================
# Rutas funcionales
data_fun <- read_excel("Analisis_funcional.xlsx", sheet = "Sheet2")

# Variables clínicas
clinica  <- read_sav("LACTOPREM_pacientes.sav")

# Abundancias de géneros (o ASVs)
abund_gen <- read_excel("DB_lactoprem.xlsx", sheet = "DB")

# Anotaciones funcionales (KEGG/MetaCyc)
ruta_funcion <- read_csv("rutas_clasificadas_completo.csv")



# ========================
# 2. PREPARAR DATOS LONG
# ========================
data_long <- data_fun %>%
  pivot_longer(cols = -(1:3),
               names_to  = "Ruta",
               values_to = "Valor") %>%
  mutate(
    Treatment = factor(Treatment, levels = c(0,1),
                       labels = c("Placebo","Lactoferrina")),
    Time      = factor(Time, levels = c(0,1),
                       labels = c("Inicio","Fin"))
  )



# ========================
# 3. MODELO MIXTO POR RUTA
# ========================
library(readxl)
library(dplyr)
library(tidyr)         # pivot_longer()
library(broom.mixed)   # tidy(..., effects="fixed")
library(purrr)
library(lme4)

# 1) Carga
data_fun <- read_excel("Analisis funcional-29-05-25.xlsx", sheet = "Sheet2")

# 2) Extrar sujeto y pivot_longer
data_long <- data_fun %>%
  mutate(
    Subject = sub("^LAC_(\\d+)_.*", "\\1", `OTU ID`)  # captura el número de neonato
  ) %>%
  pivot_longer(
    cols = -c(`OTU ID`, Subject, Treatment, Time),
    names_to  = "Ruta",
    values_to = "Valor"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c(0,1),
                       labels = c("Placebo","Lactoferrina")),
    Time      = factor(Time, levels = c(0,1),
                       labels = c("Inicio","Fin"))
  )

glimpse(data_long)  # ahora debe tener 6 columnas: OTU ID, Subject, Treatment, Time, Ruta, Valor

# 3) Ajustr modelos mixtos por ruta usando Subject como aleatorio
modelos_mixtos <- data_long %>%
  split(.$Ruta) %>%
  imap_dfr(function(df_ruta, ruta_id) {
    # Necesitamos al menos 2 sujetos para estimar var. aleatoria
    if (n_distinct(df_ruta$Subject) < 2) return(NULL)
    
    mod <- tryCatch(
      lmer(Valor ~ Treatment * Time + (1 | Subject), data = df_ruta),
      error = function(e) NULL
    )
    if (is.null(mod)) return(NULL)
    
    broom.mixed::tidy(mod, effects = "fixed") %>%
      mutate(Ruta = ruta_id)
  })

glimpse(modelos_mixtos)  # ya debe verse term, estimate, std.error, etc.

# 4) Filtrar los términos de interés
efectos_fijos <- modelos_mixtos %>%
  filter(term %in% c(
    "TreatmentLactoferrina",
    "TimeFin",
    "TreatmentLactoferrina:TimeFin"
  ))

# 5) Guardar resultados
write.csv(efectos_fijos,
          "resultados_modelo_mixto.csv",
          row.names = FALSE)

#-------------------------------------------------------------------
#Rutas significativas, distribución de coeficientes y volcano plot
# 1) Leer resultados del modelo mixto

df <- read_csv("resultados_modelo_mixto.csv")

# 2) Filtrar solo el término de interacción
inter <- df %>%
  filter(term == "TreatmentLactoferrina:TimeFin")

# 3) Aproximar p-valores desde la t-statistic
# (supone distribución normal)
inter <- inter %>%
  mutate(
    p.value = 2 * (1 - pnorm(abs(statistic)))
  )

# 4) Contar rutas según criterio
total_routes       <- nrow(inter)
sig_routes         <- sum(inter$p.value < 0.05)
pos_sig_routes     <- sum(inter$estimate >  0 & inter$p.value < 0.05)
neg_sig_routes     <- sum(inter$estimate <  0 & inter$p.value < 0.05)
non_sig_routes     <- total_routes - sig_routes

summary_df <- tibble(
  Criterio = c(
    "Total de rutas",
    "Rutas significativas (p < 0.05)",
    "Rutas positivas y significativas",
    "Rutas negativas y significativas",
    "Rutas no significativas"
  ),
  Cantidad = c(
    total_routes,
    sig_routes,
    pos_sig_routes,
    neg_sig_routes,
    non_sig_routes
  )
)

print(summary_df)

# 5) Distribución de coeficientes
ggplot(inter, aes(x = estimate)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  theme_minimal() +
  labs(
    title = "Distribución de coeficientes (interacción)",
    x     = "Estimate (TreatmentLactoferrina:TimeFin)",
    y     = "Número de rutas"
  )
ggsave("histograma_coeficientes.png", width = 6, height = 4)

# 6) Volcano plot
threshold_p <- 0.05
ggplot(inter, aes(x = estimate, y = -log10(p.value))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_hline(yintercept = -log10(threshold_p), linetype="dashed") +
  theme_minimal() +
  labs(
    title = "Volcano plot: TreatmentLactoferrina × TimeFin",
    x     = "Estimate (interacción)",
    y     = "-log10(p-valor aproximado)"
  )
ggsave("volcano_plot.png", width = 6, height = 4)
#------------------------------------------------------------


# ========================
# 4. AGRUPACIÓN FUNCIONAL
# ========================
# Comprobar nombres de columnas
print(names(ruta_funcion))

# ==============================
# 3)Renombrar columnas para el join
# ==============================
ruta_funcion <- ruta_funcion %>%
  rename(
    Ruta    = Pathway_ID,
    Funcion = Functional_Category
  )

# Vistazo rápido
glimpse(ruta_funcion)

# ==============================
# 4) Extraer log2FC de la interacción 
# ==============================
# asume que ya tienes efectos_fijos
log2FC <- efectos_fijos %>%
  filter(term == "TreatmentLactoferrina:TimeFin") %>%
  select(Ruta, estimate) %>%
  rename(log2FC_Lactoferrina = estimate)

# ==============================
# 5) Agrupación funcional
# ==============================
resultado_funcional <- log2FC %>%
  left_join(ruta_funcion, by = "Ruta") %>%
  group_by(Funcion) %>%
  summarise(
    media_log2FC         = mean(log2FC_Lactoferrina, na.rm = TRUE),
    rutas_total          = n(),
    rutas_positivas      = sum(log2FC_Lactoferrina >  0, na.rm = TRUE),
    rutas_negativas      = sum(log2FC_Lactoferrina <  0, na.rm = TRUE),
    # Por ejemplo, rutas en cuartil superior de cambio
    rutas_cuartil_sup    = sum(
      abs(log2FC_Lactoferrina) >
        quantile(abs(log2FC_Lactoferrina), 0.75),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ==============================
# 6) Resultados
# ==============================
print(resultado_funcional)
write_csv(resultado_funcional, "agrupacion_funcional.csv")



# ========================
# 5. PCA Y PERMANOVA SOBRE RUTAS
# ========================
library(dplyr)
library(tibble)
library(vegan)

# 1) Metadata con IDs únicos (igual que antes)
meta_ruta2 <- data_fun %>%
  mutate(
    Visit    = Time + 1,
    SampleID = paste0(`OTU ID`, "_V", Visit),
    Treatment = factor(Treatment, levels = c(0,1),
                       labels = c("Placebo","Lactoferrina")),
    Time      = factor(Time, levels = c(0,1),
                       labels = c("Inicio","Fin"))
  ) %>%
  select(SampleID, Treatment, Time) %>%
  column_to_rownames("SampleID")

# 2) Matriz de abundancias con SampleID único
abund_mat <- abund_gen %>%
  mutate(
    SampleID = paste0(index, "_V", Visit)
  ) %>%
  select(-index, -Grupo, -Visit, -input_sequences) %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

# 3) Alinear filas de metadata y abundancias
common_samples <- intersect(rownames(meta_ruta2), rownames(abund_mat))
meta_ruta3     <- meta_ruta2[common_samples, , drop = FALSE]
abund_mat2     <- abund_mat[common_samples, , drop = FALSE]

# 4) En lugar de scale(), usar Hellinger para que no haya valores negativos
abund_hel <- decostand(abund_mat2, method = "hellinger")

# 5) Bray–Curtis sobre Hellinger
bray2 <- vegdist(abund_hel, method = "bray")

# 6) PERMANOVA
permanova_gen <- adonis2(
  bray2 ~ Treatment * Time,
  data         = meta_ruta3,
  permutations = 999
)

# 7) Guardar y mostrar
capture.output(permanova_gen, file = "PERMANOVA_generos.txt")
print(permanova_gen)




# =================================
# 6. RED DE CO-ACTIVACIÓN (WGCNA)
# =================================

library(WGCNA)
options(stringsAsFactors = FALSE)

# 1) Preparar la matriz de expresión (cada fila = muestra, cada columna = ruta)
datExpr <- as.data.frame(ruta_scaled)

# 2) Filtrar de muestras y rutas con NAs o varianza cero
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (any(!gsg$goodGenes)) {
    message("Quitando rutas con varianza cero o NAs: ",
            paste(names(datExpr)[!gsg$goodGenes], collapse = ", "))
  }
  if (any(!gsg$goodSamples)) {
    message("Quitando muestras con NAs: ",
            paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", "))
  }
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 3) Asegurar ejecución en un único hilo para reproducibilidad
if ("disableWGCNAThreads" %in% ls("package:WGCNA")) {
  disableWGCNAThreads()
}

# 4) Elegir soft-threshold power automáticamente
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         networkType = "signed",
                         verbose     = 5)

# 4a) Visualizar ajuste a topología libre de escala
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Signed R^2 for scale-free topology",
     type = "b")
abline(h = 0.8, col = "red", lty = 2)

# 4b) Selección automática: primer power con R^2 ≥ 0.8 (o fallback a 6)
idx <- which(sft$fitIndices[,2] >= 0.8)[1]
softPower <- if (!is.na(idx)) powers[idx] else 6
message("Chosen softPower = ", softPower)

# 5) Calcular matriz de adyacencia “suave”
adjacencyMat <- adjacency(datExpr,
                          power = softPower,
                          type  = "signed")

# 6) Calcular TOM y disimilitud
TOM     <- TOMsimilarity(adjacencyMat)
dissTOM <- 1 - TOM

# 7) Clustering jerárquico de rutas
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 8) Detección de módulos con corte dinámico
dynamicMods <- cutreeDynamic(dendro           = geneTree,
                             distM            = dissTOM,
                             deepSplit        = 2,
                             pamRespectsDendro= FALSE,
                             minClusterSize   = 30)
moduleColors <- labels2colors(dynamicMods)

# 9) Visualizar dendrograma y colores de módulos
plotDendroAndColors(geneTree, moduleColors,
                    "Módulos WGCNA",
                    dendroLabels = FALSE,
                    hang         = 0.03,
                    addGuide     = TRUE,
                    guideHang    = 0.05)

# 10) Guardar resultados
saveRDS(dynamicMods,  "WGCNA_dynamicMods.rds")
saveRDS(moduleColors, "WGCNA_moduleColors.rds")
ggsave("WGCNA_dendrograma_modulos.png", width = 8, height = 6)

# 11) Calcular eigengenes de módulo
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)

# 12) Preparar datos de traits
traitData <- meta_ruta2 %>%
  mutate(Treatment_num = as.numeric(factor(Treatment, levels = c("Placebo","Lactoferrina"))),
         Time_num      = as.numeric(factor(Time,      levels = c("Inicio","Fin")))) %>%
  select(Treatment_num, Time_num)

# 13) Correlación módulo–trait y p-values
nSamples          <- nrow(datExpr)
moduleTraitCor    <- cor(MEs, traitData, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# 14) Heatmap de correlaciones módulo–trait
textMatrix <- paste0(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")"
)
labeledHeatmap(Matrix      = moduleTraitCor,
               xLabels     = names(traitData),
               yLabels     = names(MEs),
               colorLabels = FALSE,
               colors      = blueWhiteRed(50),
               textMatrix  = textMatrix,
               setStdMargins = FALSE,
               cex.text      = 0.5,
               main          = "Module–Trait Relationships")




# =================================
# 7. CORRELACIONES CON CLÍNICA
# =================================
library(purrr)

# Asegúrate de que en 'clinica' la columna de ID coincida con OTU_ID
clinica2 <- clinica %>%
  rename(OTU_ID = ID)

pca_clin <- left_join(pca_df, clinica2, by = "OTU_ID")

# 7.1) PC1 vs IL6
cor_PC1_IL6 <- cor.test(pca_clin$PC1,
                        pca_clin$IL6,
                        use = "complete.obs",
                        method = "spearman")
print(cor_PC1_IL6)

# 7.2) Correlación PC1 con otras variables
vars_clin <- c("PCR", "Dias_UCI")  # ajusta a tus nombres reales
res_corrs <- map_dfr(vars_clin, function(var){
  tt <- cor.test(pca_clin$PC1,
                 pca_clin[[var]],
                 use = "complete.obs",
                 method = "spearman")
  tibble(
    variable = var,
    rho      = tt$estimate,
    p.value  = tt$p.value
  )
})
write.csv(res_corrs, "correlaciones_PC1_clinica.csv", row.names = FALSE)


# ==============================
# 8. GUARDAR OBJETOS IMPORTANTES
# ==============================
saveRDS(pca_df,              "pca_datos.rds")
saveRDS(permanova_rutas,     "permanova_rutas.rds")
saveRDS(permanova_gen,       "permanova_generos.rds")
saveRDS(resultado_funcional, "resultado_funcional.rds")  # si ya lo creaste antes
