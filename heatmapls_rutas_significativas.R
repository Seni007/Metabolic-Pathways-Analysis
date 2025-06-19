# Cargar librerías
library(readxl)
library(dplyr)
library(pheatmap)

# Leer datos
df <- read_excel("Analisis_funcional.xlsx", sheet = "Sheet2")

# Crear BaseID
df$BaseID <- gsub("_T[03]_H", "", df$`OTU ID`)

# Convertir a data.frame (para evitar problemas de tibble)
df <- as.data.frame(df)

# Separar t0 y t1
t0 <- df[df$Time == 0, ]
t1 <- df[df$Time == 1, ]

# Encontrar columnas funcionales (todas las numéricas excepto Time y Treatment)
exclude_cols <- c("Time", "Treatment")
func_cols <- names(df)[sapply(df, is.numeric) & !names(df) %in% exclude_cols]

# Asegurar que emparejamos por BaseID
t0 <- t0[match(t1$BaseID, t0$BaseID), ]

# Calcular delta por ruta
delta <- t1[, func_cols] - t0[, func_cols]
rownames(delta) <- t1$BaseID

# Generar heatmap 
pheatmap(as.matrix(delta),
         main = "Delta funcional (Fin - Inicio)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 3)


rutas_significativas <- c("FUC-RHAMCAT-PWY", "PWY-7352", "PWY-6486", "PWY-5737", "PWY-6760", "PWY-5675", "NPGLUCAT-PWY", "GLUCUROCAT-PWY", "PWY-2201", "PWY-7022", "PWY-7323", "PWY-6395", "PWY-7483", "RHAMCAT-PWY", "PWY-5634", "PWY-5679", "GALACT-GLUCUROCAT-PWY", "PWY-7242", "PWY-5659", "GALACTUROCAT-PWY", "PWY-7456", "PWY-6507")
# Filtrar delta solo para rutas significativas
delta_signif <- delta[, rutas_significativas, drop = FALSE]

# Extraer solo esas columnas
delta_signif <- delta[, rutas_significativas, drop = FALSE]

# Evaluar variabilidad (SD) y seleccionar top 20 rutas más variables
variabilidad <- apply(delta_signif, 2, sd)
rutas_top20 <- names(sort(variabilidad, decreasing = TRUE))[1:20]

# ---- HEATMAP RUTAS SIGNIFICATIVAS (Top 20) ----
pheatmap(delta_signif[, rutas_top20],
         main = "Δ funcional — Top 20 rutas p < 0.05",
         scale = "row",
         color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6,
         fontsize_col = 5,
         border_color = NA,           # Sin bordes entre celdas
         cellwidth = NA,              # Deja que el ancho se ajuste automáticamente
         cellheight = NA,             # Altura automática para legibilidad
         treeheight_row = 20,         # Reduce tamaño del dendrograma de filas
         treeheight_col = 20,         # Idem para columnas
         legend = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         angle_col = 45)              # Inclina etiquetas de columna para legibilidad


rutas_significativas_variacion <- c("FUC-RHAMCAT-PWY", "GLUCUROCAT-PWY", "PWY-7323", "RHAMCAT-PWY", "PWY-5634", "GALACT-GLUCUROCAT-PWY", "PWY-7242", "PWY-5659", "GALACTUROCAT-PWY", "PWY-7456", "PWY-6507")
# Filtrar delta solo para rutas significativas
delta_signif_v <- delta[, rutas_significativas_variacion, drop = FALSE]

# Extraer solo esas columnas
delta_signif_v <- delta[, rutas_significativas_variacion, drop = FALSE]

# Evaluar variabilidad (SD) y seleccionar top 10 rutas más variables
variabilidad_v <- apply(delta_signif_v, 2, sd)
rutas_top20_v <- names(sort(variabilidad_v, decreasing = TRUE))[1:20]

# ---- HEATMAP RUTAS SIGNIFICATIVAS (Top 10 var) ----
# Tu lista ya está filtrada manualmente, no hace falta extraer top 10
pheatmap(delta_signif_v,
         main = "Δ funcional — rutas seleccionadas p < 0.05",
         scale = "row",
         color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 3,
         fontsize_col = 5,
         border_color = NA,
         cellwidth = NA,
         cellheight = NA,
         treeheight_row = 20,
         treeheight_col = 20,
         legend = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         angle_col = 45)
