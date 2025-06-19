# Directorio de trabajo
setwd("/Users/inesmaciasgamero/Desktop/MAster/TFM/Data")

# Cargar librerías
library(readxl)
library(tidyverse)
library(haven)
library(janitor)
library(lme4)
library(broom.mixed)

# Leer archivo funcional
fun <- read_xlsx("Analisis funcional-29-05-25.xlsx")

# Suma total por muestra
sums <- fun |>
  select(where(is.numeric)) |>
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "Sample", values_to = "total")

summary(sums$total)
hist(sums$total, breaks = 30)

# Leer datos clínicos
clin <- read_sav("LACTOPREM JULIO 14_08_2024 pacientes definitivos.sav") |>
  clean_names()

# Leer microbiota
otu_raw <- read_xlsx("DB_lactoprem_revJPD.xlsx", sheet = 2) |>
  clean_names()

names(otu_raw)
otu <- read_xlsx("DB_lactoprem_revJPD.xlsx", sheet = 2) |>
  clean_names() |>
  rename(otu_id = index) |>         # cambia aquí según el nombre correcto
  mutate(id_match = toupper(otu_id))


# Reestructurar funcional
# 1. Reestructurar datos funcionales correctamente
fun_long <- fun |>
  rename(otu_id = `OTU ID`, time = Time, tratamiento = Treatment) |>
  pivot_longer(
    cols = -c(otu_id, time, tratamiento),
    names_to = "pathway_id",
    values_to = "abund"
  ) |>
  mutate(
    # Extraer número de ID desde otu_id: "LAC_100_T0_H" → 100
    id = as.numeric(str_extract(otu_id, "(?<=LAC_)[0-9]+"))
  )

# 2. Unir con datos clínicos (por ID numérico)
tabla <- fun_long |>
  left_join(clin, by = "id")


# Reemplazo manual del log-CLR porque microbiome no está disponible
clr <- function(x) {
  log(x / exp(mean(log(x[x > 0]))))  # solo valores > 0
}

tabla <- tabla |>
  group_by(id) |>                                # cambiar id_match → id
  mutate(rel = abund / sum(abund, na.rm = TRUE)) |>
  ungroup() |>
  mutate(clr = clr(rel))

# Modelos mixtos
tabla <- tabla |> 
  rename(tratamiento = tratamiento.x)

table(tabla$tratamiento)

res_lmer <- tabla |>
  group_by(pathway_id) |>
  do({
    fit <- lmer(clr ~ tratamiento * time + (1 | id), data = .)
    tidy(fit, effects = "fixed")
  }) |>
  ungroup() |>
  filter(term %in% c("tratamiento", "time", "tratamiento:time")) |>
  mutate(q = p.adjust(p.value, "BH"))




