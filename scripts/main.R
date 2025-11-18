# =========================================================
#
#
#
#
#
#
#
#
#
#
# =========================================================

# =========================================================
# Part 0: Cargar datos, buenas prácticas y librerías
# =========================================================

# setwd("~/Desktop/Taller 1 - BigData/Taller4_PIB_Regional_italiano")
# setwd("C:/Users/Asuar/OneDrive/Escritorio/Libros Clases/Economía/Ciencia Datos y Econometria/Taller4_PIB_Regional_italiano")

rm(list = ls())

require(pacman)
p_load(sf, tidyverse, dplyr, knitr, stargazer, ggplot2, stringr, 
       spdep, spatialreg)

# Cargar datos espaciales
italia_espaciales <- st_read("stores/Reg2014_ED50g/Reg2014_ED50_g.shp")
# Cargar datos csv
GDP_ita <- read_csv("stores/Data.csv")


# =========================================================
# Part 1: Análisis descriptivo y territorial
# =========================================================

# Estadísticas descriptvas de regiones de Italia
vars <- c("GDP", "K", "L")

#  Tabla de regiones ordenadas por GDP
tabla_regiones <- GDP_ita[order(GDP_ita$GDP, decreasing = TRUE), c("Territory", vars)]
print(tabla_regiones)

# Estadísticas descriptivas generales de las variables
tabla_summary <- do.call(rbind, lapply(GDP_ita[vars], summary))

# Mostrar con nombres más limpios
tabla_summary <- as.data.frame(tabla_summary)
tabla_summary


# Mapas temáticos usando polígonos regionales

# Arreglar problema nombres previo a join
italia_norm <- italia_espaciales %>%
  mutate(
    REGIONE_raw = REGIONE,
    REGIONE = REGIONE %>%
      str_replace_all("[-–—]", " ")
  )

# Unir shapefile con datos económicos
italia_map <- italia_norm %>%
  left_join(GDP_ita, by = c("REGIONE" = "Territory"))


# --------------------------------------------------------
# Función auxiliar para mapas
plot_map <- function(data, variable, titulo) {
  ggplot(data) +
    geom_sf(aes(fill = !!sym(variable)), color = "white", size = 0.2) +
    scale_fill_viridis_c(option = "plasma", direction = -1) +
    labs(title = titulo, fill = variable) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
}

# Mapas por variables: GDP, K, L
map_GDP <- plot_map(italia_map, "GDP", "Producto Interno Bruto (GDP) por región")
map_K   <- plot_map(italia_map, "K", "Capital (K) por región")
map_L   <- plot_map(italia_map, "L", "Trabajo (L) por región")

# Mostrar en consola
map_GDP
map_K
map_L

# --------------------------------------------------------
# Exploración de correlaciones simples entre variables

pairs(GDP_ita[vars],
      main = "Matriz de dispersión entre variables económicas",
      pch = 19, col = "blue")


# =========================================================
# Part 2: Modelo base (MCO)
# =========================================================

# Estimación Cobb-Douglas en forma log-lineal
model_CB <- lm( log(GDP)~ log(K) + log(L), data=GDP_ita)
summary(model_CB)

# Suma de elasticidades
elasticidad_suma <- sum(coef(model_CB)[-1])
elasticidad_suma


# =========================================================
# Part 3: Analisis de dependencia espacial
# =========================================================

# Verificar el sistema de coordenadas
st_crs(italia_map)
## De acuerdo a este codigo, el archivo ya se encuentra en metros

# --- 3.1 Construcción de matriz de pesos espaciales (W) ---
# Obtener centroides de las regiones
coords <- st_centroid(st_geometry(italia_map))

# Probando diferentes criteros para encontrar el punto optimo donde cada punto tenga al
# menos un vecino 50, 100, 200, 300, 350, 370, 379

# CRITERIO 1: Matriz por DISTANCIA (umbral escogido = 379 km)
# Nota: Ajusta el umbral según la escala de tus datos
dist_380km <- spdep::dnearneigh(coords, 0, 379000) # 379,000 metros

# Verificar regiones sin vecinos
n_sin_vecinos_dist <- sum(card(dist_380km) == 0)
print(paste("=== MATRIZ DE PESOS POR DISTANCIA ==="))
print(paste("Umbral: 379 km"))
print(paste("Regiones sin vecinos:", n_sin_vecinos_dist))

# Si hay regiones sin vecinos, mostrar cuáles son
if(n_sin_vecinos_dist > 0) {
  regiones_sin_vecinos <- italia_map$REGIONE[which(card(dist_380km) == 0)]
  print(paste("Regiones sin vecinos:", paste(regiones_sin_vecinos, collapse = ", ")))
}

# Crear matriz de pesos W (row-standardized)
W_dist <- spdep::nb2listw(dist_380km, style = "W", zero.policy = TRUE)

# --- 3.2 Dispersión y Número promedio de vecinos
summary(dist_380km)
# Numero promedio de vecinos es 8.2
# Sparsity seria 41%

# --- 3.3 Test de Moran´s
# Primero con residuos de MCO
residuos_mco <- residuals(model_CB)
moran_dist <- spdep::moran.test(residuos_mco, W_dist, 
                                alternative = "two.sided",
                                zero.policy = TRUE)
print(moran_dist)
# P-Value = 5.17e-09 - Significativo, Si hay autocorrelación espacial, El modelo MCO **NO capturó** toda la estructura espacial
# Hay efectos espaciales (spillovers, externalidades) que el modelo ignora.

# Segundo con la variable GDP
moran_gdp <- spdep::moran.test(italia_map$GDP, W_dist, 
                               alternative = "two.sided",
                               zero.policy = TRUE)
print(moran_gdp)
# P-Value = 0.4592 - No es significativo, No hay autocorrelación espacial


# =========================================================
# Part 4: Modelos Espaciales
# =========================================================

# --- 4.1 Modelo Spatial Lag (SAR) ---
model_SAR <- lagsarlm(log(GDP) ~ log(K) + log(L), 
                      data = italia_map, 
                      listw = W_dist)
summary(model_SAR)

# --- 4.2 Modelo SARAR (SAC) ---
model_SAC <- sacsarlm(log(GDP) ~ log(K) + log(L), 
                      data = italia_map, 
                      listw = W_dist)
summary(model_SAC)

# --- 4.3 Comparación de modelos ---
# AIC y BIC para comparar
AIC(model_CB, model_SAR, model_SAC)
BIC(model_CB, model_SAR, model_SAC)

lm.LMtests(model_CB, W_dist, test = "all")

# Interpretación: 
# - La dependencia espacial está PRINCIPALMENTE en los ERRORES (error espacial)
# - También hay algo de dependencia en lag, pero es secundaria
# - SARMA significativo -> usar modelo SAC
# Modelo recomendado: SAC (SARAR)
# Razón: Captura tanto el error espacial (dominante) como algo de lag espacial


# =========================================================
# Part 5: Análisis de Impactos y Diagnóstico Territorial
# =========================================================

# --- 5.1 Impactos Directos, Indirectos y Totales ---
# El modelo SAC proporciona impactos a través de la matriz de multiplicadores espaciales
# Para SAC: Impactos = (I - ρW)^-1 * [I, β_K*I, β_L*I]

# Cálculo de impactos espaciales con bootstrap
imp_SAC <- impacts(model_SAC, listw = W_dist, R = 2000)

# Resumen
summary(imp_SAC)

# ---- Impactos medios ----
d_mean   <- imp_SAC$res$direct[,1]
ind_mean <- imp_SAC$res$indirect[,1]
tot_mean <- imp_SAC$res$total[,1]

# ---- Desviaciones estándar ----
d_sd   <- apply(imp_SAC$sres$direct,   2, sd)
ind_sd <- apply(imp_SAC$sres$indirect, 2, sd)
tot_sd <- apply(imp_SAC$sres$total,    2, sd)

# ---- Variables ----
vars <- c("log(K)", "log(L)")

# ---- Construir tabla con "media (sd)" ----
impactos_tabla <- data.frame(
  Variable   = vars,
  Directo    = sprintf("%.4f (%.4f)", d_mean,   d_sd),
  Indirecto  = sprintf("%.4f (%.4f)", ind_mean, ind_sd),
  Total      = sprintf("%.4f (%.4f)", tot_mean, tot_sd)
)

# ---- Mostrar tabla ----
knitr::kable(
  impactos_tabla,
  caption = "Impactos Directos, Indirectos y Totales con Desviaciones Estándar – Modelo SAC"
)


# --- 5.2 Identificación de Regiones con Residuos Atípicos ---

# Residuos del modelo
italia_map$Residuos_SAC <- residuals(model_SAC)

# Residuos estandarizados
italia_map$Residuos_std <- scale(italia_map$Residuos_SAC)[,1]

# Regiones con peor desempeño (muy negativo)
head(
  italia_map[order(italia_map$Residuos_std), 
               c("REGIONE", "Residuos_SAC","Residuos_std")],
  5
)

# Regiones con mejor desempeño (muy positivo)
head(
  italia_map[order(-italia_map$Residuos_std), 
             c("REGIONE", "Residuos_SAC","Residuos_std")],
  5
)


# --- 5.3 Visualización: Mapa de Residuos ---

ciudades_italia <- data.frame(
  ciudad = c(
    # Grandes ciudades
    "Roma","Milano","Napoli","Torino","Palermo","Genova","Bologna","Firenze","Bari","Catania",
    # Capitales regionales
    "L'Aquila","Catanzaro","Potenza","Ancona","Trieste","Campobasso","Trento","Bolzano",
    "Aosta","Perugia","Venezia", "Cagliari"
  ),
  lon = c(
    12.4964, 9.1900, 14.2681, 7.6869, 13.3613, 8.9463, 11.3426, 11.2558, 16.8719, 15.0873,
    13.3995, 16.6050, 15.8050, 13.5150, 13.7768, 14.6620, 11.1211, 11.3548,
    7.3151, 12.3889, 12.3155, 9.1100
  ),
  lat = c(
    41.9028, 45.4642, 40.8518, 45.0703, 38.1157, 44.4056, 44.4949, 43.7696, 41.1171, 37.5022,
    42.3499, 38.8890, 40.6395, 43.6168, 45.6495, 41.5620, 46.0707, 46.4983,
    45.7370, 43.1122, 45.4408, 39.2170
  )
)

# Converir a sf para poder usar
ciudades_sf <- st_as_sf(
  ciudades_italia,
  coords = c("lon","lat"),
  crs = 4326
)


# Mapa de residuos
mapa_residuos <- ggplot(italia_map) +
  geom_sf(aes(fill = Residuos_SAC), color = "white", size = 0.2) +
  geom_sf(data = ciudades_sf, color = "black", size = 2) +
  geom_sf_label(
    data = ciudades_sf,
    aes(label = ciudad),
    size = 3,
    label.size = 0.2,
    alpha = 0.7
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    title = "Residuos del Modelo SAC por Región",
    subtitle = "Valores positivos: el modelo subestima el PIB. Valores negativos: lo sobreestima.",
    fill = "Residuos"
  ) +
  theme_minimal()

mapa_residuos


# Mapa de residuos estandarizados
mapa_residuos_std <- ggplot(italia_map) +
  geom_sf(aes(fill = Residuos_std), color = "white", size = 0.2) +
  geom_sf(data = ciudades_sf, color = "black", size = 2) +
  geom_sf_label(
    data = ciudades_sf,
    aes(label = ciudad),
    size = 3,
    label.size = 0.2,
    alpha = 0.7
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    title = "Residuos Estandarizados (z-score)",
    subtitle = "Valores > 2 indican potenciales outliers espaciales",
    fill = "Residuos Std"
  ) +
  theme_minimal()

mapa_residuos_std

