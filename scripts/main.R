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

p_load(sf, tidyverse, dplyr, knitr, stargazer, ggplot2, stringr, spdep)

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
n_sin_vecinos_dist <- sum(card(dist_50km) == 0)
print(paste("=== MATRIZ DE PESOS POR DISTANCIA ==="))
print(paste("Umbral: 50 km"))
print(paste("Regiones sin vecinos:", n_sin_vecinos_dist))

# Si hay regiones sin vecinos, mostrar cuáles son
if(n_sin_vecinos_dist > 0) {
  regiones_sin_vecinos <- italia_map$REGIONE[which(card(dist_50km) == 0)]
  print(paste("Regiones sin vecinos:", paste(regiones_sin_vecinos, collapse = ", ")))
}

# Crear matriz de pesos W (row-standardized)
W_dist <- spdep::nb2listw(dist_50km, style = "W", zero.policy = TRUE)

# --- 3.2 Dispersión y Número promedio de vecinos
summary(dist_380km)
# Numero promedio de vecinos es 8.2
# Sparcity seria 41%

# --- 3.3 Test de Moran´s
# Primero con residuos de MCO
residuos_mco <- residuals(model_CB)
moran_dist <- spdep::moran.test(residuos_mco, W_dist, 
                                alternative = "two.sided",
                                zero.policy = TRUE)
print(moran_dist)
# P-Value = 5.17e-09 - Significativo 

# Segundo con la variable GDP
moran_gdp <- spdep::moran.test(italia_map$GDP, W_dist, 
                               alternative = "two.sided",
                               zero.policy = TRUE)
print(moran_gdp)
# P-Value = 0.4592 - Significativo al 95%


