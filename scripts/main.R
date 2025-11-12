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

# setwd("C:/Users/Asuar/OneDrive/Escritorio/Libros Clases/Economía/Ciencia Datos y Econometria/Taller4_PIB_Regional_italiano")

rm(list = ls())

require(pacman)

p_load(sf, tidyverse, dplyr, knitr, stargazer, ggplot2, stringr)

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







