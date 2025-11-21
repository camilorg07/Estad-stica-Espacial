```{r}
#| context: server 
#| output: false

proj4string(puntajes) <- CRS("EPSG:3310")


g_CP1<-gstat(id = "f1", form = f1 ~ 1,
             data = puntajes, model = exp_pc1T)
g_CP2<-gstat(id = "f2", form = f2 ~ 1,
             data = puntajes, model = exp_pc2T)
g_CP3<-gstat(id = "f3", form = f3 ~ 1,
             data = puntajes, model = exp_pc3T)

#Kriging ------------------------------------------------

invisible(predic_CP1 <- predict(g_CP1, newdata = new, nmin = 20))
invisible(predic_CP2 <- predict(g_CP2, newdata = new, nmin = 20))
invisible(predic_CP3 <- predict(g_CP3, newdata = new, nmin = 20))

pred_pca <- cbind(
  predic_CP1@data["f1.pred"], 
  predic_CP2@data["f2.pred"], 
  predic_CP3@data["f3.pred"]
)

names(pred_pca) <- c("PC1.pred", "PC2.pred", "PC3.pred")

coef_predichos_variacion <- as.matrix(V_matrix) %*% t(pred_pca) 

# Paso 2: Coeficientes totales (k filas x N_new columnas)
# Sumar el vector de la media a CADA columna de la matriz de variación.
coef_predichos_totales <- sweep(x = coef_predichos_variacion, # Matriz de variación (k x N_new)
                                MARGIN = 1,                     # Sumar por FILAS (cada fila es un coeficiente de base)
                                STATS = mean_coefs,             # Vector de la media (k x 1)
                                FUN = "+")

mean_coefs <- PCA$meanfd$coefs
mean_coefs <- as.vector(mean_coefs)

M<-NULL
MM<-NULL
for (i in 1:24) {
  M<-coef_predichos[i,]+mean_coefs
  MM<-cbind(MM,M)
}

fd_temp_pred <- fd(coef = coef_predichos_totales, basisobj = PCA$harmonics$basis)

plot(fd_temp_pred)

matplot(temp_pred_eval, type="l")

matplot(fd_temp_pred$coefs, type="l")


t_points <- seq(0, 168, by=1)
temp_pred_eval <- eval.fd(0:168, fd_temp_pred)

time_h <- 0:168

inicio <- as.POSIXct("2024-03-01 00:00:00", tz = "America/Bogota")
inicio <- inicio + (0:744)*3600

predic_CP1@coords

temp_df <- data.frame(time = inicio, temp_pred_eval)

temp_df_clean <- as.data.frame(lapply(temp_df, function(x) as.numeric(x)))
temp_df_clean$time <- temp_df$time


p <- ggplot() +
  labs(x = "Fecha", y = "Temperatura") +
  theme_light()

for (i in 2:ncol(temp_df)) {
  
  df_i <- data.frame(
    time  = temp_df$time,
    value = temp_df[[i]]
  )
  
  p <- p + geom_line(
    data = df_i,
    aes(x = time, y = value),
    color = cols[i - 1],
    alpha = 0.8
  )
}

# 1. Convert to Long Format
# Assuming 'time' is your only ID variable and other columns are temperatures
temp_df_long <- temp_df %>%
  pivot_longer(
    cols = -time, # Pivot all columns except 'time'
    names_to = "Series_Name", # New column for the old column names (e.g., 'temp_A', 'temp_B')
    values_to = "Temperature" # New column for the temperature values
  )


rownames(temp_pred_eval)<-inicio

matplot(temp_pred_eval)

# 2. Plot in One Call
p <- ggplot(
  data = temp_df_long,
  aes(
    x = time,
    y = Temperature,
    group = Series_Name, # Group lines by the original column name
    color = Series_Name # Color based on the original column name
  )
) +
  geom_line(alpha = 0.8) +
  labs(x = "Fecha", y = "Temperatura") +
  theme_light()
temp_df$time <- as.POSIXct(temp_df$time,
                            format = "%m/%d/%y %H:%M:%S",
                            tz = "America/Bogota")


p
temp_long <- temp_df %>%
  pivot_longer(
    cols = -time,
    names_to = "serie",
    values_to = "valor"
  )


ggplot(temp_long, aes(x = time, y = valor, color = serie)) +
  geom_line(alpha = 0.8) +
  labs(x = "Fecha", y = "Temperatura") +
  theme_light()

ggplot()+
  geom_line(aes(x=temp_df$time, y = temp_df[,2]))

cols <- rainbow(ncol(temp_df)-1)

p <- ggplot() +
  labs(x = "Fecha", y = "Temperatura") +
  theme_light()

for (i in 2:ncol(temp_df)) {
  
  df_i <- data.frame(
    time = temp_df$time,
    value = temp_df[[i]]
  )
  
  p <- p + geom_line(
    data = df_i,
    aes(x = time, y = value),
    color = cols[i-1],
    alpha = 0.8
  )
}




temp_long <- temp_df %>%
  pivot_longer(
    cols = -time,
    names_to = "curve",
    values_to = "value"
  )

# 4. Graficar
ggplot(temp_long, aes(x = time, y = value, color = curve)) +
  geom_line(alpha = 0.9) +
  labs(x = "Fecha", y = "Temperatura", color = "Curva") +
  theme_light()



temp_long <- temp_long %>%
  mutate(
    DateTime = inicio + time * 3600   # cada fila es una hora
  )



matplot(temp_pred_eval, type)

# 1. Inicializar predic_spdf con los resultados de Kriging
# (predic_CP1 ya contiene las coordenadas y CRS de la grilla)
predic_spdf <- as(predic_CP1, "SpatialPointsDataFrame")

# 2. Agregar las 24 columnas de temperatura AL predic_spdf
hora_names <- paste0("Temp_h", 0:23)
predic_spdf@data[, hora_names] <- temp_pred_eval

# 3. CONVERTIR A SpatialPixelsDataFrame AHORA, después de tener TODAS las columnas
predic_spixdf <- as(predic_spdf, "SpatialPixelsDataFrame") # ¡CORRECCIÓN!

raster_list <- list()
CRS_ORIGEN <- "+init=epsg:3310" 
CRS_WGS84 <- "+init=epsg:4326"

for (i in 0:23) {
  pred_layer_name <- hora_names[i + 1] 
  
  # r_temp_utm ahora encontrará la columna en predic_spixdf
  r_temp_utm <- raster(predic_spixdf, layer = pred_layer_name) 
  crs(r_temp_utm) <- CRS_ORIGEN
  
  # ... (resto del código de proyección)
  r_temp_wgs84 <- raster::projectRaster(
    raster::mask(r_temp_utm, sh_mundos_sp_utm_simple),
    crs = CRS_WGS84,  
    method = "bilinear"
  )
  
  raster_list[[i + 1]] <- r_temp_wgs84
}

temp_raster_stack <- stack(raster_list)

mins_por_capa <- cellStats(temp_raster_stack, stat='min')
maxs_por_capa <- cellStats(temp_raster_stack, stat='max')

# 2. Calcular el valor MÍNIMO GLOBAL (el valor más bajo de todo el vector)
min_global <- min(mins_por_capa) 
# O más directo: min_global <- min(cellStats(temp_raster_stack, stat='min'))

# 3. Calcular el valor MÁXIMO GLOBAL
max_global <- max(maxs_por_capa)
# O más directo: max_global <- max(cellStats(temp_raster_stack, stat='max'))

# 4. Recalcular los bins con los nuevos valores de longitud 1
pred_range_global <- c(min_global, max_global)
step_pred_global <- (max_global - min_global) / 6
bins_global <- round(seq(min_global, max_global, by = step_pred_global), 1)

step_pred_global <- (max_global - min_global) / 6
bins_global <- round(seq(min_global, max_global, by = step_pred_global), 1)

pal_pred_global <- colorBin(
  palette = "viridis",
  domain = pred_range_global,
  bins = bins_global,
  na.color = "transparent"
)

ui <- fluidPage(
  # Título
  h2("🌡️ Predicción de Temperatura por Hora (Kriging Funcional)"),
  hr(),
  
  # Layout con sidebar para el slider
  sidebarLayout(
    sidebarPanel(
      # Slider para seleccionar la hora (0 a 23)
      sliderInput("hora_seleccionada", 
                  "Seleccionar Hora del Día (0 a 23):", 
                  min = 0, 
                  max = 23, 
                  value = 12, # Valor inicial al mediodía
                  step = 1,
                  animate = animationOptions(interval = 1000) # Opción de play automático
      ),
      # Display de la hora seleccionada
      textOutput("hora_display")
    ),
    
    # Panel principal para el mapa
    mainPanel(
      leafletOutput("temp_map", height = 600)
    )
  )
)

ui <- fluidPage(
  # Título
  h2("🌡️ Predicción de Temperatura por Hora (Kriging Funcional)"),
  hr(),
  
  # Layout con sidebar para el slider
  sidebarLayout(
    sidebarPanel(
      # Slider para seleccionar la hora (0 a 23)
      sliderInput("hora_seleccionada", 
                  "Seleccionar Hora del Día (0 a 23):", 
                  min = 0, 
                  max = 23, 
                  value = 12, # Valor inicial al mediodía
                  step = 1,
                  animate = animationOptions(interval = 1000) # Opción de play automático
      ),
      # Display de la hora seleccionada
      textOutput("hora_display")
    ),
    
    # Panel principal para el mapa
    mainPanel(
      leafletOutput("temp_map", height = 600)
    )
  )
)

# ===================================================================
# PASO 3: LÓGICA DEL SERVIDOR (SERVER)
# ===================================================================

server <- function(input, output, session) {
  
  # 1. Renderizar el Mapa Base (Solo una vez)
  output$temp_map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      addLegend(
        pal = pal_pred_global,  
        values = pred_range_global,
        title = "Temperatura Predicha (°C)",
        position = "topright"
      )
  })
  
  # 2. Observador para actualizar el Raster cuando el slider cambia
  observe({
    # Asegura que el valor del slider esté disponible
    req(input$hora_seleccionada) 
    
    hora_sel <- input$hora_seleccionada # Hora seleccionada (0 a 23)
    raster_index <- hora_sel + 1        # Índice del stack (1 a 24)
    
    # Obtiene el raster de la hora seleccionada
    raster_sel <- temp_raster_stack[[raster_index]] 
    
    # Actualiza el texto de la hora
    output$hora_display <- renderText({
      paste("Mostrando predicción para la Hora:", hora_sel)
    })
    
    # Actualiza la capa de raster en el mapa
    leafletProxy("temp_map", data = raster_sel) %>%
      clearImages() %>% # Limpia cualquier capa de imagen (raster) anterior
      addRasterImage(
        raster_sel,
        colors = pal_pred_global,
        opacity = 0.7,
        project = FALSE # Es importante si el raster ya está en WGS84
      )
  })
}

shinyApp(ui = ui, server = server)
```