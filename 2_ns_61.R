################ MALDI-TOF ANALISIS OTAMENDI ###################################
################ 1) NO SUPERVISADO - 61 INDIVIDUOS  ############################
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################

library("readBrukerFlexData")
library("binda")
library("fs")
library("readxl")
library("MALDIquant")
library("MALDIquantForeign")
library("MALDIrppa")
library("tidyverse")
library("dplyr")
library("clValid")
library(cluster)
library(factoextra)
library(ggplot2)
library(gridExtra)
library(ggrepel)

### CARGA DE ARCHIVOS ##########################################################
################################################################################

# Detectar el usuario actual del sistema
usuario <- Sys.info()[["user"]]

# Definir la ruta del proyecto de acuerdo al usuario
ruta_proyecto <- ifelse(usuario == "urtea", 
                        "C:/Users/urtea/OneDrive/Documents/Proyectos/MALDI_Otamendi_plasma",
                        "C:/Users/Facundo/Desktop/MALDI_Otamendi_plasma")

# Crear la ruta de datos usando file.path
ruta_datos <- file.path(ruta_proyecto)

# Cargar los archivos Rdata con la ruta definida
archivos_a_cargar <- c("matint_61_ind_dico.Rdata", "matint_61_ind.Rdata")

# Cargar cada archivo en la ruta de datos
for (archivo in archivos_a_cargar) {
  load(file.path(ruta_datos, archivo))
}

df_metadata <- df_metadata_prom_mue

### SELECCIÓN DE PICOS #########################################################
################################################################################

# Selección de picos para binary discriminant analysis (BDA)
factor_tipo <- factor(df_metadata_prom_mue$estado)
is.binaryMatrix(matint_61_ind_dico) # TRUE
br <- binda.ranking(matint_61_ind_dico, factor_tipo, verbose = FALSE)

# Gráfico de picos vs score 
nueva_columna <- c()
matriz <- matrix(br, nrow = 496, ncol = 4) #496 es la cantidad de picos
for (i in 1:496) {
  nuevo_valor <- colnames(matint_61_ind_dico)[br[i]]
  nueva_columna<- c(nueva_columna, nuevo_valor)
}
matriz <- cbind(matriz, nueva_columna)
df_br <- data.frame(matriz)
plot(df_br$nueva_columna, df_br$V2, 
     xlab = "m/z", ylab = "Score", 
     main = "Ranking de picos de los espectros (61 ind-dia)")
# Crear un gradiente de colores (por ejemplo, de azul a rojo)
colores <- colorRampPalette(c("green4", "red2"))(496)
# Agregar puntos con colores en forma de gradiente
for (i in 1:496) {
  points(df_br$nueva_columna[i], df_br$V2[i], col = colores[i]) 
}
# Agregar puntos con relleno de colores en forma de gradiente
for (i in 1:496) {
  points(df_br$nueva_columna[i], df_br$V2[i], pch = 19, col = colores[i]) 
}

# Selección de picos mas preponderantes
top.b5 <- br[1:5]  ## primeros 5 picos
top.b10 <- br[1:10]  ## primeros 10 picos
top.b15 <- br[1:15]  ## primeros 15 picos
top.b20 <- br[1:20]  ## primeros 20 picos
top.b30 <- br[1:30]  ## primeros 30 picos 
top_actual <- top.b10

# Elección de mejores algoritmos de clustering
comparacion <- clValid(
  obj        = matint_61_ind_dico[, top_actual],
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("stability", "internal")
)
summary(comparacion)
optimalScores(comparacion) #Se puede ir probando con distintos top picos


### ALGORITMO DE CLUSTERING ####################################################
################################################################################

# PAM clustering con top5 y 2 clusters

top_actual <- br[setdiff(1:15, 4:9)]
K.num <- 2 # clusters
var2 = 0.95

pam.top8 <- pam(matint_61_ind[, top_actual], 
                K.num, metric = "manhattan")

# Recalcular el clustering y el gráfico con elipse
cluster_pam8 <- fviz_cluster(pam.top8, 
                             ellipse.type = "convex",  # Cambiar a "norm" o "t" para ver la elipse
                             data = matint_61_ind[, top_actual],
                             outlier.color = "black",
                             ggtheme = theme_grey(),
                             ellipse.level = var2, # Nivel de confianza de la elipse
                             show.clust.cent = FALSE,  # Mostrar el centroide de los clusters
                             geom = "point",
                             pointsize = 0.5,
                             main = "CTRL/SEPSIS - ind/dia - PAM - 9 picos - 2 cl") +
  scale_colour_manual(values = c("gray", "gray")) +
  scale_fill_manual(values = c("gray", "gray")) 

cluster_pam8 <- cluster_pam8 + 
  geom_point(data = cluster_pam8$data, 
             aes(x = x, y = y, color = df_metadata$grupo, size = as.factor(df_metadata$dia_abs))) +
  # geom_line(data = cluster_pam10$data, 
  #           aes(x = x, y = y, group = df_metadata$paciente), 
  #           color = "black", size = 0.5, alpha = 0.6) +  # Agrega líneas conectando los puntos de cada paciente
  scale_color_manual(values = c("gray","gray","olivedrab3","khaki4","steelblue2","mediumpurple2")) +
  scale_size_manual(values = c("0" = 1, "2" = 2.5, "4" = 4, "6" = 5.5)) +
  labs(color = "Grupo", size = "Día abs") +
  theme(legend.position = "right")

# Mostrar el gráfico actualizado con elipse
print(cluster_pam8)

# PCA para el biplot

matint_61_top <- matint_61_ind[, top_actual]

# Realizar PCA para obtener los componentes principales
pca_res <- prcomp(matint_61_top, center = TRUE, scale. = TRUE)

# Obtener las coordenadas de los vectores (las variables)
var_coords <- pca_res$rotation[, 1:2]  # Coordenadas de las variables en los dos primeros PCs
var_names <- rownames(var_coords)  # Nombres de las variables

# Escalar las coordenadas de los vectores (para reducir la longitud de los vectores)
scale_factor <- 0.6  # Ajustar este factor para reducir más o menos el tamaño de los vectores
var_coords_scaled <- var_coords * scale_factor

# Agregar los vectores y sus etiquetas ajustadas al gráfico de clustering
cluster_pam8 <- cluster_pam8 +
  geom_segment(data = as.data.frame(var_coords_scaled),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = as.data.frame(var_coords_scaled),
                  aes(x = PC1, y = PC2, label = var_names),
                  color = "black",
                  size = 3,          # Tamaño de las etiquetas de los vectores
                  box.padding = 0.35, # Aumentar separación entre etiquetas
                  point.padding = 0.3,
                  segment.size = 0.2)

# Mostrar el gráfico con los vectores ajustados
print(cluster_pam8)




### GRAFICA DE ESPECTROS Y LUGARES DE INTERÉS ##################################
################################################################################


# Cargar librería necesaria
library(scales)  # Para agregar transparencia a los colores

# Definir los puntos de referencia y tolerancia
highlight_positions <- as.integer(round(as.numeric(colnames(matint_61_top))))
tolerance <- 10

# Definir los límites del eje x
x_lim <- c(2000, 8000)

# Definir los colores con transparencia (50% de opacidad)
colors <- c(alpha("blue", 0.3), alpha("blue", 0.3), alpha("blue", 0.3), alpha("blue", 0.3), alpha("blue", 0.3),
            alpha("blue", 0.3), alpha("blue", 0.3), alpha("blue", 0.3), alpha("blue", 0.3))

# Espectros que quiero graficar
selected_spectra <- c(1, 6, 20, 50)

# Ajustar la cantidad de espectros por gráfico
spectra_per_plot <- 2
total_spectra <- length(selected_spectra)

# Calcular el número de gráficos necesarios
num_plots <- ceiling(total_spectra / spectra_per_plot)

# Loop para generar los gráficos
for (plot_idx in 1:num_plots) {
  
  # Definir los índices de los espectros que van en este gráfico
  start_idx <- (plot_idx - 1) * spectra_per_plot + 1
  end_idx <- min(plot_idx * spectra_per_plot, total_spectra)
  
  # Configurar un layout de 2 filas y 1 columna (2 subplots por gráfico)
  par(mfrow = c(spectra_per_plot, 1), mar = c(4, 4, 2, 2))  # Márgenes ajustados
  
  # Graficar cada espectro en el rango de este gráfico
  for (i in start_idx:end_idx) {
    spec_idx <- selected_spectra[i]  # Obtener el índice correcto
    mass <- Spectra_list_prom_muestra[[spec_idx]]@mass
    intensity <- Spectra_list_prom_muestra[[spec_idx]]@intensity
    
    # Graficar el espectro individual
    plot(mass, intensity, type = "l", col = "black",
         xlab = "Mass/Charge (m/z)", ylab = "Intensity",
         main = paste("Espectro ", spec_idx,"- grupo: ", df_metadata[spec_idx,7]),
         xlim = x_lim)
    
    # Dibujar las barras semitransparentes con tolerancia
    for (j in 1:length(highlight_positions)) {
      rect(highlight_positions[j] - tolerance, par("usr")[3],
           highlight_positions[j] + tolerance, par("usr")[4],
           col = colors[j], border = NA)
    }
    
    # Dibujar el espectro encima para que la barra quede de fondo
    lines(mass, intensity, col = "black")
  }
}

### GRAFICA DE TASA DE ACIERTO Y CÁLCULO DE MÉTRICAS ###########################
################################################################################

print("Resultados de PAM.top8.k2: ")

# Obtener los nombres de las muestras y los clusters asignados
sample_names <- rownames(matint_61_ind)
cluster_assignments <- pam.top8$cluster

# Contar la cantidad total de muestras por categoría
total_CTL <- sum(df_metadata$estado == "CTL")
total_SEPSIS <- sum(df_metadata$estado == "SEP")

# Contar la cantidad de aciertos en la clasificación
acierto_CTL <- sum(df_metadata$estado == "CTL" & cluster_assignments == 1)
error_CTL <- sum(df_metadata$estado == "CTL" & cluster_assignments == 2)
acierto_SEPSIS <- sum(df_metadata$estado == "SEP" & cluster_assignments == 2)
error_SEPSIS <- sum(df_metadata$estado == "SEP" & cluster_assignments == 1)

# Totales por cluster
total_cluster1 <- sum(cluster_assignments == 1)
total_cluster2 <- sum(cluster_assignments == 2)

# Cálculo de la tasa de acierto
tasa_acierto_CTL <- acierto_CTL / total_CTL
tasa_acierto_SEPSIS <- acierto_SEPSIS / total_SEPSIS

# Crear un dataframe con los resultados para el gráfico
datos_grafico <- data.frame(
  Categoria = c("CTL", "SEPSIS"),
  Aciertos = c(acierto_CTL, acierto_SEPSIS),
  Total = c(total_CTL, total_SEPSIS),
  Tasa_Acierto = c(tasa_acierto_CTL, tasa_acierto_SEPSIS)
)

# Gráfico de tasas de acierto
grafico_acierto <- ggplot(datos_grafico, aes(x = Categoria, y = Tasa_Acierto, fill = Categoria)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Tasa_Acierto)), vjust = -0.5) +
  labs(title = "Tasa de Acierto - PAM Clustering 2 clusters",
       x = "Categoría",
       y = "Tasa de Acierto") +
  scale_fill_manual(values = c("CTL" = "khaki4", "SEPSIS" = "mediumpurple4")) +
  theme_minimal()

plot(grafico_acierto)

# 1. Silhouette Score
silhouette_score <- silhouette(cluster_assignments, dist(matint_61_ind[, top_actual]))
fviz_silhouette(silhouette_score)

# 2. Within-cluster Sum of Squares (WCSS)
# Obtener los medoides de cada cluster
medoids <- pam.top8$medoids  

# Obtener los clusters asignados
clusters <- pam.top8$clustering  

# Calcular WCSS manualmente sumando distancias dentro de cada cluster
wcss <- sum(pam.top8$clusinfo[, "av_diss"] * pam.top8$clusinfo[, "size"])
print(wcss)



# 3. Between-cluster Sum of Squares (BCSS)

# Calcular el centroide general (promedio de los datos)
centroid_global <- colMeans(matint_61_ind[, top_actual])

# Calcular la distancia de cada medoide al centroide global
bcss <- sum(rowSums((medoids - centroid_global)^2) * pam.top8$clusinfo[, "size"])
print(bcss)


# Crear la matriz de confusión
conf_matrix <- table(Real = df_metadata$estado, Predicho = cluster_assignments)

# Extraer valores de la matriz de confusión
TP <- conf_matrix["SEP", "2"]  # Verdaderos positivos (SEPSIS correctamente clasificados)
FN <- conf_matrix["SEP", "1"]  # Falsos negativos (SEPSIS clasificados como CTL)
FP <- conf_matrix["CTL", "2"]  # Falsos positivos (CTL clasificados como SEPSIS)
TN <- conf_matrix["CTL", "1"]  # Verdaderos negativos (CTL correctamente clasificados)

# Calcular precisión, sensibilidad (recall) y F1-score
precision <- TP / (TP + FP)
sensibilidad <- TP / (TP + FN)
f1_score <- 2 * (precision * sensibilidad) / (precision + sensibilidad)

# Imprimir los resultados
cat("Precisión:", round(precision, 4), "\n")
cat("Sensibilidad:", round(sensibilidad, 4), "\n")
cat("F1-score:", round(f1_score, 4), "\n")