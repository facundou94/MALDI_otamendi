################ MALDI-TOF OTAMENDI PLASMA #####################################
################ 1) PREPROCESAMIENTO  ##########################################
################################################################################
#
# Autor: Bioing. Facundo Urteaga (IBB-CONICET)
#
#
### CARGA DE LIBRERIAS #########################################################
################################################################################


library(binda)
library(here)
library(dplyr)
library(readBrukerFlexData)
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(stringr)
library(lubridate)


### CARGA DE ESPECTROS #########################################################
################################################################################


# Detectar el usuario actual del sistema
usuario <- Sys.info()[["user"]]

# Definir la ruta del proyecto de acuerdo al usuario
ruta_proyecto <- ifelse(usuario == "urtea", 
                        "C:/Users/urtea/OneDrive/Documents/Proyectos/MALDI_Otamendi_plasma",
                        "C:/Users/Facundo/Desktop/MALDI_Otamendi_plasma")

# Crear la ruta de datos usando file.path
ruta_datos <- file.path(ruta_proyecto)

# Importar espectros
Spectra_list <- importBrukerFlex(file.path(ruta_datos), verbose=FALSE)


### CREACIÓN DE METADATA #######################################################
################################################################################


# Extraer los nombres de archivo
file_names <- sapply(Spectra_list, function(x) x@metaData[["name"]])

# Función para procesar cada nombre
df <- data.frame(index = 1:length(file_names), name = file_names) %>%
  mutate(
    name = str_remove(name, "\\..*$"), # Remover extensión
    name = str_remove(name, "ok$"),     # Remover "ok" al final si está
    name = ifelse(index == 298, str_replace(name, "PJJ", "OJJ"), name),
    
    parts = str_split(name, "_+|(?<=[A-Za-z])(?=\\d{4,6})|(?<=\\d{4,6})(?=[A-Za-z])", simplify = TRUE),
    num_parts = ncol(parts),
    
    paciente = parts[,1],
    dia = case_when(
      str_detect(parts[,2], "^\\d{6}$") ~ str_c(substr(parts[,2], 1, 2), "-", substr(parts[,2], 3, 4)),
      str_detect(parts[,2], "^\\d{4}$") ~ str_c(substr(parts[,2], 1, 2), "-", substr(parts[,2], 3, 4)),
      index <= 45 ~ "00-00", # Asignar fecha 0 a las primeras 45 filas
      TRUE ~ NA_character_
    ),
    
    replica_b = ifelse(num_parts >= 3, parts[,num_parts-2], NA),
    replica_t = ifelse(num_parts >= 3, parts[,num_parts-1], NA)
  ) %>%
  select(index, name, paciente, dia, replica_b, replica_t)


# Crear una copia del DataFrame original para modificar solo los intervalos problemáticos
df_metadata <- df

# Recorrer solo los índices problemáticos
for (i in 1:nrow(df_metadata)) {
  if (i >= 1 && i <= 45) {  # Intervalo 1-45
    parts <- unlist(strsplit(df_metadata$name[i], "_"))
    df_metadata$paciente[i] <- parts[1]
    df_metadata$dia[i] <- "00-00"
    df_metadata$replica_b[i] <- parts[2]
    df_metadata$replica_t[i] <- parts[3]
  } else if (i >= 46 && i <= 54 || i >= 299 && i <= 307) {  # Intervalos 46-54 y 299-307
    name_parts <- gsub("(\\d{6})([A-Za-z]+)(\\d+)_(\\d+)", "\\1_\\2_\\3_\\4", df_metadata$name[i])
    parts <- unlist(strsplit(name_parts, "_"))
    df_metadata$dia[i] <- paste0(substr(parts[1], 1, 2), "-", substr(parts[1], 3, 4))
    df_metadata$paciente[i] <- parts[2]
    df_metadata$replica_b[i] <- parts[3]
    df_metadata$replica_t[i] <- parts[4]
  }
}

# Definir los patrones
patron_ctl_sano <- "CTL_SANO"
patron_ctl_uti <- "CTL_UTI"
patron_sepsis <- "SEPSIS"
patron_shock <- "SHOCK"

# Crear una nueva columna para el estado
df_metadata$grupo <- NA
df_metadata$estado <- NA

# Recorrer las filas del dataframe
for (i in 1:nrow(df_metadata)) {
  # Obtener el nombre del archivo
  nombre <- Spectra_list[[i]]@metaData$file
  
  # Comprobar los patrones y asignar el estado correspondiente
  if (str_detect(nombre, patron_ctl_sano)) {
    df_metadata$grupo[i] <- "CTL_SANO"
    df_metadata$estado[i] <- "CTL"
  } else if (str_detect(nombre, patron_ctl_uti)) {
    df_metadata$grupo[i] <- "CTL_UTI"
    df_metadata$estado[i] <- "CTL"
  } else if (str_detect(nombre, patron_sepsis)) {
    df_metadata$grupo[i] <- "SEPSIS"
    df_metadata$estado[i] <- "SEP"
  } else if (str_detect(nombre, patron_shock)) {
    df_metadata$grupo[i] <- "SHOCK"
    df_metadata$estado[i] <- "SEP"
  }
  
  # Guardo metadata en el espectro
  Spectra_list[[i]]@metaData$grupo <- df_metadata$grupo[i]
  Spectra_list[[i]]@metaData$paciente <- df_metadata$paciente[i]
  Spectra_list[[i]]@metaData$dia <- df_metadata$dia[i]
  Spectra_list[[i]]@metaData$replica_b <- df_metadata$replica_b[i]
  Spectra_list[[i]]@metaData$replica_t <- df_metadata$replica_t[i]
  
}

# Suponiendo que tu DataFrame se llama df_metadatos
# Corregir las filas 321 a 331 con la fecha correcta
df_metadata$dia[321:331] <- "05-06"

# Convertir la columna "dia" a formato de fecha (usando año 2024)
df_metadata$fecha <- as.Date(paste0("2024-", df_metadata$dia), format = "%Y-%d-%m")

# Ahora, agrupar por paciente y calcular el "día 0"
df_metadata <- df_metadata %>%
  group_by(paciente) %>%
  mutate(
    # Asignar "dia_abs = 0" si la fecha es NA
    dia_abs = case_when(
      is.na(fecha) ~ 0,  # Si la fecha es NA, asignar dia_abs = 0
      TRUE ~ as.numeric(fecha - min(fecha))  # El resto, calcular los días desde el día 0
    )
  ) %>%
  ungroup()

  # Creación de factores de agrupamiento para su uso posterior
  df_metadata$f_pac_dabs <- paste0(df_metadata$paciente, "_", df_metadata$dia_abs)
  # Creación de factores de agrupamiento para su uso posterior
  df_metadata$f_pac_dabs_rep_b <- paste0(df_metadata$paciente, "_", df_metadata$dia_abs, "_", df_metadata$replica_b)


### CONTROL DE CALIDAD Y LIMPIEZA DE ESPECTROS #################################
################################################################################


# Screening inicial: Detección de espectros de baja calidad
sc.results <- screenSpectra(Spectra_list, meta = df_metadata)
summary(sc.results)
plot(sc.results, labels = TRUE)

# plot(Spectra_list[[253]]) # Ploteo de espectros ruidosos (no hay)

# Descartamos espectros defectuosos
Spectra_list_f1 <- sc.results$fspectra # Filtramos espectros
df_metadata_f1 <- sc.results$fmeta # Filtramos metadatos


### FILTRADO Y TRANSFORMACIÓN DE ESPECTROS #####################################
################################################################################

# Parámetros de procesamiento de espectros
thScale <- 2.5 # Smoothing
ite <- 105 # Baseline correction
SigNoi <- 2 # Peak extraction
hws <- 20 # Peak extraction
tol <- 0.03 # Peak binning

# Transformación/filtrado/corrección de espectros con parámetros definidos
# 1) Transformación de intensidad por medio de función sqrt
Spectra_list_f1 <- transfIntensity(Spectra_list_f1, fun = sqrt)
plot(Spectra_list_f1[[30]])
# 2) Suavizado del espectro mediante el método Wavelet
Spectra_list_f1 <- wavSmoothing(Spectra_list_f1, method = "Wavelet", n.levels = 4)
plot(Spectra_list_f1[[30]])
# Detección de la linea de base
baseline <- estimateBaseline(Spectra_list_f1[[30]], method = "SNIP",
                             iterations = ite)
plot(Spectra_list_f1[[30]])
lines(baseline, col="red", lwd=2)
# 3) Remoción de linea de base mediante método SNIP
Spectra_list_f2 <- removeBaseline(Spectra_list_f1, method = "SNIP",
                                  iterations = ite)
plot(Spectra_list_f2[[30]])
# 4) Calibración de intensidad mediante método PQN
Spectra_list_f2 <- calibrateIntensity(Spectra_list_f2, method = "PQN")
plot(Spectra_list_f2[[30]])
# 5) Alineación de espectros
Spectra_list_f3 <- alignSpectra(Spectra_list_f2,
                                halfWindowSize=hws,
                                SNR=SigNoi,
                                tolerance=0.02, # Parámetro sensible
                                warpingMethod="lowess")
plot(Spectra_list_f3[[30]])


### PROMEDIO DE LECTURAS DE UNA MISMA RÉPLICA TÉCNICA ##########################
################################################################################


# Promedio de lecturas de una misma well
Spectra_list_prom_rep <- averageMassSpectra(Spectra_list_f3,
                                            labels = factor(df_metadata_f1$f_pac_dabs_rep_b),
                                            method = "mean")

# Creo la nueva metadata de los espectros promediados
df_metadata_prom_rep <- df_metadata_f1 %>%
  distinct(df_metadata_f1$f_pac_dabs_rep_b, .keep_all = TRUE)

# Promedio de wells de una misma muestra
Spectra_list_prom_muestra <- averageMassSpectra(Spectra_list_prom_rep,
                                                labels = factor(df_metadata_prom_rep$f_pac_dabs),
                                                method = "mean")

# Creo la nueva metadata de los espectros promediados
df_metadata_prom_mue <- df_metadata_prom_rep %>%
  distinct(df_metadata_prom_rep$f_pac_dabs, .keep_all = TRUE)


### EXTRACCIÓN DE PICOS Y ALINEACIÓN ###########################################
################################################################################

# A partir de acá probamos trabajar con Spectra_list_prom_rep

# Análisis de la SNR en espectros para chequear que utilizamos el valor correcto
noise <- estimateNoise(Spectra_list_prom_rep[[20]])
plot(Spectra_list_prom_rep[[20]], xlim=c(4000, 20000), ylim=c(0, 0.002))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue") # Se ve que es correcto el 2

# Detección de picos a partir del umbral definido de SNR
peaks <- detectPeaks(Spectra_list_prom_rep,
                     SNR = SigNoi,
                     halfWindowSize = hws)

# Ploteo de picos detectados en un espectro de ejemplo
plot(Spectra_list_prom_rep[[20]], xlim=c(4000, 20000), ylim=c(0, 0.002))
points(peaks[[20]], col="red", pch=4)

# Alineado de picos
peaks <- alignPeaks(peaks,
                    minFreq = 0.8,
                    tolerance = tol)

#summaryPeaks(peaks[1:10])  # resumen estadistico de picos (primeros 10)

# Conteo de picos por perfil
cP <- countPeaks(peaks)

# Gráfico de picos
plot(cP, type = "n")
text(cP, label = 1:length(cP))


# Patrones de picos
peakPatterns(peaks)

# Filtrado de picos de baja frecuencia de aparición
picos_filtrados <- filterPeaks(peaks,
                               minFreq = 0.25,
                               labels = df_metadata_prom_rep$paciente ) #labels

# Patrones de picos
peakPatterns(picos_filtrados)

# Conteo de picos por perfil
cP2 <- countPeaks(picos_filtrados)

# Gráfico
plot(cP2, type = "n")
text(cP2, label = 1:length(cP2))

# Fusión de picos de la misma muestra
picos_fusion_muestra <- mergeMassPeaks(picos_filtrados,
                                       labels = df_metadata_prom_rep$f_pac_dabs,
                                       method = "median")

# Patrones de picos
peakPatterns(picos_fusion_muestra)


### CREACIÓN DE MATRIZ DE INTENSIDADES Y DICOTÓMICA ############################
################################################################################


# Matriz de intensidades 19 individuos
matint_61_ind <- intensityMatrix(picos_fusion_muestra,
                                 Spectra_list_prom_muestra) # sin valores NA

# Matriz de intensidades de 80 muestras
# matint_na_51 <- intensityMatrix(picos_fusion_muestra) # con valores NA
matint_179_mue <- intensityMatrix(picos_filtrados,
                                 Spectra_list_prom_rep) # sin valores NA

# Definición de umbrales
thr1 <- optimizeThreshold(matint_61_ind,
                          df_metadata_prom_mue$grupo,
                          verbose = T)
thr2 <- optimizeThreshold(matint_179_mue,
                          df_metadata_prom_rep$grupo,
                          verbose = T)

# Dicotomización
matint_61_ind_dico <- dichotomize(matint_61_ind, thr1)
matint_179_mue_dico <- dichotomize(matint_179_mue, thr2)

# Agrego nombres a las filas de cada df
rownames(matint_61_ind_dico) <- df_metadata_prom_mue$factor_dia
rownames(matint_179_mue_dico) <- df_metadata_prom_rep$factor_rep
rownames(matint_61_ind) <- df_metadata_prom_mue$factor_dia
rownames(matint_179_mue) <- df_metadata_prom_rep$factor_rep


### GUARDAR DATOS ##############################################################
################################################################################

# Establecer el directorio de trabajo en la ubicación del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Se puede mejorar, subcarpeta

# Guardar matrices y metadata asociada como archivo .Rdata
save(matint_61_ind_dico,df_metadata_prom_mue, file = "matint_61_ind_dico.Rdata")
save(matint_179_mue_dico, df_metadata_prom_rep, file = "matint_179_mue_dico.Rdata")
save(matint_61_ind,df_metadata_prom_mue, Spectra_list_prom_muestra, file = "matint_61_ind.Rdata")
save(matint_179_mue, df_metadata_prom_rep,  Spectra_list_prom_rep, file = "matint_179_mue.Rdata")


# Guardar matrices y metadata asociada como archivo .csv
write.csv(matint_61_ind_dico, "matint_61_ind_dico.csv", row.names = TRUE)
write.csv(matint_61_ind, "matint_61_ind.csv", row.names = TRUE)
write.csv(matint_179_mue_dico, "matint_179_mue_dico.csv", row.names = TRUE)
write.csv(matint_179_mue, "matint_179_mue.csv", row.names = TRUE)
write.csv(df_metadata_prom_mue, "df_61.csv", row.names = TRUE)
write.csv(df_metadata_prom_rep, "df_179.csv", row.names = TRUE)

# #
# #
# #
