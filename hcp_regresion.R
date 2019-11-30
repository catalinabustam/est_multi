library(MVN)
library(dplyr)
library(lmtest)
library(car)
library(randtests)
library(normtest)
library(corrplot)
library(het.test)
library(mvShapiroTest)
library(ResourceSelection)
library(psych)
library(caret)
library(magrittr)
library(plyr)
library(moonBook)
library(rsq)
library(stringr)
source("fun_regresion.R")



#################################################################

datos_join_i<-read.csv(file=file.choose(),header=T,sep=',')
head(datos_join_i)
attach(datos_join_i)

datos_hpc_data_trans <- read.csv(file=file.choose(),header=T,sep=',')
datos_clean = datos_hpc_data_trans
head(datos_clean)
attach(datos_clean)

#Se seleccionan los datos de interés
age_o<-Age
gender_o <-Gender
data <- datos_clean
r_list = c()
# Se eliminan los outliers
for(i in 1:length(data)){
  li = which(data[,i] %in% boxplot.stats(data[,i])$out)
  r_list <- c(r_list, li)
}
length(r_list)

new_data = data[-r_list,]
age <- age_o[-r_list]
gender <- gender_o[-r_list]

attach(new_data)

# Se escalan los datos
data_rescale <- new_data %>% 
                mutate_if(is.numeric, funs(as.numeric(scale(.))))

# Prueba de normalidad de de las variables
p <- list()
for (i in names(data.frame(data_rescale))) {
  column = get(i)
  p[i] <- as.numeric(unlist(ks.test(column,pnorm,mean(column),sd(column))["p.value"]))
}

################################# MANOVA #####################################################

# Se seleccionas solo las variables normales para el manova
norm_data <- data_rescale[names(which(p>0.05))]
p_norm <- list()
for (i in names(norm_data)) {
  column = get(i)
  p_norm[i] <- as.numeric(unlist(ks.test(column,pnorm,mean(column),sd(column))["p.value"]))
}

#Mutivariada

#**
datos_trans_matriz_no<-as.matrix(norm_data)
mvShapiro.Test(datos_trans_matriz_no)
#***
result <- mvn(datos_trans_matriz_no, mvnTest = "royston")
result
#****
resulta <- mvn(datos_trans_matriz_no, mvnTest = "mardia")
resulta$multivariateNormality
# No se logra probar la normalidad multivariada 

# gráficos iniciales

graphics.off()
par(mar=c(1,1,1,1))

par(mfrow=c(4,4))

for(i in 1:length(norm_data)){
  plot(norm_data[,i] ~ age)
  title(names(norm_data[i]))
}

par(mfrow=c(4,4))

for(i in 1:length(norm_data)){
  plot(norm_data[,i] ~ gender)
  title(names(norm_data[i]))
}

#Test de manova en las variables con el sexo y edad
beh = c("ProcSpeed_AgeAdj","CardSort_AgeAdj","AngAffect_Unadj","PercStress_Unadj")
beh_data = norm_data[beh]
dvars <- names(norm_data) %in% beh
vol_var = norm_data[!dvars ]

test_manova_sex_vol <- manova(data.matrix(norm_data) ~ gender)
summary(test_manova_sex_vol)
summary = summary.aov(test_manova_sex_vol)
summary

# Se corre la funcions para haer test de supuestos en manova
supuestos_manova_sex = supuestos_manova(test_manova_sex_vol,norm_data,gender )

## Se hace un anova para cada variable

anova_list_gender = data.frame()
for(i in seq(length(norm_data))){
  name = names(norm_data)[i]
  anova_1 = aov(norm_data[,i] ~ gender)
  supuestos = supuestos_manova(anova_1,norm_data[,i,drop=F],gender )
  row = data.frame(name= name, p_s = unlist(summary(anova_1)[[1]]$'Pr(>F)'<0.05)[[1]], assu = supuestos$all_assu )
  anova_list_gender <- rbind.fill(anova_list_gender,row)
}

anova_list_age = data.frame()
for(i in seq(length(norm_data))){
  name_2 = names(norm_data)[i]
  anova_2 = aov(norm_data[,i] ~ age)
  supuestos_2 = supuestos_manova(anova_2,norm_data[,i,drop=F],age )
  row2 = data.frame(name= name_2, p_s = unlist(summary(anova_2)[[1]]$'Pr(>F)'<0.05)[[1]], assu = supuestos_2$all_assu )
  anova_list_age<- rbind.fill(anova_list_age,row2)
}

################################# REGRESIÓN ##############################################


# Se prueba la regresión para todos los valores, buscando un significado lógico 

reg_frame<-data.frame()

columns_o = colnames(norm_data)


for(column_name in columns_o){
  newRow<-list()
  y_data <- column_name
  newRow$y_data <- column_name
  p_value <- ks.test(get(y_data),pnorm,mean(get(y_data)),sd(get(y_data)))["p.value"]
  newRow$p_value <- p_value
  newRow$norm <- (p_value>0.05)*1
  x_data <- subset(data_rescale, select =-get(y_data))
  x_data$age <- age
  x_data$gender <- gender

  #Estimo el modelo 

  f <- paste(y_data, "~",paste(names(x_data), collapse=" + "))
  MLG <- glm(as.formula(f))

  newRow$names = paste(names(which(summary(MLG)$coef[,"Pr(>|t|)"]<0.05)), collapse = ",")
  newRow$iac = summary(MLG)$aic
  newRow$r = rsq(MLG)
  newRow$r2 = rsq(MLG, adj= TRUE)
  
  reg_frame <- rbind.fill(reg_frame, as.data.frame(newRow))
}

write.csv(reg_frame, file="glm_hcp_2.csv",row.names=TRUE)

# Se selecciona la columna PercStress_Unadj por su interpretabilidad

y_data <- "PicVocab_AgeAdj_log"
p_value <- ks.test(get(y_data),pnorm,mean(get(y_data)),sd(get(y_data)))["p.value"]
x_data <- subset(datos_clean, select =-get(y_data))
x_data <- x_data[-(1)]
x_data$age <- age_o
x_data$gender <- gender_o

f <- paste(y_data, "~",paste(names(x_data), collapse=" + "))
MLG <- glm(as.formula(f), data=x_data)
summary(MLG)
plot(MLG)

# para el R

rsq(MLG) #PARA EL R
rsq(MLG, adj= TRUE) # PARA EL R2

# VALIDACIÓN DE SUPUESTOS DE LOS DATOS

#1. NO MULTICOLINEALIDAD 

num_data <- select_if(x_data, is.numeric)
corr_matrix_MM<-cor(num_data)
corrplot.mixed(corr_matrix_MM, lower="number", upper="ellipse")
  
# Como no pasa supuestos se hace análisis de factor analysis
## Se realiza factor analysis
nocorr_data_x<- get_factors(x_data)
out_names = colnames(nocorr_data_x)

#********nueva Regresión*********

x_data_n <- nocorr_data_x

f_n <- paste(y_data, "~",paste(names(x_data_n), collapse=" + "))
MLG_n <- glm(as.formula(f_n), data=x_data_n)
summary(MLG_n)
par(mfrow = c(2, 2))

rsq(MLG_n) #PARA EL R
rsq(MLG_n, adj= TRUE) # PARA EL R2
  

## Se seleccionan las variables significativas
#Numéricas
num_data_n <- data.frame(select_if(x_data_n, is.numeric))
rel_col <- names(which(summary(MLG_n)$coef[-1,"Pr(>|t|)"]<0.05))
col_names <- intersect(rel_col, names(num_data_n))
no_col_names <- setdiff(rel_col, names(num_data_n))
new_num <- x_data_n[col_names]

#Categóricas
factor_c <- data.frame(select_if(x_data_n, is.factor))
rel_factor <- lapply(names(factor_c), function(col) str_detect(no_col_names, col))
re_factor <- lapply(rel_factor, function(r) any(r))
new_factor <- factor_c[unlist(re_factor)]

sig_data <- data.frame(new_num,new_factor)

f_s <- paste(y_data, "~",paste(names(sig_data), collapse=" + "))
MLG_s <- glm(as.formula(f_s), data=sig_data)
summary(MLG_s)
par(mfrow = c(2, 2))
plot(MLG_s)

rsq(MLG_s) #PARA EL R
rsq(MLG_s, adj= TRUE) # PARA EL R2


#validación de supuestos
supuestos_glm(MLG_s, sig_data, get(y_data))


x_array <- transform_x(new_num)
y_array <- transform_y(get(y_data))

out_trans <- trans_supuestos_glm(y_array, x_array, factor= array("gender_o"))
  
par(mfrow = c(2, 2))
plot(MLG )

#**********************************************************************
#         REGRESION DE CADA VARIABLE DE COMPORTAMIENTO CON LOS VOLÚMENES
#**********************************************************************

beh = c("PicVocab_AgeAdj_log","Flanker_AgeAdj_T", "VSPLOT_TC", "ListSort_AgeAdj_T", "FearSomat_Unadj_log", "Sadness_Unadj_log", "MeanPurp_Unadj_T", "ProcSpeed_AgeAdj","CardSort_AgeAdj","AngAffect_Unadj","PercStress_Unadj")
vol = datos_clean[c(17:33,36:63)]
test_all_y = data.frame()

for(b in beh){
  # Se inicializan row dataframe de salida
  b= "PicVocab_AgeAdj_log"
  out_y = data.frame(y=b)
  out_assu <- data.frame(p_reset = NA, p_gqtest=NA, p_white =NA, p_runs=NA, p_ks=NA, p_jb=NA, r = NA, r2= NA, assu=NA)
  
  # Se obtiene la variable seleccionada
  be = get(b)
  
  # Se agregan variables categóricas
  vol$ageF <- age_o
  vol$genderF <- gender_o
  
  # Se unen las variables
  vol_v = cbind(vol, be)
  colnames(vol_v)[colnames(vol_v)=="be"] <- b
  
  # Se eliminan los valores atípicos
  r_list = c()
  for(i in 1:length(vol_v)){
    li = which(vol_v[,i] %in% boxplot.stats(vol_v[,i])$out)
    r_list <- c(r_list, li)
  }
  out_data = vol_v[-r_list,]
  
  # Se escalan los datos
  new_data <- out_data %>% 
    mutate_if(is.numeric, funs(as.numeric(scale(.))))
  
  # Se Selecciona la nueva y y se hace prueba de normalidad
  y_data = new_data[,c(48)]
  p_value <- ks.test(y_data,pnorm,mean(y_data),sd(y_data))["p.value"]
  out_y$y_norm_p  = p_value[[1]]
  out_y$y_norm = p_value[[1]] > 0.05
  
  #Se seleccionan las nuevas variables explicatorias
  x_data <- subset(new_data, select =-get(b))
  
  # Se aplica la función get factors para obtener factores por
  # Factor análisis en el caso de que exista multicolinealidad
  
  nocorr_data_x<- get_factors(x_data)
  out_names = colnames(nocorr_data_x)
  nocorr_data = cbind(nocorr_data_x, new_data[,c(48)])
  colnames(nocorr_data)[colnames(nocorr_data)=="new_data[, c(48)]"] <- b
  
  # Se hace un modelo lineal general
  
  f <- paste(b, "~",paste(out_names, collapse=" + "))
  MLG <- glm(as.formula(f), data= nocorr_data )
  summary(MLG)
  r = rsq(MLG)
  r2= rsq(MLG, adj= TRUE) 
  print(r)
  ## Se seleccionan las variables significativas
  
  rel_col <- names(which(summary(MLG)$coef[-1,"Pr(>|t|)"]<0.055))
  out_y$n_sig= length(rel_col)
  
  # Se sigue el proceso si se encuentra al menos una variabale significativa
  
  if(length(rel_col)>0){
    
    #Numéricas
    col_names <- intersect(rel_col, out_names)
    no_col_names <- setdiff(rel_col, out_names)
    no_c_name = c()
    for(n in no_col_names){
      n_name = paste(strsplit(n, "F")[[1]][1], "F", collapse="", sep = "")
      no_c_name = append(no_c_name, n_name)
    }
    
    f_names = unique(no_c_name)
    
    new_s_names <- append(col_names, f_names)
    
    # Se hace de nuevo el modelo con las variabales significativas
    
    f_s <- paste(b, "~",paste(new_s_names, collapse=" + "))
    MLG_s <- glm(as.formula(f_s), data= nocorr_data)
    summary(MLG_s)
    par(mfrow = c(2, 2))
    plot(MLG_s)
    
    x_data = nocorr_data[, c(new_s_names), drop=F]
    y_data = nocorr_data[, c(b)]
    
    # Se corre la función que verifica los supuestos para un mlg
    
    out_assu = supuestos_glm(MLG_s, x_data, y_data)
  }
  out_y = cbind(out_y, out_assu)
  test_all_y = rbind(test_all_y,out_y)
}
#*************************************************************************
