library(lmtest)
library(normtest)
library(randtests)
library(corrplot)


supuestos_manova <- function(manova, data, group) {
  #1. Residuos aleatorios: gráficos
  #Wald-Wolfowitz 
  residuos<-residuals(manova)
  plot(residuos)+ abline(h=0, col="red")
  p_runs <- runs.test(residuos)$p.value
  
  #2. Normalidad de residuos
  p_jb = jb.norm.test(residuos, nrepl=2000)$p.value
  res<-as.matrix(residuos)
  
  out <- data.frame(p_runs= p_runs, p_jb= p_jb)
  
  if(dim(res)[2] > 1){
    mv_res <- mvn(res, mvnTest = "hz") 
    out$p_multi_hz= mv_res$multivariateNormality["p value"]
    out$no_normal_uni <- paste(which(mv_res$univariateNormality$Normality== "   NO    "),collapse="-")
  }

  out$p_ks = ks.test(res,pnorm,mean(res),sd(res))$p.value
  p_ks = out$p_ks
  #3. Homoscedasticidad
  
  p_leve <- list()
  
  for (i in names(data)) {
    column = get(i)
    p_leve[i] <- leveneTest(y = column, group = group)["Pr(>F)"][[1]][1]
  }
  
  #out$p_bar<-bartlett.test(manova$residuals ~ group)$p.value
  
  out <- cbind(out,p_leve)
  out$all_assu = (p_runs>0.05)&(p_jb>0.05|p_ks>0.05)&(all(unlist(p_leve)>0.05))
  return(out)
}

supuestos_glm <- function(glm, data, y) {
  ##1. Especificaci?n
  p1<- reset(glm)$p.value
  
  #2. Homosedasticidad
  #Goldfeld-Quandt Test
  p2 <- gqtest(glm)$p.value
  
  #White's test
  num_data <- select_if(data, is.numeric)
  dataset <- data.frame(y, num_data)
  model1 <- VAR(dataset, p = 1)
  p21<-whites.htest(model1)$p.value
  
  #3. No autocorrelaci?n: residuos independientes
  res_MLG<-residuals(glm)
  
  p3<-runs.test(res_MLG)$p.value
  
  #4. Normalidad de los residuos
  
  p4<- ks.test(res_MLG,pnorm,mean(res_MLG),sd(res_MLG))$p.value
  p41 <- jb.norm.test(res_MLG, nrepl=2000)$p.value

  r = rsq(MLG)
  r2= rsq(glm, adj= TRUE) 

  assu = (p1>0.05)&(p2>0.05|p21>0.05)&(p3>0.05)&(p4>0.05|p41>0.05)
  
  outrow <- data.frame(p_reset = p1, p_gqtest=p2, p_white =unname(p21), p_runs=p3, p_ks=p4, p_jb=p41, r = r, r2= r2, assu=unname(assu))
  return(outrow)
}


transform_x = function(data){
  k = length(data); n=nrow(data); m = 3
  big_array = array(NA, c(n,m,k))
  for(i in seq(k)){
    big_array[,1,i]= new_num[,i]
    big_array[,2,i]= 1/new_num[,i]
    big_array[,3,i]= 1/(new_num[,i]+1)
  }
  return(big_array)
}

transform_y = function(y_val){
  n = length(y_val);
  y_array = array(NA, c(n,8))
  y_array[,1] = y_val
  y_array[,2] = sqrt(y_val)
  y_array[,3] = sqrt(y_val)+ sqrt(y_val +1)
  y_array[,4] = log(y_val)
  y_array[,5] = log(y_val+1)
  y_array[,6] = 1/y_val
  y_array[,7] = 1/(y_val+1)
  y_array[,8] = 1/sin(sqrt(y_val))
  return(y_array)
}

trans_supuestos_glm = function(y_array, x_array, factor_v= c(), current_index=1, iter_list = vector(), out_array = array(),out_names=vector(), output = data.frame(), number_of_loops=0, range_list= c(), stop_f=0) {
  #https://stackoverflow.com/questions/7186518/function-with-varying-number-of-for-loops-python
  
  if(!length(iter_list)){
    n= dim(x_array)[1]
    number_of_loops = dim(x_array)[3]  + 1
    range_list <- vector("list", number_of_loops)
    range_list[[1]] = seq(dim(y_array)[2])
    for(d in seq(dim(x_array)[3])){range_list[[d+1]]<- seq(dim(x_array)[2])}
    iter_list = c(1:number_of_loops)*0
    out_array = array(NA, c(n,(number_of_loops-1)))
    out_names = c((1:(number_of_loops-1)))*0
  }
  
  if(current_index == number_of_loops){
    for(i in range_list[[current_index]]){
      iter_list[current_index]<-i 
      y_i = iter_list[1]
      x_i = iter_list[-1]
    
      for(j in seq(length(x_i))){
        out_array[,j]= x_array[,x_i[j],j]
        out_names[j] = paste("out_array[,",j, "]", collapse="")
      }
      y_name = paste("y_array[,c(", y_i,")]")
      names_with_factor <- as.vector(rbind(out_names,factor_v)) 
      names_with_factor <-  names_with_factor[!is.na( names_with_factor)]
      out_name_w_f <- paste(names_with_factor, collapse=" + ")
      f <- paste(y_name, "~",out_name_w_f)
      MLG_i <- glm(as.formula(f))
        
      #1. Especificaci?n
      p1<- reset(MLG_i)$p.value

      #2. Homosedasticidad
      #Goldfeld-Quandt Test
      p2 <- gqtest(MLG_i)$p.value
        
      #3. No autocorrelaci?n: residuos independientes
      res_MLG<-residuals(MLG_i)
        
      p3<-runs.test(res_MLG)$p.value
        
      #4. Normalidad de los residuos

      p4<- ks.test(res_MLG,pnorm,mean(res_MLG),sd(res_MLG))$p.value
      r = rsq(MLG_i)
      r2= rsq(MLG_i, adj= TRUE) 
    
      if(all(c(p1,p2,p3,p4)>0.05)){
        r = rsq(MLG_i)
        r2= rsq(MLG_i, adj= TRUE) 
        outrow <- data.frame(happy_list = paste(iter_list, collapse =""), p_reset = p1, p_gqtest=p2, p_runs=p3, p_ks=p4, r = r, r2= r2)
        output <- rbind.fill(output,outrow)
      }
    }
    if(stop_f==1){
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  else{
    for(i in range_list[[current_index]]){
      iter_list[current_index]<-i
      output <- trans_supuestos_glm(y_array, x_array, factor_v= factor_v, current_index = current_index+1, iter_list = iter_list, out_array = out_array,out_names = out_names,output=output, number_of_loops= number_of_loops,range_list= range_list)
    }
  }
  return(output)
}
  
#generate_cordenates(y_array, big_array, factor= array("gender"))


get_factors <- function(x_data){
  #1. NO MULTICOLINEALIDAD 
  factor <- data.frame(select_if( x_data, is.factor))
  num_data <- select_if(x_data, is.numeric)
  corr_matrix_MM<-cor(num_data)
  
  corre_var = findCorrelation(corr_matrix_MM, cutoff = 0.5, verbose = FALSE, names = TRUE)
  
  if(corre_var > 2){
    df_corretated <- num_data
    corr_matrix_co<-cor(df_corretated)
    parallel_beh <- fa.parallel(df_corretated, fm='minres', fa='fa')
    n_f = parallel_beh$nfact
    factorAnalysis<-fa(df_corretated,nfactors=n_f,fm="ml",scores="regression",rotate="varimax")
    fa.diagram(factorAnalysis$loadings)
    scores_MM<-factorAnalysis$scores
    scores_df = as.data.frame(scores_MM)
    return(cbind(scores_df,factor))
  }
  else{
    return(cbind(x_data,factor))
  }
}
  
