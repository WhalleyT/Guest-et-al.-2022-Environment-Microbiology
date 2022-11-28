## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

library(dplyr) 
library(stringr)
library(spatstat)
library(broom)
library(data.table)
library(tidyr)
library(diptest)
library(pbapply)
library(multimode)
library(clustertend)
library(parallel)


## ----get files, message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------------------------------------------------------------


files <- list.files("new_csvs/", full.names = T, pattern = "\\.green.csv$")
grn   <- lapply(files, read.csv)

annotate_data <- function(df, name){
  stain <- str_match_all(name, "(?<=\\().+?(?=\\))")[[1]][,1]
  split <- strsplit(name, " ")[[1]]
  
  if(nrow(df) > 0){
    df$stain <- stain
    df$strain <- strsplit(split[1], "//")[[1]][1]
    df$biocide <- split[3]
    df$time_growth <- split[5]
    df$treatment_time <- split[3]
    df$file_name <- name
    df$sample <- strsplit(split[6], "\\.")[[1]][1]
  }
  return(df)
}

annotate_data <- function(df, name){

  split <- strsplit(name, " ")[[1]]
  stain <- gsub(".csv", "", split[4])
  stain <- gsub("^.*\\.","", stain)
  
  if(nrow(df) > 0){
    df$stain <- stain
    df$strain <- strsplit(split[1], "//")[[1]][1]
    df$biocide <- split[3]
    df$time_growth <- split[5]
    df$treatment_time <- split[3]
    df$file_name <- name
    df$sample <- strsplit(split[6], "\\.")[[1]][1]
    return(df)
  }
  
}


grn_dfs <- mapply(annotate_data, grn, files, SIMPLIFY = FALSE)

grn_data <- do.call("rbind", grn_dfs[c(2,3,4,5,6,7,9,10,11,12,13,14)])


## ----fix names, message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------------------------------------------------------------
grn_data$name <- gsub(" green.csv", "", grn_data$file_name)
grn_data$name <- gsub("data/", "", grn_data$name)


## ----quadrat------------------------------------------------------------------------------------------------------------------------------------------------------------------------

power_dip <- function(df, nsim=1000, sig = 0.05){
  
  if("XM" %in% colnames(df) & "YM" %in% colnames(df)){
    numerics <- df %>% select(XM, YM)
    pvals    <- replicate(n = nsim, broom::tidy(dip.test(as.matrix(numerics), simulate.p.value = T))$p.value)
    name <- unique(df$file_name)
    name <- gsub("new_csvs/", "simulation_tmp/dip/", name)
         
    write.csv2(pvals, name)
    return(sum(ifelse(pvals <= sig, T, F)) / length(pvals))
  }else{
    return(NA)
  }
}

power_quad <- function(df, nsim=1000, sig = 0.05){
  win  <- owin(xrange = c(0, 1050), yrange = c(0, 1050))
  obj  <- ppp(df$XM, df$YM, window = win)
  pvals <- replicate(n= nsim, broom::tidy(quadrat.test(obj))$p.value)
  
  name <- unique(df$file_name)
  name <- gsub("new_csvs/", "simulation_tmp/quadrat/", name)
  
  if(length(name) != 0){
     write.csv2(pvals, name)
  }
 
  return(sum(ifelse(pvals <= sig, T, F)) / length(pvals))
}


power_silverman <- function(df, nsim=1000, sig = 0.05){
    if("XM" %in% colnames(df) & "YM" %in% colnames(df)){
    numerics <- df %>% select(XM, YM)
    pvals    <- replicate(n = nsim, broom::tidy(modetest(as.matrix(numerics), method = "SI"))$p.value)
    
    name <- unique(df$file_name)
    name <- gsub("new_csvs/", "simulation_tmp/silverman/", name)
         
    if(length(name) != 0){
       write.csv2(pvals, name)
    }
         
    return(sum(ifelse(pvals <= sig, T, F)) / length(pvals))
  }else{
    return(NA)
  }
}

power_hopkins <- function(df, nsim=1000, sig = 0.05){
     if("XM" %in% colnames(df) & "YM" %in% colnames(df)){
       numerics <- df %>% select(XM, YM)
       if(nrow(na.omit(numerics)) > 10){
         nsamp  <- round(nrow(na.omit(numerics)) * 0.1)
         scores <- replicate(n= nsim, hopkins(numerics, n=nsamp)$H)
         
         name <- unique(df$file_name)
         name <- gsub("new_csvs/", "simulation_tmp/hopkins/", name)
         
         if(length(name) != 0){
            write.csv2(scores, name)
         }
         
         return(sum(scores<qbeta(sig, nsamp, nsamp)) / nsim)
       }else{
         return(NA)
       }
     }else{
       return(NA)
     }
}


## ----call functions-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#apply and cbind
number_sims <- 10
significance <- 0.01
print("calling functions")
dip_powers     <- pblapply(grn_dfs, power_dip, nsim=number_sims, sig = significance)
quadrat_powers <- pblapply(grn_dfs, power_quad, nsim=number_sims, sig = significance)
hopkins_powers <- pblapply(grn_dfs, power_hopkins, nsim=number_sims,  sig = significance)
silverman_powers <- lapply(grn_dfs, power_silverman, nsim=number_sims,  sig = significance)



## ----output-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
outdf <- do.call(rbind, Map(data.frame, 
                            "Hopkins (1000 runs, significance = 0.01)" = hopkins_powers, 
                            "Quadrat (1000 runs, significance = 0.01)" = quadrat_powers,
                            "Dip (1000 runs, significance = 0.01)" = dip_powers,
                            "Silverman (1000 runs, significance = 0.01)" = silverman_powers))


outdf$file <- gsub("new_csvs/", "", files)
write.csv(outdf, "power_estimates.csv")

