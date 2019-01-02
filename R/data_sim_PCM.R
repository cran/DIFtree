#' @name data_sim_PCM
#' @title Simulated Data Set with Polytomous Items  
#' @description The data set is simulated from a Partial Credit Model  
#' where some items exhibit differential item functioning. 
#' Existing differences in item difficulties are simulated by step-functions. 
#' The true, simulated DIF structure is described in Bollmann et al. (2017), Section 4.3. 
#' @docType data
#' @usage data(data_sim_PCM)
#' @format A data frame containing 500 observations on 4 variables:
#' 
#'  \bold{Y}   matrix with categorical responses (3-point scale)
#'  
#'  \bold{x1}  binary covariate 
#'
#'  \bold{x2}  ordinal covariate
#' 
#'  \bold{x3}  numeric covariate 
#' 
#' 
#' @references 
#' Bollmann, Stella, Berger, Moritz & Tutz, Gerhard (2018): Item-Focussed Trees for the Detection 
#' of Differential Item Functioning in Partial Credit Models, Educational and Psychological Measurement 78(5), 781-804.
#' 
#' @examples 
#' data(data_sim_PCM)
#'  
#' Y <- data_sim_PCM[,1]
#' X <- data_sim_PCM[,-1]
#' 
#' apply(Y,2,table)
#' summary(X)
#'   
#'  
NULL
