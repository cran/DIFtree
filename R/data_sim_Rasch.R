#' @name data_sim_Rasch
#' @title Simulated Data Set with Dichotomous Items  
#' @description The data set is simulated from a Rasch model 
#' where some items exhibit differential item functioning. 
#' Existing differences in item difficulties are simulated by step-functions. 
#' The true, simulated DIF structure is described in Tutz and Berger (2015), Section 4.2. 
#' @docType data
#' @usage data(data_sim_Rasch)
#' @format A data frame containing 500 observations on 5 variables:
#' 
#'  \bold{Y}   matrix with binary 0/1 response for 20 items
#'  
#'  \bold{x1}  binary covariate 1
#'
#'  \bold{x2}  metric covariate 1
#' 
#'  \bold{x3}  binary covariate 2
#' 
#'  \bold{x4}  metric covariate 2
#' 
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2016): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed TreesJournal of Educational and Behavioral Statistics 41(6), 559-592.
#' 
#' Tutz, Gerhard and Berger, Moritz (2016): Item focussed Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika 81(3), 727-750. 
#' 
#' @examples 
#' data(data_sim_Rasch)
#'  
#' Y <- data_sim_Rasch[,1]
#' X <- data_sim_Rasch[,-1]
#' 
#' hist(rowSums(Y), breaks = 0:19 + 0.5)
#' summary(X)
#'   
#'  
NULL
