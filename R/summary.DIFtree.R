#' Summary for fitted Item focussed Trees 
#' 
#' @description
#' The function takes an object of class \code{"DIFtree"} and returns an useful summary 
#' with an overiew of all executed splits during the estimation procedure.
#' 
#' @param object Object of class \code{\link[DIFtree]{DIFtree}}
#' @param x Object of class \code{\link[DIFtree]{summary.DIFtree}}
#' @param ... Further arguments passed to or from other methods 
#' 
#' @return Object of class \code{"summary.DIFtree"}. 
#' An object of class \code{"summary.DIFtree"} is a list containing the following components:
#' 
#' \item{stats}{Useful overview of detected DIF items, responsible variables and executed splits}
#' \item{nosplits}{Total number of executed splits during the estimation procedure}
#' 
#' @author Moritz Berger <moritz.berger@imbie.uni-bonn.de> \cr \url{http://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2016): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed Trees, Journal of Educational and Behavioral Statistics 41(6), 559-592.
#' 
#' Bollmann, Stella, Berger, Moritz & Tutz, Gerhard (2018): Item-Focussed Trees for the Detection 
#' of Differential Item Functioning in Partial Credit Models, Educational and Psychological Measurement 78(5), 781-804.
#' 
#' Tutz, Gerhard and Berger, Moritz (2016): Item focussed Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika 81(3), 727-750. 
#' 
#' @seealso \code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{plot.DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}
#' 
#' @examples 
#' data(data_sim_Rasch)
#'  
#' Y <- data_sim_Rasch[,1]
#' X <- data_sim_Rasch[,-1]
#'  
#' \dontrun{
#'  
#' mod <- DIFtree(Y=Y,X=X,model="Logistic",type="udif",alpha=0.05,nperm=1000,trace=TRUE)
#'  
#' summary(mod)
#' }
#' 
#
#' @method summary DIFtree
#' @exportClass summary.DIFtree 
#' @export

summary.DIFtree <-
function(object, # object of class DIFtree
                            ...){
  
  to_return <- list(call=object$call)
  
  model <- which(c("Rasch","Logistic","PCM")==object$model)
  if(model==2){
    type <- which(c("udif","dif","nudif")==object$type)
  } else{
    type <- 1
  }
  nitems <- object$items
  overview <- infos_summary(object$splits,1:nitems,model,type)
  
  if(model==2 & type==2){
    nos <- nrow(rbind(object$splits[[1]],object$splits[[2]]))
  } else{
    nos <- nrow(object$splits)
  }
  to_return$stats  <- overview 
  to_return$nosplits <- nos  
  to_return$model <- object$model
  to_return$type  <- object$type

  class(to_return) <- "summary.DIFtree"
  to_return
  
}

#' @rdname summary.DIFtree 
#' @method print summary.DIFtree
#' @export

print.summary.DIFtree <-
  function(x, # object of class summary.DIFtree 
           ...){
    
    model <- which(c("Rasch","Logistic","PCM")==x$model)
    if(model==2){
      type <- which(c("udif","dif","nudif")==x$type)
    } else{
      type <- 1
    }
    cat("\n")
    if(model==1){
      cat("Item focussed Trees based on the Rasch Model:\n")
    }
    if(model==2){
      cat("Item focussed Trees based on the Logistic Regression Approach",c("(uniform","(uniform and non-uniform","(non-uniform")[type],"DIF):\n")
    }
    if(model==3){
      cat("Item focussed Trees based on the PCM:\n")
    }
    cat("\n")
    cat("Call:\n",paste(deparse(x$call)),"\n")
    cat("\n")
    cat("----------\n")
    cat("\n")
    cat("Overview:\n")
    cat("\n")
    print(x$stats)
    cat("\n")
    cat("Total number of Splits:", x$nosplits)
  }
