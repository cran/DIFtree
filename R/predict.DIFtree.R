#' Prediction of item difficulties of Item focussed Trees 
#' 
#' @description
#' The function returns predictions of item difficulties 
#' based on the estimated trees by item focussed recursive partitioning.  
#' 
#' @param object object of class \code{\link[DIFtree]{DIFtree}}
#' @param items number of items, for which the prediction shall be returned  
#' @param X new data.frame, for which the prediction shall be returned 
#' @param ... further arguments passed to or from other methods 
#' 
#' @author Moritz Berger <moritz.berger@@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @references Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Department of Statistics, LMU Munich
#' 
#' @seealso \code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{plot.DIFtree}}
#' 
#' @examples 
#' 
#' data(data_sim)
#'  
#' Y <- data_sim[,1]
#' X <- data_sim[,-1]
#' 
#' Xnew <- data.frame("x1"=c(0,1),"x2"=c(-1.1,2.5),"x3"=c(1,0),"x4"=c(-0.2,0.7))
#'  
#' \dontrun{
#'  
#' mod <- DIFtree(Y=Y,X=X,alpha=0.05,nperm=100,trace=TRUE)
#'  
#' predict(mod,items=c(1,2,3),Xnew)
#' }
#' 
#
#' @method predict DIFtree
#' @export


predict.DIFtree <-
function(object,items,X,...){
  
  if(!is.data.frame(X)){
    stop("X must be of class 'data.frame'")
  }
  if(ncol(X)!=ncol(object$X)){
    stop("X does not have the right dimension")
  }
  if(!all(names(X) %in% names(object$X))){
    stop("One or more variables are missing")
  }
  
  nitems       <- object$items
  
  dif_items    <- unique(object$splits[,"item"])
  if(!is.null(dif_items)){
    no_dif_items <- seq(1,nitems)[-dif_items]
  } else{
    no_dif_items <- seq(1,nitems)
  }
  
  one_betas <- function(x) {
    betas_hat <- numeric(length=length(items))
    
    for(i in 1:length(items)){
      if(items[i] %in% no_dif_items){
        betas_hat[i] <- object$beta_hat_nodif[paste0("beta",items[i])]
      } else{
        string <- c()
        info <- object$splits[object$splits[,"item"]==items[i],]
        help1 <- strsplit(names(object$beta_hat_dif[[paste(items[i])]]),"")
        help2 <- sapply(1:length(help1), function(j) paste0(help1[[j]][which(help1[[j]]=="_")+1],collapse=""))
        done <- FALSE
        j <- 1 
        while(!done){
          row <- info[info$number==j,]
          if(nrow(row)==1){
            var <- row$variable
            thres <- row$threshold
            if(x[var] <= thres){
              string <- paste0(string,"l")
              j <- row$left
            } else{
              string <- paste0(string,"u")
              j <- row$right
            }
          } else{
            done <- TRUE 
          }
        }
        which_coef <- which(string==help2)
        betas_hat[i] <- object$beta_hat_dif[[paste(items[i])]][which_coef]
      }
    }
  
    return(betas_hat)
  }

  betas_hat <- matrix(t(sapply(1:nrow(X),function(j) one_betas(X[j,]))),nrow=nrow(X),ncol=length(items))
  rownames(betas_hat) <- 1:nrow(X)
  colnames(betas_hat) <- paste("Item", items)
  return(betas_hat)
  
  invisible(object)
  
}
