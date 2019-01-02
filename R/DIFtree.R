#' Item focussed Trees for the Identification of Items in Differential Item Functioning 
#' 
#' @description
#' A function to estimate item focussed trees for simultaneous selection of items and variables 
#' that induce DIF (Differential Item Functioning) in dichotomous or polytomous items. DIF detection can be
#' based on the Rasch Model (dichotomous case), the Logistic Regression Approach (dichotomous case) or the Partial Credit Model (polytomous case).
#' The basic method of item focussed recursive partitioning in Rasch Models is described in Tutz and Berger (2015).
#' 
#' @param Y Matrix or Data.frame of binary 0/1 or categorical response (rows correspond to persons, columns correspond to items)
#' @param X Data.frame of (not scaled) covariates (rows correspond to persons, columns correspond to covariates)
#' @param model Type of model to be fitted; can be \code{"Rasch"}, \code{"Logistic"} or \code{"PCM"}.
#' @param type  Type of DIF to be modelled; one out of \code{"udif"}, \code{"dif"} and \code{"nudif"}. 
#' For \code{"Rasch"} and \code{"PCM"} only uniform DIF can be modelled and therefore \code{type} will be ignored.
#' @param alpha Global significance level for the permutation tests
#' @param nperm Number of permutations used for the permutation tests
#' @param trace If true, information about the estimation progress is printed
#' @param penalize If true, a small ridge penalty is added to ensure existence of model parameters; only for \code{"Rasch"}.
#' @param x Object of class \code{"DIFtree"}
#' @param ... Further arguments passed to or from other methods
#' 
#' @details 
#' The methods require 0/1 coded answers on binary items (\code{"Rasch"} and \code{"Logistic"}) or categorical answers on polytomous items (\code{"PCM"}). 
#' Items with DIF are gradually identified by recursive partitioning.
#' 
#' For \code{"Rasch"} one yields a model with linear predictors 
#' \deqn{eta_{pi}=theta_p-tr_i(x_p),}
#' where \eqn{theta_p} correspond to the ability and \eqn{x_p} correspond to the covariate vector of person p. 
#' 
#' For \code{"Logistic"} one yields a model with linear predictors 
#' \itemize{
#' \item Uniform DIF, \code{type="udif"}
#' \deqn{eta_{pi}=S_p beta_i+tr_i(x_p),}
#' where \eqn{S_p} corresponds to the test score and \eqn{x_p} corresponds to the covariate vector of person p.
#' \item DIF and Non-Uniform DIF, \code{type="dif", "nudif"}
#' \deqn{eta_{pi}=tr_i(x_p)+tr_i(S_p,x_p),}
#' where \eqn{S_p} corresponds to the test score and \eqn{x_p} corresponds to the covariate vector of person p. 
#' }
#' 
#' For \code{"PCM"} one yields a model with linear predictors 
#' \deqn{eta_{pir}=theta_p-tr_{ir}(x_p),}
#' where \eqn{theta_p} correspond to the ability and \eqn{x_p} correspond to the covariate vector of person p.
#' 
#'  
#' 
#' Significance of each split is verified by permutation tests. The result of the permutation tests 
#' can strongly depend on the number of permutations \code{nperm}.
#' In the case of pure terminal nodes estimates of the model do not exist. If \code{penalize=TRUE} 
#' a small ridge penalty is added during estimation to ensure existence of all parameters. 
#'
#' @return Object of class \code{"DIFtree"}. 
#' An object of class \code{"DIFtree"} is a list containing the following components:
#' 
#' \item{splits}{Matrix with detailed information about all executed splits during the estimation process}
#' \item{coefficients}{List of estimated coefficients for items with and without DIF. 
#' Structure of \code{coefficients} depends on \code{model} and \code{type}.}
#' \item{pvalues}{P-values of each permutation test during the estimation process}
#' \item{devs}{Maximal value statistics \eqn{T_j} of the selected variables in each iteration during the estimation process}
#' \item{crit}{Critical values of each permutation test during the estimation process}
#' \item{Y}{Response matrix used in the estimation process}
#' \item{X}{Model matrix used in the estimation process}
#' \item{persons}{Number of persons} 
#' \item{items}{Number of items} 
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
#' Swaminathan, Hariharan and Rogers, H Jane (1990): Detecting differential item functioning 
#' using logistic regression procedures, Journal of Educational Measurements 27(4), 361-370.
#' 
#' Tutz, Gerhard and Berger, Moritz (2016): Item focussed Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika 81(3), 727-750.
#' 
#' @seealso \code{\link[DIFtree]{plot.DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}, \code{\link[DIFtree]{summary.DIFtree}}
#' 
#' @examples 
#' data(data_sim_Rasch)
#' data(data_sim_PCM)
#'  
#' Y1 <- data_sim_Rasch[,1]
#' X1 <- data_sim_Rasch[,-1]
#' 
#' Y2 <- data_sim_PCM[,1]
#' X2 <- data_sim_PCM[,-1]
#'  
#' \dontrun{
#'  
#' mod1 <- DIFtree(Y=Y1,X=X1,model="Logistic",type="udif",alpha=0.05,nperm=1000,trace=TRUE)
#' print(mod1)
#' 
#' mod2 <- DIFtree(Y=Y2,X=X2,model="PCM",alpha=0.05,nperm=100,trace=TRUE)
#' print(mod2)
#' }
#' 
#
#' @exportClass DIFtree
#' @export
#' @importFrom penalized penalized 
#' @importFrom VGAM vglm acat
#' @importFrom stats binomial coef deviance formula glm predict quantile var coefficients na.omit


DIFtree <-
function(Y,
                    X,
                    model=c("Rasch","Logistic","PCM"),
                    type=c("udif","dif","nudif"),
                    alpha=0.05,
                    nperm=1000,
                    trace=FALSE,
                    penalize=FALSE,
                    ...){
  UseMethod("DIFtree")
}

#' @rdname DIFtree 
#' @method print DIFtree
#' @export

print.DIFtree <-
  function(x, # object of class DIFtree 
           ...){
    
    model <- which(c("Rasch","Logistic","PCM")==x$model)
    if(model==2){
      type <- which(c("udif","dif","nudif")==x$type)
    } else{
      type <- 1
    }
    
    npersons  <- x$persons
    nitems    <- x$items
    
    if(model==2 & type==2){
      dif_items <- unique(c(x$splits[[1]][,"item"],x$splits[[2]][,"item"]))
    } else{
      dif_items <- unique(x$splits[,"item"])
    }
    
    if(model==2 & type==2){
      splits <- lapply(1:2,function(j) x$splits[[j]][,c("item","variable","threshold")])
    } else{
      splits    <- x$splits[,c("item","variable","threshold")]
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
    cat("Number of persons:", npersons, "\n")
    cat("Number of items:", nitems, "\n")
    if(!is.null(dif_items)){
      cat("DIF items:", dif_items, "\n")
    } else{
      cat("DIF items: no DIF item\n")
    }
    cat("\n")
    cat("Overview of executed splits:\n")
    if(!is.null(dif_items)){
      cat("\n")
      if(model==2 & type==2){
        cat("Intercept:\n")
        if(!is.null(splits[[1]])){
          print(splits[[1]])
        } else{
          cat("no split performed")
        }
        cat("\n")
        cat("Slope:\n")
        if(!is.null(splits[[2]])){
          print(splits[[2]])
        } else{
          cat("no split performed")
        }
      } else{
        print(splits)
      }
    } else{
      cat("no split performed")
    }
    invisible(x)
  }  
