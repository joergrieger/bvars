#' @export
#' @title set cholesky identification
#' @return returns an S3 object of the class chol

set_identification_cholesky <- function(){

  id <- structure(list(
    identification = "cholesky"

  ),class="chol")

  return(id)
}

#'
#' @title identify a vector-autoregressive model using cholesky identification
#' @param id_obj an S3 object containing information about the information
#' @param Alpha draw of coefficients
#' @param Sigma draw of variance-covariance matrix.
#' @return returns a K x K matrix with the identified variance-covariance


identify.chol <- function(id_obj,Alpha,Sigma){

  return(t(chol(Sigma)))

}

#' @export
#' @title set identification via sign restrictions
#' @param restrictions the sign restrictions
#' @return returns an S3 object of the class chol
#'
set_identification_sign <- function(restrictions){

  id <- structure(list(identification = "sign",
                       restrictions = restrictions),
                  class = "sign")
  return(id)
}

#' @title identify a vector-autoregressive model using sign-restrictions
#' @param id_obj an S3 object containing information about the information
#' @param Alpha draw of coefficients
#' @param Sigma draw of variance-covariance matrix.
#' @return returns a K x K matrix with the identified variance-covariance

identify.sign <- function(id_obj,Alpha,Sigma){

  K <- ncol(Sigma)
  SignRestriction <- FALSE
  cholsigma <- t(chol(Sigma))

  while(!SignRestriction){

    qrmatrix <- matrix(rnorm(K * K),nrow=K)
    qrdecomp <- qr(qrmatrix)
    qrdecomp <- qr.Q(qrdecomp)
    testmatrix <- qrdecomp %*% cholsigma
    SignRestriction <- !CheckSign(id_obj$restrictions,testmatrix)

  }

  Sigma <- testmatrix

  return(Sigma)

}

CheckSign <- function(RestrictionMatrix,TestMatrix){

  # Check if Restriction-Matrix and Test-Matrix are of the same sign

  Test1 <- dim(as.matrix(RestrictionMatrix))
  Test2 <- dim(as.matrix(TestMatrix))
  Test  <- Test1 == Test2

  if(!Test[1]){

    stop("Matrix with sign restrictions and test matrix do not have the same size")

  }

  if(!Test[2]){

    stop("Matrix with sign restrictions and test matrix do not have the same size")

  }
  n1 <- Test1[1]
  n2 <- Test1[2]


  # Loop over variance-covariance matrix
  TestFail=FALSE

  for(ii in 1:n1){

    for(jj in 1:n2){

      if(!is.na(RestrictionMatrix[ii,jj])){

        # Check if sign in [ii,jj] in Restrictio and Test-Matrix are bot the same

        if(RestrictionMatrix[ii,jj]<0){

          if(TestMatrix[ii,jj]>0){

            TestFail=TRUE

          }
        }
        if(RestrictionMatrix[ii,jj]>0){

          if(TestMatrix[ii,jj]<0){

            TestFail=TRUE

          }
        }
      }
    }
  }

  return(TestFail)

}
