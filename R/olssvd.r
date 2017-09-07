olssvd <- function(y,ly){
  duv          <- svd(t(ly)%*%ly)
  x.inv        <- duv$v%*%diag(1/duv$d)%*%t(duv$u)
  x.pseudo.inv <- x.inv %*% t(ly)
  b <- x.pseudo.inv%*%y
  return(b)
}