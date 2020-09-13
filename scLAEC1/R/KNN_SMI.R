KNN_SMI <- function(W,K=10,t=50) {
  normalize <- function(X) X / rowSums(X)
  newW =matrix(0,dim(W)[1], dim(W)[2])
  nextW =matrix(0,dim(W)[1], dim(W)[2])
  W = normalize(W);
  W = (W+t(W))/2;

  newW = (.dominateset(W,K))

  nextW=W;
  for (i in 1:t) {
    nextW=newW %*% (nextW) %*% t(newW);
  }

  W = nextW + diag(nrow(W));
  W = (W + t(W))/2;


  return(W)
}