bivar_quantile <- function(n, p, a, E.a, y, l.full, K, V.inv, theta.tilde, test.size, nrep, M, r){
  ## theta.tilde := (only the selected psi's, delta)
  if(p <= 60){
    comp1 <- lapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%y[[i]])
    comp2 <- lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l.full[[i]])%*%V.inv[[i]]%*%y[[i]])
    comp3 <- lapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%l.full[[i]])
    comp4 <- lapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l.full[[i]])) #(c(a[[i]]-E.a[[i]])*l[[i]]))
    comp5 <- lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l.full[[i]])%*%V.inv[[i]]%*%l.full[[i]])
    comp6 <- lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l.full[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l.full[[i]])) #(c(a[[i]]-E.a[[i]])*l[[i]]))
    Z <- lapply(1:n, function(i) c(c(comp1[[i]]), c(comp2[[i]]),
                                   comp3[[i]][upper.tri(comp3[[i]], diag=TRUE)],
                                   comp4[[i]][upper.tri(comp4[[i]], diag=TRUE)],
                                   comp5[[i]][upper.tri(comp5[[i]], diag=TRUE)],
                                   comp6[[i]][upper.tri(comp6[[i]], diag=TRUE)]))
    Z.bar <- Reduce("+", Z)/n
    S.boot.i <- lapply(1:nrep, function(B)
      lapply(1:n, function(i) r[[B]][i]*(Z[[i]] - Z.bar)/sqrt(n)))
    S.boot <- lapply(1:nrep, function(B) Reduce("+", S.boot.i[[B]]))
    S.boot.I <- unlist(lapply(1:nrep, function(B) max(abs(S.boot[[B]][1:(2*K)]))))
    S.boot.II <- unlist(lapply(1:nrep, function(B) max(abs(S.boot[[B]][(2*K+1):(2*K+4*(K+K*(K-1)/2))]))))
  } else {
    comp1 <- sapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%y[[i]])
    comp2 <- sapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l.full[[i]])%*%V.inv[[i]]%*%y[[i]])
    comp3 <- lapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%l.full[[i]])
    comp4 <- lapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l.full[[i]])) #(c(a[[i]]-E.a[[i]])*l[[i]]))
    comp5 <- lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l.full[[i]])%*%V.inv[[i]]%*%l.full[[i]])
    comp6 <- lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l.full[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l.full[[i]])) #(c(a[[i]]-E.a[[i]])*l[[i]]))
    Z.I <- rbind(comp1, comp2)
    Z.II <- rbind(sapply(1:n, function(i) c(comp3[[i]][upper.tri(comp3[[i]], diag=TRUE)])),
                  sapply(1:n, function(i) c(comp4[[i]][upper.tri(comp4[[i]], diag=TRUE)])),
                  sapply(1:n, function(i) c(comp5[[i]][upper.tri(comp5[[i]], diag=TRUE)])),
                  sapply(1:n, function(i) c(comp6[[i]][upper.tri(comp6[[i]], diag=TRUE)])))
    S.boot.I <- sapply(1:nrep, function(B)
      max(abs(apply(apply(Z.I, 1, function(x) r[[B]]*(x - mean(x))/sqrt(length(x))), 2, sum))))
    S.boot.II <- sapply(1:nrep, function(B)
      max(abs(apply(apply(Z.II, 1, function(x) r[[B]]*(x - mean(x))/sqrt(length(x))), 2, sum))))
  }
  c <- quantile(S.boot.I, probs = 1 - test.size)
  univ.quant <- quantile(S.boot.II + sum(abs(theta.tilde))*S.boot.I, probs = 1 - test.size)
  quant2 <- univ.quant - sum(abs(theta.tilde))*c

  return(list(C_G = c, C_W = quant2))
}
