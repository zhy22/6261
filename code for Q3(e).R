library(ggplot2)

logsumexp=function(loga, logb) {
  if (is.infinite(loga) && is.infinite(logb)) {
    return(-Inf)
  }
  max.log=max(loga, logb)
  return(max.log+log(exp(loga-max.log)+exp(logb-max.log)))
}

nebc=function(N,delta, ph=0.75) {
  pt=1-ph
  results=data.frame(delta=double(length(delta)), normH=double(length(delta)))
  for (i in seq_along(delta)) {
    d=delta[i]
    logcumprob=-Inf  
    logcumseq=-Inf
    for (k in seq(N, 0, by=-1)) {
      logpk=k*log(ph)+(N-k)*log(pt)
      logCNK=lgamma(N+1)-lgamma(k+1)-lgamma(N-k+1)
      logtpk=logCNK+logpk
      logcumprob=logsumexp(logcumprob, logtpk)
      logcumseq=logsumexp(logcumseq, logCNK)
      cumprob=exp(logcumprob)
      if (cumprob >= 1-d) {
        H.delta=logcumseq/log(2)
        normH=H.delta/N
        results$delta[i]<-d
        results$normH[i]<-normH
        break
      }
    }
  }
  return(results)
}

Nvalues=c(10, 100, 500, 1000, 1800)
delta=seq(0.0001, 0.99, length.out=100)

all.results=data.frame()

for (N in Nvalues) {
  results=nebc(N,delta)
  results$N=as.factor(N)
  all.results=rbind(all.results,results)
}

ggplot(all.results, aes(x=delta, y=normH, color=N)) +
  geom_line(size=1) +
  geom_hline(yintercept=0.811, linetype="dashed", color="red", size=0.8) +
  annotate("text", x=0.02, y=0.78, label="H(X) = 0.811",color="red",hjust = 0) +
  labs(x='Delta', y='Normalized Essential Bit Content', title='Normalized Essential Bit Content - δ') +
  theme_minimal() +
  theme(
    legend.position.inside=c(0.1, 0.1),  
    legend.justification=c(0, 0),    
    legend.background=element_rect(fill="white", color=NA),  
    legend.title=element_text(size=10),
    legend.text=element_text(size=9)
  )


