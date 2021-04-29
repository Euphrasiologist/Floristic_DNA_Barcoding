library(MCMCglmm)

# load backup files
load(file = "../data/Backups/vcv5d.RData")
load(file = "../data/Backups/Ainv.vcv5d.RData")

#prior
prior.vcv5d <-list(R=list(V=diag(1), fix=1), 
                   G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                          G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

# model
mcmc.vcv5d <- MCMCglmm(y ~ dist2 + Hectads_shared + Annual_Perennial + Genus.Size, 
                       family = "threshold",
                       random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
                       ginverse = list(Sp1 = Ainv.vcv5d,
                                       Sp2 = Ainv.vcv5d),
                       data= vcv5d,
                       prior=prior.vcv5d,
                       pr=TRUE, 
                       nitt = 13000*100,
                       thin = 10*100,
                       burnin = 3000*100)