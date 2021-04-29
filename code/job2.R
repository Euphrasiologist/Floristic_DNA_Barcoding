library(MCMCglmm)

# load backup files
load(file = "../data/Backups/vcv6d.RData")
load(file = "../data/Backups/Ainv.ploidd.RData")

#prior
prior.h1d <-list(R=list(V=diag(1), fix=1), 
                 G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                        G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
# model
mcmc.ploidd <- MCMCglmm(y ~ dist2 + Cross_Ploid + Hectads_shared + Annual_Perennial + Genus.Size, 
                        family = "threshold",
                        random = ~ mm(Sp1 + Sp2) + mm(T1 + T2),
                        ginverse = list(Sp1 = Ainv.ploid,
                                        Sp2 = Ainv.ploid),
                        data= vcv6d,
                        prior=prior.h1d,
                        pr=TRUE, 
                        nitt = 13000*100,
                        thin = 10*100,
                        burnin = 3000*100)