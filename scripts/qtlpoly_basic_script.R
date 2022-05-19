############## START of QTLPOLY #######

library(qtlpoly)

pheno<-read.csv("phenotype_file.csv",header=T, row.names = 1)
genoprob<-readRDS("genoprob_file.RDS")

data <- read_data(ploidy = 4, geno.prob = genoprob, pheno = pheno, step = 1)

#  IF you are running a new map and have not had genomewide significance estimated 
#  please run the followwing lines.
#
#   Very time intensive and after running save the data.sim, min,pvl, and score.null 
#   to load in later
# 
# data.sim <- simulate_qtl(data = data, mu = 0, h2.qtl = NULL, var.error = 1,
#                          n.sim = 1000, missing = TRUE, seed = 123)
# score.null <- null_model(data = data.sim$results, n.clusters = 10, plot = NULL)
# min.pvl <- numeric(length(score.null$results))
# for (p in 1:length(score.null$results)) {
#   min.pvl[p] <- score.null$results[[p]]$pval[which.max(score.null$results[[p]]$stat)]
# }
# 
# saveRDS(data.sim, "data.sim.RDS")
# saveRDS(min.pvl, "min.pvl.RDS")
# saveRDS(score.null, "score.null.RDS")
# 
# 
# quantile(sort(min.pvl), c(0.2, 0.05))
##          20%           5% 
## 0.0020829605 0.0003276568

#load back in the genome wide significances
data.sim<-readRDS("data.sim.RDS")
min.pvl<-readRDS("min.pvl.RDS")
score.null<-readRDS("score.null.RDS")


# this line does the QTL mapping using multiple interval mapping method if you have a powerful computer you can change n.clusters
remim.mod <- remim(data = data, w.size = 15, score.null = score.null ,sig.fwd = 0.2, sig.bwd = 0.05,   d.sint = 1.5, n.clusters = 4)

# these lines tell you your peak and upper and lower confidence intervals for your QTL
print(remim.mod)
print(remim.mod, sint = "lower")
print(remim.mod, sint = "upper")


# these lines estimate roughly how much phenotypic variation your QTL is accounting for in the h2 column
fitted.mod <- fit_model(data = data, model = remim.mod, probs = "joint", polygenes = "none")
summary(fitted.mod)


# these plot graphics of your QTL
plot_profile(data = data, model = remim.mod, grid = TRUE, sup.int = T)
plot_sint(data = data, model = remim.mod)

# these plot the graphics of the estimated effects of your QTL
est.effects <- qtl_effects(ploidy = 4, fitted = fitted.mod)
plot(est.effects,  p1="Parent 1 name", p2="Parent 2 name",pheno.col = 2)