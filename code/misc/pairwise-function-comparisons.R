install.packages("pairwiseAdonis")
# install.packages("phylomatic/ecole")
install.packages("RVAideMemoire")
install.packages("vegan")

library(pairwiseAdonis)
# library(ecole)
library(RVAideMemoire)
library(vegan)

data(dune)
data(dune.env)

# PERMANOVA
vegan::adonis2(dune ~ Management, 
               data = dune.env,
               method = "euclidean")

vegan::adonis2(dist(dune, "euclidean") ~ Management, 
               data = dune.env)

# method 1 from pairwiseAdonis
pairwiseAdonis::pairwise.adonis2(dist(dune, "euclidean") ~ Management, 
                                 data = dune.env)

# method 2 from ecole
# ecole::permanova_pairwise(x = dist(dune, "euclidean"), 
#                           grp = dune.env$Management,
#                           padj = "BH")

# methods 1 and 2 have the same SumsofSquares, Fstatistic, R2, p-values
# (though no adjustment for pairwise.adonis2())

# method 3 from RVAideMemoire
# note p-values are different between this and method 2
test <- RVAideMemoire::pairwise.perm.manova(resp = dist(dune, "euclidean"),
                                    fact = dune.env$Management,
                                    p.method = "BH")

write_rds(x = test, file = "test.rds")
