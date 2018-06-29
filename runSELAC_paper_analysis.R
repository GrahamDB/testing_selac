
library(selac)
#source("selac.R")
tree<-read.tree('SalichosRokas.tre')
#tree.pruned <- drop.tip(tree,'Calb')

#result <- SelacOptimize(codon.data.path = 'concat.rokas/', phy = tree.pruned,
#             edge.length = 'optimize', optimal.aa = 'none', data.type='codon',
#              nuc.model = "GTR", include.gamma = FALSE, ncats = 4, numcode = 1,
#              diploid = FALSE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
#              n.cores = 60, max.restarts = 1, parallel.type='by.site')
random.seed = sample(1:100000, 100)
print(random.seed[2])
set.seed(random.seed[22])
result <- SelacOptimize(codon.data.path = '', phy = tree, n.partitions=3,
                        edge.length = 'optimize', optimal.aa = 'none', data.type='codon',
                        codon.model = 'GY94', nuc.model = 'GTR', include.gamma = FALSE, gamma.type='quadrature', ncats = 4, numcode = 1,
                        diploid = FALSE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
                        n.cores = 3, max.restarts = 1, max.evals=20)

save(result,file='yeastSalRokSelacGTRG_quad.Rdata')
