from phylopandas import read_newick
from phylovega import VegaTree
import dendropy

pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open("bacterias_cambio_256.dat"), delimiter=",")

nj_tree = pdm.nj_tree()
print nj_tree.as_string("newick")
