# Introgression Clustering
This clustering algorithm differentiates genetic histories of incomplete lineage sorting and introgression with whole genome phylogenetic data.

Build a distance matrix for locus trees across a chromosome.
```
data <- get_distances(trees_real,minor,major,hybrid,outgroup,out)
```

Visualize first 3 dimensions of PCA for distance matrix.
```
PCAdistancematrix_allcombinations(data$reslist[[1]], nf = 3, out, trees_real, minor[1], major[1], hybrid[1])
```

Calculate clusters of locus-tree topologies (excluding branch length information). Evaluate lieklihood of each cluster to contain introgression histories. 
```
subsets <-check_subsets(trees_real,outgroup,data$branches_clusters,data$uni,data$roundedpcoli,data$pco_related,data$combinations,data$reslist,0,data$df,catTree,x,y,z)
```

Using the community with the most likely introgression histories, identify  additional members of that cluster using the silhouette score. Locus trees with score larger than 0 are returned as introgression histories.
```
chosen_community <- community(-1,subsets$branches_clusters,trees_real,data$uni,data$roundedpcoli,data$combinations,data$reslist,data$pco_related,paste0(out),x,y,z))
```

