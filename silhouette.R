#Silhouette Functions

#Load packages###########
suppressMessages(library(rgl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(distory))
suppressMessages(library(plotly))
suppressMessages(library(devtools))
suppressMessages(library(htmltools))
suppressMessages(library("ggplot2"))
suppressMessages(library("ape"))
suppressMessages(library("hash"))
suppressMessages(library("ggtree"))
suppressMessages(library("ggpubr"))
suppressMessages(library("dplyr"))
install_github("thibautjombart/treespace")
suppressMessages(library(treespace))

#tree_identity- given a minor, major, and hybrid parent, it identifies if a set of trees match introgression, ILS, or species tree topology
tree_identity <- function(trees,minor,major,hybrid) {
  treeidentity <- data.frame(id = rep(0, length(trees)))
  for(i in 1:length(trees)){
    tree_sub <- keep.tip(trees[[i]],c(minor,major,hybrid))
    if(is.monophyletic(tree_sub, c(minor,hybrid))){
      treeidentity$id[i] = "INT"
    }else if(is.monophyletic(tree_sub, c(major,hybrid))){
      treeidentity$id[i] = "SP"
    }else{
      treeidentity$id[i] = "ILS"
    }
  }
  trees_sub <- keep.tip(trees,c(minor,major,hybrid))
  treeidentity$branchlengths=NA
  for(i in 1:length(trees_sub)){
    treeidentity$branchlengths[i] = dist.nodes(trees_sub[[i]])[4,5]
  }
  return(treeidentity)
}

#Add information about the real trees
tree_identity_simulation <- function(windows,ids) {
  #add line to ids to include $realid
  return(treeidentity)
}

#Gets clusters of trees based on their related distance to the introgression tree
related_dist <- function(trees,minor,major,hybrid, outgroup,out) {
  df<- cbind(c(rep("Minor",length(minor)),rep("Major",length(major)),rep("Hybrid",length(hybrid)),rep("Out",1)),
             c(minor,major,hybrid,outgroup))
  x <- relatedTreeDist(trees,df)
  pco <- dudi.pco(x, scannf =F,nf = 10)
  ids <- tree_identity(trees_real,minor[1],major[1],hybrid[1])
  roundedpcoli <- round(pco$li,5)
  uni <- unique(roundedpcoli[which(ids$id == "INT"),])
  
  fig <- plot_ly(
    x = pco$li[, 1],
    y = pco$li[, 2],
    z = pco$li[, 3],
    color = ids$id,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    layout(title = list(text = "Related Tree Dists", font = list(size = 18)  ))
  htmlwidgets::saveWidget(fig, paste0(out, "_relatedtreedistance.html"), selfcontained = TRUE)
  
  return(list(uni=uni,roundedpcoli=roundedpcoli,ids=ids, fig=fig,df=df,relateddist=x,pco=pco))
}

evaluate_subsets <- function(uni,roundedpcoli,ids, out, simulations){
  ggs <- list()
  if(simulations==1){
    branches_clusters = data.frame(unival = c(1:dim(uni)[1]),numgenes=rep(0,dim(uni)[1]),meanbl=rep(0,dim(uni)[1]), avgINT = rep(0,dim(uni)[1]))
  }else{
    branches_clusters = data.frame(unival = c(1:dim(uni)[1]),numgenes=rep(0,dim(uni)[1]),meanbl=rep(0,dim(uni)[1]))
  }
  for(i in 1:dim(uni)[1]){
    vals = which(roundedpcoli$A1 == uni[i,1])
    hist_data <- ids$branchlengths[vals]
    ggs[[i]] <- ggplot(data.frame(branchlengths = hist_data), aes(x = branchlengths)) +
      geom_histogram()+
      ggtitle(paste("Community", i,as.character(uni[i, 1]))) +
      theme_minimal()
    branches_clusters$numgenes[i] = length(ids$branchlengths[vals])
    branches_clusters$meanbl[i] = mean(ids$branchlengths[vals])
    if(simulations == 1){
      branches_clusters$avgINT[i] = mean(ids$realid)
    }
  }
  
  ggs[[dim(uni)[1]+1]] <- ggplot(data=ids, aes(x = branchlengths)) +
    geom_histogram() +
    ggtitle("Branch Length Distribution")
  theme_minimal()
  
  pdf(paste0(out,"_preliminaryhistograms.pdf"))
  for (i in seq_along(ggs)) {
    print(ggs[[i]])
  }
  dev.off()
  return(list(branches_clusters=branches_clusters,ggs=ggs))
}

triplet_distances <- function(trees,combinations,outgroup){
  reslist <- list()
  for(i in 1:dim(combinations)[1]){
    hyb = as.character(combinations[i,"hybrid"])
    maj =as.character(combinations[i,"major"])
    min = as.character(combinations[i,"minor"])
    trees_sub <- keep.tip(trees, c(hyb,maj,min,outgroup))
    D <- dist.multiPhylo(trees_sub)
    if (!ade4::is.euclid(D)) {
      warning("Distance matrix is not Euclidean; making it Euclidean using ade4::cailliez")
      D <- ade4::cailliez(D, print = FALSE)
    }
    reslist[[i]] <- D
  }
  return(reslist)
}

#1.For triplet sets of taxa, gets distance matricies
#2. evaluates the gene tree subset communities
get_distances <- function(trees,minor,major,hybrid,outgroup,out) {
  plotlys <- list()
  numplotlys=1
  
  simulations=0
  rel <- related_dist(trees,minor,major,hybrid, outgroup,out)
  branch <- evaluate_subsets(rel$uni,rel$roundedpcoli,rel$ids, out,simulations)
  combinations <- expand.grid(major= major,minor=minor,hybrid=hybrid)
  reslist <- triplet_distances(trees,combinations,outgroup)
  
  #Saving data
  data <- list(df = rel$df,relateddist = rel$relateddist,pco_related = rel$pco,roundedpcoli = rel$roundedpcoli,uni = rel$uni,branches_clusters = branch$branches_clusters,plotlys = rel$fig,ggs = branch$ggs,combinations = combinations,
               reslist = reslist)
  save(data,file=paste0(out,"_distancedata.Rdata"))
  return(data)
}

#Plot 3D colored by an ID
create_3d_plot <- function(data, x_col, y_col, z_col, ids) {
  plot_ly(
    x = data[, x_col],
    y = data[, y_col],
    z = data[, z_col],
    color = ids$id,
    type = "scatter3d",
    mode = "markers",
    text = paste("Index:", seq_len(nrow(data))),  # Add index as hover text
    hoverinfo = "text+color",
    marker = list(size = 3)
  ) %>%
    layout(
      title = paste("3D Plot: Columns", x_col, y_col, z_col),
      scene = list(
        xaxis = list(title = colnames(data)[x_col]),
        yaxis = list(title = colnames(data)[y_col]),
        zaxis = list(title = colnames(data)[z_col])
      )
    )
}

#If you want to investigate what the principal components of treespace look like then you can 
#give this function one of the distance matricies and it will plot all combinations of the 
#first 10 eigenvectors.
PCAdistancematrix_allcombinations <- function(D, nf = 3, out = "output", trees_real, min, maj, hyb) {
  pco_1 <- dudi.pco(D, scannf = FALSE, nf = nf)
  column_combinations <- combn(1:10, 3)
  
  # Create a list of plots
  plots <- lapply(seq_len(ncol(column_combinations)), function(i) {
    cols <- column_combinations[, i]
    create_3d_plot(pco_1$tab, cols[1], cols[2], cols[3], tree_identity(trees_real, min, maj, hyb))
  })
  
  # Wrap plots in a tagList for HTML output
  html_page <- tagList(plots)
  
  # Save to HTML
  save_html(html_page, file = paste0(out, "_PCAaxes.html"))
}


average_distance_matrix <- function(reslist, combinations, trees_real,ids){
  distmat <- as.matrix(reslist[[1]])
  
  #Get average distance matrix
  if((dim(combinations)[1])>1){
    for(i in 2:dim(combinations)[1]){
      hyb = as.character(combinations[i,"hybrid"])
      maj =as.character(combinations[i,"major"])
      min = as.character(combinations[i,"minor"])
      ids_temp <- tree_identity(trees_real,min,maj,hyb)
      ids$id = ifelse(ids_temp$id == "INT", "INT",ids$id)
    }
    for(b in 2:length(reslist)){
      res <- reslist[[b]]
      distmat <- as.matrix(res) + distmat
    }
  }
  distmat <- distmat/length(reslist)
}

plot_subset <- function(g,sil,x,y,z,pco,sil_allcombos_justint){
  color_range <- quantile(sil$avg, probs = c(0.05, 0.95),na.rm=T)  # 5th and 95th percentiles
  p <- plot_ly(
    x = pco$tab[, x],
    y = pco$tab[, y],
    z = pco$tab[, z],
    color = sil$avg,
    type = "scatter3d",
    mode = "markers",
    text = paste("Index:", seq_len(nrow(pco$li))),  # Add index as hover text
    hoverinfo = "text+color",
    marker = list(
      size = 3,
      colorbar = list(title = "Silhouette Avg"),  # Optional: Add colorbar title
      cmin = color_range[1],  # Minimum for color mapping
      cmax = color_range[2]   # Maximum for color mapping
    )
  )
  q = ggplot(data=sil) + geom_histogram(aes(x=avg))+ggtitle(g) + theme_bw()
  pp <- ggplot() + geom_histogram(data=sil_allcombos_justint[which(sil_allcombos_justint$avg < 0),],aes(x=branchlengths),fill="red",alpha=0.5)+ geom_histogram(data=sil_allcombos_justint[which(sil_allcombos_justint$avg > 0),],aes(x=branchlengths),fill="blue",alpha=0.5) +theme_bw()+ ggtitle(g)
  return(list(p,q,pp))
}

plot_multiple_arrangements <- function(plots,out,chunk_size){
  num_chunks <- ceiling(length(plots) / chunk_size)

  for (chunk in seq_len(num_chunks)) {
    chunk_plots <- plots[((chunk - 1) * chunk_size + 1):(chunk * chunk_size)]
    chunk_file <- tagList(
      tags$head(
        tags$style(HTML("
        .tab { display: none; }
        .tab.active { display: block; }
        .tab-button { cursor: pointer; margin-right: 10px; background-color: #f1f1f1; padding: 5px; border: 1px solid #ccc; }
        .tab-button.active { background-color: #ddd; }
      "))
      ),
      tags$body(
        h1(paste("Toggleable 3D Scatter Plots - Chunk", chunk)),
        div(
          id = "tabs",
          do.call(
            tagList,
            lapply(seq_along(chunk_plots), function(i) {
              tags$button(
                class = "tab-button",
                `data-tab` = paste0("tab-", i),
                paste("Plot", i)
              )
            })
          )
        ),
        div(
          id = "tab-content",
          do.call(
            tagList,
            lapply(seq_along(chunk_plots), function(i) {
              div(
                class = "tab",
                id = paste0("tab-", i),
                chunk_plots[[i]]
              )
            })
          )
        ),
        tags$script(HTML("
        document.querySelectorAll('.tab-button').forEach(button => {
          button.addEventListener('click', () => {
            document.querySelectorAll('.tab').forEach(tab => tab.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(btn => btn.classList.remove('active'));
            button.classList.add('active');
            document.getElementById(button.dataset.tab).classList.add('active');
          });
        });
        document.querySelector('.tab-button').click(); // Activate the first tab by default
      "))
      )
    )
    save_html(chunk_file, file = paste0(out, "_chunk_", chunk, ".html"))
  }
}

#calculate the silhouette score
get_sil <- function(D,identit1,int_outliers,sim){
  distmat <- as.matrix(D)
  intindices = which(identit1[,1] == "INT") #ils or int
  ilsindices <- which(identit1[, 1] == "ILS")
  spindices <- which(identit1[, 1] == "SP")
  mean_distances_ILS <- numeric(length(intindices))
  median_distances_ILS <- numeric(length(intindices))
  mean_distances_SP <- numeric(length(intindices))
  median_distances_SP <- numeric(length(intindices))
  mean_distances_all <- numeric(length(intindices))
  median_distances_all <- numeric(length(intindices))
  mean_distances_outlier <- numeric(length(intindices))
  median_distances_outlier <- numeric(length(intindices))
  if(is.character(int_outliers)){
    int_outliers=intindices
  }
  for(i in seq_along(intindices)){
    distances_to_ils <- distmat[intindices[i],ilsindices]
    mean_distances_ILS[i] <- mean(distances_to_ils)
    median_distances_ILS[i] <- median(distances_to_ils)
    distances_to_sp <- distmat[intindices[i],spindices]
    mean_distances_SP[i] <- mean(distances_to_sp)
    median_distances_SP[i] <- median(distances_to_sp)

    distances_to_all <- distmat[intindices[i],c(spindices,ilsindices)]
    mean_distances_all[i] <- mean(distances_to_all)
    median_distances_all[i] <- median(distances_to_all)

    if(i %in% int_outliers){y= int_outliers[int_outliers != i]}else{y=int_outliers}
    distances_to_outlier <- distmat[intindices[i],y]
    mean_distances_outlier[i] <- mean(distances_to_outlier)
    median_distances_outlier[i] <- median(distances_to_outlier)
  }
  if(sim !="sim"){
    distances_between_ilspoints_intpoint <- data.frame(
      point_index = intindices,
      mean_distance_ils = mean_distances_ILS,
      median_distance_ils = median_distances_ILS,
      mean_distance_sp = mean_distances_SP,
      median_distance_sp = median_distances_SP,
      mean_distance_all = mean_distances_all,
      median_distance_all = median_distances_all,
      min_mean = pmin(mean_distances_ILS,mean_distances_SP),
      min_median = pmin(median_distances_ILS,median_distances_SP),
      mean_distance_outlier = mean_distances_outlier,
      median_dist_outlier = median_distances_outlier
      # pint=identit1$propintro[intindices],
      # pils=identit1$propils[intindices]
    )
  }else{
    distances_between_ilspoints_intpoint <- data.frame(
      point_index = intindices,
      mean_distance_ils = mean_distances_ILS,
      median_distance_ils = median_distances_ILS,
      mean_distance_sp = mean_distances_SP,
      median_distance_sp = median_distances_SP,
      mean_distance_all = mean_distances_all,
      median_distance_all = median_distances_all,
      min_mean = pmin(mean_distances_ILS,mean_distances_SP),
      min_median = pmin(median_distances_ILS,median_distances_SP),
      mean_distance_outlier = mean_distances_outlier,
      median_dist_outlier = median_distances_outlier,
      pint=identit1$propintro[intindices],
      pils=identit1$propils[intindices]
    )
  }
  distances_between_ilspoints_intpoint$sil= (distances_between_ilspoints_intpoint$min_mean - distances_between_ilspoints_intpoint$mean_distance_outlier)/pmax(distances_between_ilspoints_intpoint$min_mean, distances_between_ilspoints_intpoint$mean_distance_outlier)
  return(distances_between_ilspoints_intpoint)
}

#Choose the subset community with the most likely introgression histories
check_subsets <- function(trees_real,outgroup,branches_clusters,uni,roundedpcoli,pco,combinations,reslist,drop,df,catTree,xax,yax,zax){
  p <- list();pp <- list();q <- list()
  #get concordance and minimum sil score
  branches_clusters$meansil = NA
  branches_clusters$minsil = NA
  branches_clusters$concordmean=NA
  branches_clusters$ilsblmean=NA
  branches_clusters$intblmean=NA
  hyb = as.character(combinations[1,"hybrid"])
  maj =as.character(combinations[1,"major"])
  min = as.character(combinations[1,"minor"])
  ids <- tree_identity(trees_real,min,maj,hyb)
  ids$branchlengths=NA
  trees_sub <- keep.tip(trees_real, c(hyb,maj,min,outgroup))
  for(j in 1:length(trees_sub)){
    if(is.monophyletic(trees_sub[[j]],c(hyb,min))){
      ids$branchlengths[j] = dist.nodes(trees_sub[[j]])[4,5]
    }
  }

  distmat <- average_distance_matrix(reslist, combinations, trees_real,ids)
  
  #Assess silhouette score per group
  for(g in 1:dim(uni)[1]){ #identify the silhouette score of those points to the rest
    if(branches_clusters$numgenes[g] != 1){
      sil <- data.frame(tree=c(1:length(trees_real)))
      trees2 <- keep.tip(trees_real[which(roundedpcoli$A1 == uni[g,1])],df[,2])
      #print(trees2[[1]]$tip.label)
      if(drop != 0){    indTree1 <- drop.tip(trees2,drop) }else{indTree1 <-trees2}
      concordances <- sapply(indTree1, function(x) treeConcordance(catTree,x,df))
      branches_clusters$concordmean[g] = mean(concordances)
      
      ys_simall = which(roundedpcoli$A1 == uni[g,1])
      #Get score of affinity to own community or external community
      data <- data.frame(tree=ys_simall, affinitytoINT=NA, affinitytonearest=NA,nearest=NA,sil=NA)
      for(i in seq_along(ys_simall)){
        #affinity to community
        y=which(ids[, 1] == "INT") #ys_simall[-(i)]
        data$affinitytoINT[i] = mean(distmat[ys_simall[i],y])
        #affinity to other communities
        ilsindices <- which(ids[, 1] == "ILS")
        spindices <- which(ids[, 1] == "SP")
        c <- mean(distmat[ys_simall[i],spindices])
        d <- mean(distmat[ys_simall[i],ilsindices])
        if(min(c,d) == c){ data$nearest[i] = "SP" }else if(min(c,d) == d){ data$nearest[i] = "ILS" }else{ data$nearest[i] = -1 }
        data$sil[i] <- (min(c,d) - data$affinitytoINT[i])/max(min(c,d),data$affinitytoINT[i])
        
      }
      
      branches_clusters$meansil[g] = mean(data$sil)
      branches_clusters$minsil[g] = min(data$sil)
      
      sil_dylan <- get_sil(distmat,ids,ys_simall,0)
      sil_dylan$branchlengths=ids$branchlengths[which(ids$id == "INT")]
      sil_dylan_2 <- sil_dylan[,c("point_index","sil","branchlengths")]
      colnames(sil_dylan_2)[2] = "avg"
      sil <- merge(sil,sil_dylan_2,by.x = "tree",by.y = "point_index",all.x = TRUE)
      #print(sil)
      sil_allcombos_justint = sil[which(!is.na(sil$avg)),]
      plots <- plot_subset(g,sil,xax,yax,zax,pco,sil_allcombos_justint)
      p[[g]] <- plots[[1]];q[[g]] <- plots[[2]];pp[[g]] <- plots[[3]]
      branches_clusters$ilsblmean[g] = mean(sil_allcombos_justint[which(sil_allcombos_justint$avg < 0),"branchlengths"],na.rm=T)
      branches_clusters$intblmean[g] = mean(sil_allcombos_justint[which(sil_allcombos_justint$avg > 0),"branchlengths"],na.rm=T)
    }
  }
  pdf(paste0(out,"_gs_bl_histograms.pdf"),width=5,height=3)
  for (i in seq_along(pp)) {print(pp[[i]])}
  for (i in seq_along(q)) {print(q[[i]])}
  dev.off()
  
  plot_multiple_arrangements(p,paste0(out,"_gs_plotly_dists.html"),15)
  branches_clusters$difference=branches_clusters$intblmean-branches_clusters$ilsblmean
  branches_clusters$intvils = ifelse(branches_clusters$intblmean>branches_clusters$ilsblmean,1,0)
  branches_clusters$combined = branches_clusters$intvils*branches_clusters$concordmean * branches_clusters$minsil
  branches <- list(branches_clusters=branches_clusters,sil = sil,sil_allcombos_justint = sil_allcombos_justint,g_bl_hist = pp,g_plotly = p,sil_hist = q)
  write.csv(branches_clusters, paste0(out,"_subsetdata.csv"),row.names = F)
  save(branches,file=paste0(out,"_subsetdata.Rdata"))
  return(branches)
}

#given a subset, identify other introgression community members
community <- function(g,branchdata,trees_real,uni,roundedpcoli,combinations,reslist,pco, out,x,y,z){
  if(g== -1){
    g=branchdata$unival[which(branchdata$combined == max(branchdata$combined,na.rm = T))]
    print(paste0("Best community value: ",g))
  }
  ys = which(roundedpcoli$A1 %in% uni[g,1])
  sil <- data.frame(tree=c(1:length(trees_real)))
  sil_stats <- data.frame(tree=c(1:length(trees_real)))
  if(dim(combinations)[1] == 1){
    i=1
    hyb = as.character(combinations[i,"hybrid"])
    maj =as.character(combinations[i,"major"])
    min = as.character(combinations[i,"minor"])
    ids <- tree_identity(trees_real,min,maj,hyb)
    res <- reslist[[i]]
    
    sil_dylan <- get_sil(res,ids,ys,0)
    sil_dylan$branchlengths=ids$branchlengths[which(ids$id == "INT")]
    sil_stats[,paste(hyb,maj,min,"ID")] = ids[,1]
    sil_stats[,paste(hyb,maj,min,"bl")] = ids[,2]
    sil_dylan_2 <- sil_dylan[,c("point_index","sil")]
    colnames(sil_dylan_2)[2] = paste0(hyb, maj, min)
    sil <- merge(
      sil,
      sil_dylan_2,
      by.x = "tree",
      by.y = "point_index",
      all.x = TRUE
    )
    sil$avg = (sil[,2])
    
  }else{
    for(i in 1:dim(combinations)[1]){
      hyb = as.character(combinations[i,"hybrid"])
      maj =as.character(combinations[i,"major"])
      min = as.character(combinations[i,"minor"])
      ids <- tree_identity(trees_real,min,maj,hyb)
      res <- reslist[[i]]
      
      sil_dylan <- get_sil(res,ids,ys,0)
      sil_dylan$branchlengths=ids$branchlengths[which(ids$id == "INT")]
      sil_stats[,paste(hyb,maj,min,"ID")] = ids[,1]
      sil_stats[,paste(hyb,maj,min,"bl")] = ids[,2]
      sil_dylan_2 <- sil_dylan[,c("point_index","sil")]
      colnames(sil_dylan_2)[2] = paste0(hyb, maj, min)
      sil <- merge(
        sil,
        sil_dylan_2,
        by.x = "tree",
        by.y = "point_index",
        all.x = TRUE
      )
    }
    #see how they all correlate with each other?
    sil$avg = rowMeans(sil[,2:dim(sil)[2]],na.rm=T)
  }
  sil$branchlength_temp = sil_stats[,3]
  sil_allcombos_justint = sil[which(!is.na(sil$avg)),]
 # head(sil)
#  head(sil_allcombos_justint)
  plots <- plot_subset(g,sil,x,y,z,pco,sil_allcombos_justint)
#  p <- plots[[1]];q <- plots[[2]];pp <- plots[[3]]
  write.table(sil_allcombos_justint, paste0(out,"_INTtrees.txt"),quote=F, row.names = F)
  community <- list(sil = sil,sil_allcombos_justint = sil_allcombos_justint,g_bl_hist = plots[[3]],g_plotly = plots[[1]],sil_avg = plots[[2]])
  save(community,file=paste0(out,"_communitydata.Rdata"))
  return(community)
}




args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
name2 <- args[2]
minor <- strsplit(args[3], ",")[[1]]  # If passed as a comma-separated string
major <- strsplit(args[4], ",")[[1]]  # If passed as a comma-separated string
hybrid <- strsplit(args[5], ",")[[1]]  # If passed as a comma-separated string
start <- args[6]
bestcomm <- args[7]
# You can then print or use these variables
print(out)
print(name2)
print(minor)
print(major)
print(hybrid)
print(start)
catTree <- read.tree(text="(Out,(Major,(Minor,Hybrid)));")

x=1;y=2;z=3
outgroup="Lampropeltis-triangulum-SR1453"
trees_real <- root(read.tree(paste0(name2,".treefile")),outgroup,resolve.root=T)
if(start == 0){
 	system.time(data <- get_distances(trees_real,minor,major,hybrid,outgroup,out))
 	system.time(PCAdistancematrix_allcombinations(data$reslist[[1]], nf = 3, out, trees_real, minor[1], major[1], hybrid[1]))
	system.time(branches <- check_subsets(trees_real,outgroup,data$branches_clusters,data$uni,data$roundedpcoli,data$pco_related,data$combinations,data$reslist,0,data$df,catTree,x,y,z))
	#Examine subsetdata.csv file to chose the best community OR continue with the best community
	system.time(community1 <- community(-1,branches$branches_clusters,trees_real,data$uni,data$roundedpcoli,data$combinations,data$reslist,data$pco_related,paste0(out),x,y,z))
}else if(start > 0){
  load(paste0(out,"_distancedata.Rdata"))
  load(paste0(out,"_subsetdata.Rdata"))
  community1 <- community(bestcomm,trees_real,data$uni,data$roundedpcoli,data$combinations,data$reslist,data$pco_bhv,paste0(out),x,y,z)
}
