############## CellFindR Functions 
# asking if the grouping is a cluster
# input: tenx = tenx object
# thresh_genes = threshold of genes at thresh_val
# thresh_val = value of the threshold in log space
# pval = cut off of pval for the signficance
is_cluster <- function(tenx, thresh_genes = 10, thresh_val = log(2), pval = 1e-4){
  val = 0 # groups that does not satisfy threshold genes
  counter = 0 # groups that satisfy threshold genes 
  # loop through the identitiy
  matrix_output <- data.frame(row.names = row.names(tenx))
  
  for (j in sort(unique(tenx@active.ident))){
    markers <- FindMarkers(tenx, ident.1 = j, min.pct = 0.25)
    markers <- markers[markers$p_val_adj < pval,]
    #find if the 10th biggest is less than log2, sum 
    print(sort(markers$avg_logFC, decreasing = TRUE)[thresh_genes])
    # if less than 10 significant genes
    if (sum(tenx@active.ident == j) < 5){
      return(FALSE)
    }
    if (length((markers$avg_logFC)) < 10){
      val <- val + 1
    } else if (sort(markers$avg_logFC, decreasing = TRUE)[thresh_genes] < thresh_val){
      #print(val)
      val <- val + 1
    } else{
      counter = counter + 1
    }
  }
  if (val > 1){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

# finds resolution that satisfy
# input: tenx object
# initial resolution of starting clustering
# how much to increment up 
# threshold of genes
# value of the threshold 
find_res <- function(tenx, initial_res = 0.1, jump = 0.1, thresh_genes = 10, thresh_val = log(2) ) {
  RES_POST <- initial_res # keeping
  RES_IT <- initial_res # iterative
  
  while(TRUE){
    print(paste('Trying',RES_IT, sep = ' '))
    tenx <- FindNeighbors(tenx, dims = 1:10)
    tenx <- FindClusters(tenx, resolution = RES_IT)
    
    # also check if theres only 1 cluster/ then can go up higher es
    # Find number of clusters
    length_group <- length(unique(tenx@active.ident))
    # if only one group then need to look deeper
    if (length_group == 1){
      print(paste('noclusterfound',RES_IT, sep = ' '))
      # still not groups at 0.7 res stop and just step as 1
      if (RES_IT == 0.7){
        print(paste('noclusterfound',RES_IT, sep = ' '))
        break
      }
    } else{
      testing <- is_cluster(tenx)
      if (testing == FALSE){ # if not real group
        print(paste('broke', RES_IT, sep = ' '))
        RES_IT <- RES_IT - jump
        RES_POST <- RES_IT
        print(RES_POST)
        break
      } else{ # valid groups
        RES_POST <- RES_IT
        print(paste('ok',RES_IT, sep = ' '))
      }
    }
    RES_IT <- RES_IT + jump
  }
  return(RES_POST)
}

# getsubclustering
# input: tenx object
# location of output folder
# project_name
sub_clustering <- function(tenx, output_folder = '.', proj_name = 'proj_name', thresh_genes = 10, thresh_val = log(2)){
  # set cell directory
  print('running through clusters to find subclusters')
  res_keep <- data.frame('cluster'= NA,'res'= NA, 'num_clusters' =NA)
  celltocluster <- data.frame(row.names = colnames(tenx))
  celltocluster$cellnames <- colnames(tenx)
  res_counter <- 1
  
  for (j in sort(unique(tenx@active.ident))){
    sub_tenx <- subset(tenx, idents = toString(j))
    print(paste('clustering ', j, sep = ''))
    
    sub_tenx <- FindVariableFeatures(sub_tenx, selection.method = "vst", nfeatures = 2000)
    sub_tenx <- FindNeighbors(sub_tenx, dims = 1:10)
    sub_tenx <- RunUMAP(sub_tenx, dims = 1:10, n.neighbors = 10)
    
    # get resolution
    set_res <- find_res(sub_tenx)
    if (set_res > 0){
      #### label everything
      sub_tenx <-FindClusters(sub_tenx,pc.use = 1:10, resolution = set_res)
    }
    
    # only one group
    num_groups <- length(unique(sub_tenx@active.ident))
    if (num_groups == 1){
      print(paste(j, ' has no subclusters', sep = ''))
      framer <- data.frame(row.names = colnames(sub_tenx))
      framer$cellnames <- row.names(framer)
      framer$clusterid <- j
      celltocluster <- merge(celltocluster, framer, all = TRUE)
    } else{
      print(paste(j, ' has ', num_groups, ' subgroups', sep = ''))
      
      # plot data
      gen_matrix_plot(sub_tenx, output_folder, j)
      
      # layer 2: X.X
      for (k in sort(unique(sub_tenx@active.ident))){
        print(paste('clustering ', j,'.',k, sep = ''))
        sub2_tenx <- subset(sub_tenx, idents = toString(k))
        
        sub2_tenx <- FindVariableFeatures(sub2_tenx, selection.method = "vst", nfeatures = 2000)
        sub2_tenx <- FindNeighbors(sub2_tenx, dims = 1:10)
        sub2_tenx <- RunUMAP(sub2_tenx, dims = 1:10, n.neighbors = 10)
        
        # get resolution
        set_res <- find_res(sub2_tenx)
        if (set_res > 0){
          #### label everything
          sub2_tenx <-FindClusters(sub2_tenx,pc.use = 1:10, resolution = set_res)
        }
        
        num_groups <- length(unique(sub2_tenx@active.ident))
        if (num_groups == 1){
          print(paste(j,'.', k, ' has no subclusters', sep = ''))
          framer <- data.frame(row.names = colnames(sub2_tenx))
          framer$cellnames <- row.names(framer)
          framer$clusterid <- paste(j,k, sep = '.')
          celltocluster <- merge(celltocluster, framer, all = TRUE)
        } else{
          print(paste(j,'.', k, ' has ', num_groups, ' subgroups', sep = ''))
          
          # plot data
          gen_matrix_plot(sub2_tenx, output_folder, paste(j, k, sep = '.'))
          
          # layer 3: X.X.X
          for (l in sort(unique(sub2_tenx@active.ident))){
            sub3_tenx <- subset(sub2_tenx, idents = toString(l))
            framer <- data.frame(row.names = colnames(sub3_tenx))
            framer$cellnames <- row.names(framer)
            framer$clusterid <- paste(j,k,l, sep = '.')
            celltocluster <- merge(celltocluster, framer, all = TRUE)
          }
        }
      }
    }
  }
  celltocluster <- na.omit(celltocluster)
  write.csv(celltocluster, paste(output_folder, "/cell_labels.csv", sep = ''), row.names = TRUE)
  celltocluster <- celltocluster[order(celltocluster$"cellnames"),]
  tenx@meta.data <- tenx@meta.data[order(rownames(tenx@meta.data)),]
  tenx@meta.data$CellfindR <- celltocluster$clusterid
  tenx <-SetIdent(tenx,cells = rownames(tenx@meta.data), value = tenx@meta.data$CellfindR)
  levels(tenx) <-str_sort(levels(tenx), numeric = TRUE)
  
  ggsave(paste(output_folder, '/', 'UMAP_CellfindR.pdf', sep = ''),    
         DimPlot(tenx, label = TRUE),width = 8, height = 8)
  return(tenx)
}

# generate matrix and plots
gen_matrix_plot <-function(tenx, output_folder = '.', proj_name = 'proj_name'){
  file_create <-paste(output_folder,'/', proj_name,sep = '')
  print(file_create)
  ### Output data destinations
  # create folder
  ggsave(paste(file_create, '_umap.pdf', sep = ''), DimPlot(tenx, label = TRUE))
  markers <-FindAllMarkers(tenx,only.pos = TRUE,min.pct = 0.25,thresh.use = 0.25)
  markers_filtered <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) 
  genes <- unique(markers_filtered$gene)
  write.csv(markers,paste(file_create, '_matrix.csv', sep = ''), row.names = TRUE)
  
  #create subdirectories
  file_create2 <- paste(output_folder, '/',proj_name, sep ='')
  dir.create(file_create2)
  dir.create(paste(file_create2, 'Cluster', sep = '/'))
  dir.create(paste(file_create2, 'Violin', sep = '/'))
  for (i in genes) {
    # cluster maps
    ggsave(paste(paste(file_create, 'Cluster', '', sep = '/'),i,'.pdf',sep = ''),
           FeaturePlot(tenx, features = c(i),pt.size = 2), width = 8, height = 8)
    #violin plot
    ggsave(paste(paste(file_create, 'Violin', '', sep = '/'), i, '.pdf', sep = ''),    
           VlnPlot(tenx, c(i)),width = 6, height = 4)
  }
}

# generate matrix
get_matrix <- function(tenx){
  print("getting matrix")
  avg_expression <- AverageExpression(tenx)
  matrix_all <- data.frame(row.names = rownames(avg_expression$RNA))
  for (i in levels(tenx@active.ident)){
    print(i)
    markers <- FindMarkers(tenx, ident.1 = i)
    avg_val <- avg_expression$RNA[i]
    avg_diff <- markers[rownames(avg_expression$RNA),]$avg_logFC
    avg_diff[is.na(avg_diff)] <-0
    p_val <- markers[rownames(avg_expression$RNA),]$p_val_adj
    p_val[is.na(p_val)] <-1
    
    matrix_all <- cbind(matrix_all, avg_val)
    matrix_all <- cbind(matrix_all, avg_diff)
    matrix_all <- cbind(matrix_all, p_val)
  }
  name_col <- c()
  for (k in levels(tenx@active.ident)){
    print(k)
    name_col <- c(name_col,(c(paste(k,'Mean', sep = '_'),paste(k,'Avg_diff', sep = '_') , paste(k,'Pval', sep = '_'))))
  }
  colnames(matrix_all) <- name_col
  matrix
  return(matrix_all)
}

# get stats
get_stats <- function(tenx, num_genes = 50){
  aoe <- c("Group", "cell_number", "avg_read", "avg_umi")
  for (i in 1:num_genes){
    aoe <- c(aoe, paste('top_',i, sep = ""))
  }
  df <- data.frame(aoe)
  #initialize matrix
  for (groups in levels(tenx@active.ident)){
    subgroup <-subset(tenx, idents = groups)
    # group name
    aod <- c(groups)
    # cell number
    aod <- c(aod, length(subgroup@meta.data$nCount_RNA))
    # avg_read
    aod <- c(aod, mean(subgroup@meta.data$nCount_RNA))
    # avg_umi
    aod <- c(aod, mean(subgroup@meta.data$nFeature_RNA))
    # top 10 diff genes
    markers <- FindMarkers(tenx, groups)
    top_markers <- row.names(markers)[1:num_genes]
    for (topm in top_markers){
      aod <- c(aod, topm)
    }
    df[groups] <-aod
  }
  return(df)
}

# getting plots
get_plots<- function(tenx, output_folder = '.'){
  dir_creater <- paste(output_folder, '/plots', sep = '')
  dir.create(dir_creater)
  for (groups in levels(tenx@active.ident)){
    markers <-FindMarkers(tenx, groups)
    maxer <- min(50, length(rownames(markers)))
    for (gene in rownames(markers)[1:maxer]){
      ggsave(paste(dir_creater,'/', gene, '_cluster.pdf', sep = ''),    
             FeaturePlot(tenx, features = gene), width = 6, height = 6)
      ggsave(paste(dir_creater,'/', gene, '_violin.pdf', sep = ''),    
             VlnPlot(tenx, features = gene, slot = "counts", log = TRUE), width = length(levels(tenx@active.ident)), height = 5)
    }
  }
}

# output metrics
metrics_output <- function(tenx, output_folder = '.', species = 'mouse'){
  tenx@active.assay
  # mito percent
  if (species == 'mouse'){
    ggsave(paste(output_folder, '/', 'percent_mito.pdf', sep = ''),    
           VlnPlot(tenx, features = 'percent.mt'), width =  length(levels(tenx@active.ident)), height = 5)
  }
  if (species == 'human'){
    ggsave(paste(output_folder, '/', 'percent_mito.pdf', sep = ''),    
           VlnPlot(tenx, features = 'percent.MT'), width =  length(levels(tenx@active.ident)), height = 5)
  }
  # umi
  print('test')
  ggsave(paste(output_folder, '/', 'uMI.pdf', sep = ''),    
         VlnPlot(tenx, features = 'nFeature_RNA'), width =length(levels(tenx@active.ident)), height = 5)
}
