
# Start ---------------------------------------------------------------------

#old_path <- Sys.getenv("PATH")
#Sys.setenv(PATH = paste(old_path, "/home/XXX/", sep = ";")) # use if files are in another folder

#import libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(LPCM)
library(foreach)

# Function definitions ----------------------------------------------------

#function to calculate distance from cell to all other cells in a set
Euclidean<-function(data, point1){
  distances<-numeric()
  data_comp<-anti_join(data, point1, by = c("POSITION_X", "POSITION_Y"))
  #calculate and return distances
  for(i in 1:nrow(data_comp)){
    distances[i]<-dist(rbind(data_comp[i,],point1))
  }
  return(distances)
}

#function to find nearest neighbor distances
nearest_neighbor_dist<-function(celldata,im_size1,im_size2){
  coord_data<-celldata
  ln<-nrow(coord_data)
  top_distances<-data.frame(firstneighbor=numeric(),secondneighbor=numeric(),thirdneighbor=numeric(),fourthneighbor=numeric())
  #search for nearest neighbor within 400 pixel square
  for(i in 1:nrow(coord_data)){
    point <- as.numeric(coord_data[i,1:2])
    high_bounds_x<-as.numeric(point[1] + 200)
    low_bounds_x<-as.numeric(point[1] - 200)
    high_bounds_y<-as.numeric(point[2] + 200)
    low_bounds_y<-as.numeric(point[2] - 200)
    search_data<-subset(coord_data, (POSITION_X > low_bounds_x) & POSITION_X < high_bounds_x &
                          (POSITION_Y > low_bounds_y & POSITION_Y < high_bounds_y))
    #if four neighbors not found, expand square
    if(nrow(search_data)<4){
      high_bounds_x<-as.numeric(point[1] + 400)
      low_bounds_x<-as.numeric(point[1] - 400)
      high_bounds_y<-as.numeric(point[2] + 400)
      low_bounds_y<-as.numeric(point[2] - 400)
      search_data<-subset(coord_data, (POSITION_X > low_bounds_x) & POSITION_X < high_bounds_x &
                            (POSITION_Y > low_bounds_y & POSITION_Y < high_bounds_y))
      if(nrow(search_data)<4){
        high_bounds_x<-as.numeric(point[1] + 600)
        low_bounds_x<-as.numeric(point[1] - 600)
        high_bounds_y<-as.numeric(point[2] + 600)
        low_bounds_y<-as.numeric(point[2] - 600)
        search_data<-subset(coord_data, (POSITION_X > low_bounds_x) & POSITION_X < high_bounds_x &
                              (POSITION_Y > low_bounds_y & POSITION_Y < high_bounds_y))
        if(nrow(search_data)<4){
          search_data<-coord_data
        }}}
    #find distance to nearest neighbors
    d<-Euclidean(search_data,point1=coord_data[i,])
    sorted_d<-sort(d)
    #remove cells on edges of image
    if(point[1]<100 | point[2]<100
       | (im_size1-point[1])<100 | (im_size2-point[2] < 100)){
      top_distances[i,1]<-NA
      top_distances[i,2]<-NA
      top_distances[i,3]<-NA
      top_distances[i,4]<-NA
    } else{
      top<-sorted_d[1:4]
      #if no neighbor found, set to maximum distance
      top[is.na(top)]<-im_size1
      top_distances[i,1:4]<-top[1:4]
    }
  }
  return(top_distances)
}

#function to calculate number of neighbors in a given boundary
Number_of_Neighbors<-function(celldata, distance_bound){
  coord_data<-celldata[c(1:2)]
  ln<-nrow(coord_data)
  neighbor_counts<-data.frame(neighbor_n=numeric(),boundary_distance=numeric())
  #define boundary box for neighbor search
  for(i in 1:nrow(coord_data)){
    point<-as.numeric(coord_data[i,1:2])
    high_bounds_x<-as.numeric(point[1] + distance_bound)
    low_bounds_x<-as.numeric(point[1] - distance_bound)
    high_bounds_y<-as.numeric(point[2] + distance_bound)
    low_bounds_y<-as.numeric(point[2] - distance_bound)
    #find the cells within the neighborhood
    search_data<-subset(coord_data, (POSITION_X > low_bounds_x) & POSITION_X < high_bounds_x &
                          (POSITION_Y > low_bounds_y & POSITION_Y < high_bounds_y))
    d<-Euclidean(search_data,point1=coord_data[i,])
    #first check if cell is on the edge, if so remove it
    if(point[1]<distance_bound | point[2]<distance_bound
       | (1215-point[1])<distance_bound | (1092-point[2] < distance_bound)){
      neighbor_counts[i,1]<-NA
      neighbor_counts[i,2]<-NA
    }
    else{
      #if no neighbors found, set to 0
      if(anyNA(d)){
        neighbor_counts[i,1]<-0
        neighbor_counts[i,2]<-distance_bound
      }else{
        cells_in_bounds<-length(d[d<distance_bound])
        neighbor_counts[i,1]<-cells_in_bounds
        neighbor_counts[i,2]<-distance_bound
      }
    }
  }
  return(neighbor_counts)
}

#function to use mean-shift clustering to calculate cluster sizes
cluster_sizes<-function(cell_data, h){
  cluster<-data.frame()
  cell_count<-nrow(cell_data)
  #for each image, check if more then one cell
  if(cell_count>1){
    coord_data<-data.matrix(cell_data[c(1:2)], rownames.force = NA)
    #perform mean-shift clustering
    ms<-ms(coord_data, scaled=0, plot=F, h=h)
    labels<-as.data.frame(ms$cluster.label) %>%
      rename(Cluster_Label=`ms$cluster.label`)
    #record the cluster labels
    label_counts<-labels %>%
      group_by(Cluster_Label) %>%
      summarize(n()) %>%
      rename(label_count='n()') %>%
      mutate(Image_Cell_Count=cell_count)
  }else{
    labels<-as.data.frame(T)
    label_counts<-labels %>%
      mutate(Image_Cell_Count=cell_count) %>%
      rename(SingleCell='T')
  }
  cluster<-bind_rows(cluster, label_counts)
  return(cluster)
}

# Image data analysis --------------------------------------------------------

### File names to read data from:
filenames <- c('B7_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'B8_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'B9_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'B10_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'B11_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'C7_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'C8_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'C9_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'C10_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'C11_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'D7_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'D8_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'D9_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'D10_-1_2_1_Stitched[imaging_TagBFP 390,447]_001',
               'D11_-1_2_1_Stitched[imaging_TagBFP 390,447]_001')
### Image sizes for each image file above (in the same order):
dims.data <- c(1208,1082, 1220,1083, 1215,1089, 1216,1088, 1218,1083,
               1215,1090, 1212,1074, 1222,1093, 1216,1098, 1208,1097,
               1216,1087, 1208,1096, 1220,1088, 1227,1081, 1204,1092)
fnamesmat <- matrix(filenames,nrow=15,byrow=TRUE)
dimsmat <- matrix(dims.data,nrow=15,byrow=TRUE)


### Use for loop for sequential run, use foreach for parallel 
# for(i in 1:nrow(fnamesmat)){ 
foreach(i=1:nrow(fnamesmat)) %do% {
  
  fname = fnamesmat[i]
  rfname = paste(fname,'.csv',sep="")
  
  imsize1 = dimsmat[i,1]
  imsize2 = dimsmat[i,2]
  
  image_dat<-read.csv(file = rfname, stringsAsFactors = F)
  cell_count_dat<-image_dat %>%
    select(POSITION_X,POSITION_Y)
  
  cell_count_dat = cell_count_dat[4:nrow(cell_count_dat),] # discard the labels
  cellsnum = nrow(cell_count_dat)
  
  cell_count_tot<-cell_count_dat
  cell_count_tot$size<-(imsize1*imsize2) # image size, differs a little for each
  cell_count_tot$null<-.5/sqrt(nrow(cell_count_dat)/cell_count_tot$size)
  
  
  ###find distances to nearest neighbors
  total_cell_distances<-cell_count_tot[1:cellsnum,] %>%
    do(data.frame(., e=nearest_neighbor_dist(.,imsize1,imsize2)))
  #bind to cell data
  IF_cells_distance<-(full_join(cell_count_tot[1:cellsnum,], total_cell_distances)) %>%
    rename(FirstNeighbor_Dist=e.firstneighbor,SecondNeighbor_Dist=e.secondneighbor,
           ThirdNeighbor_Dist=e.thirdneighbor,FourthNeighbor_Dist=e.fourthneighbor)
  #calculate normalized distance metrics
  IF_cells_distance$Normalized_First_Neighbor_Dist<-IF_cells_distance$FirstNeighbor_Dist/IF_cells_distance$null
  IF_cells_distance$Normalized_Second_Neighbor_Dist<-IF_cells_distance$SecondNeighbor_Dist/IF_cells_distance$null
  IF_cells_distance$Normalized_Third_Neighbor_Dist<-IF_cells_distance$ThirdNeighbor_Dist/IF_cells_distance$null
  IF_cells_distance$Normalized_Fourth_Neighbor_Dist<-IF_cells_distance$FourthNeighbor_Dist/IF_cells_distance$null
  
  
  ###find nearest neighbors
  total_neighborhood_numbers<-cell_count_tot[1:cellsnum,] %>%
    do(data.frame(., e=Number_of_Neighbors(., distance_bound=100)))
  #rename neighbor measure
  IF_cells_neighbors<-(full_join(cell_count_tot[1:cellsnum,], total_neighborhood_numbers)) %>%
    rename(number_neighbors=e.neighbor_n)
 
  
  ###bind all cell data together
  total_distance_data<-full_join(IF_cells_distance, IF_cells_neighbors)
  #output distance csv
  rfname = paste(fname,'_Cell_File.csv',sep="")
  write_csv(total_distance_data, rfname)
  
  
  ###mean shift clustering analysis
  total_cluster<-cluster_sizes(cell_data=cell_count_tot[1:cellsnum,], h=30)
  #calculate proportion of cells in cluster > n=8
  image_cluster_proportions<-total_cluster %>%
    filter(label_count>8) %>%
    mutate(cells_in_large_cluster=sum(label_count))
  cluster_bind<-full_join(image_cluster_proportions, total_cluster)
  cluster_bind[is.na(cluster_bind$cells_in_large_cluster),]$cells_in_large_cluster<-0
  #calculate image means for all measures
  image_props<-cluster_bind %>%
    mutate(proportion_in_large_cluster=max(cells_in_large_cluster)/mean(Image_Cell_Count)) %>%
    summarise(proportion_in_large_cluster=mean(proportion_in_large_cluster, na.rm = T))
  #output cluster csv
  rfname = paste(fname,'_Cluster_File.csv',sep="")
  write_csv(cluster_bind, rfname)
  
  
  ### calculate summaries/means for the metrics
  dist_means<-total_distance_data %>%
    summarize_if(is.numeric, funs(mean, .args = list(na.rm=T)))
  clust_means<-cluster_bind %>%
    summarize_if(is.numeric, funs(mean, .args = list(na.rm=T)))
  #output summaries
  write_csv(dist_means, paste(fname,'_Dist_means.csv',sep=""))
  write_csv(clust_means, paste(fname,'_Clust_means.csv',sep=""))
  write_csv(image_props, paste(fname,'_image_props.csv',sep=""))
}