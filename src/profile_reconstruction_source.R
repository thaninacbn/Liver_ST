library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(readxl)
library(tidyverse)
library(reshape2)
library(ggvenn)
library(mclust)
library(zoo)
library(dplyr)

# Function that shows the profile of expression for one region ----------------

SingleRegion <- function(sample, 
                         starting_point,
                         dist = 400,
                         GOI,
                         trig = F){
  
  
  sample@meta.data$imagerow <- sample@images$slice1@coordinates[ , "imagerow" ]
  sample@meta.data$imagecol <- sample@images$slice1@coordinates[ , "imagecol" ]
  df_coords =  sample@images$slice1@coordinates[starting_point,]
  row = df_coords$row/(1/sqrt(0.75)) # rescaling  so one dimension unit == 100 microns
  col = df_coords$col/2
  imrow = df_coords$imagerow
  imcol = df_coords$imagecol
  dx = as.integer(dist/2.4)
  dy = as.integer(dist/3.6)
  
  # Creates a square around the starting point proportional to the desired distance
  Region = sample[,which(sample$imagerow < (imrow + dx))]
  Region = Region[,which(Region$imagerow > (imrow - dx))]
  Region = Region[,which(Region$imagecol < (imcol + dx))]
  Region = Region[,which(Region$imagecol > (imcol - dx))]
  
  
  # rescaling  so one dimension unit == 100 microns
  Region_rescaled = as.data.frame(Region@images$slice1@coordinates[,2:3])
  pt = c(row,col)
  Region_rescaled$col = Region_rescaled$col/2 
  Region_rescaled$row = Region_rescaled$row/(1/sqrt(0.75))
  
  # Calculate Euclidean distances of pt to all points in df
  dist_to_pt <- as.matrix(dist(rbind(Region_rescaled, pt)))[nrow(Region_rescaled) + 1, 1:nrow(Region_rescaled)]
  dist_to_pt = round(dist_to_pt, 3)
  
  dist_list = seq(0,(dist/100 + 1),by=1) + .00001
  #prepare results dataframe
  res <- data.frame(matrix(ncol = 4, nrow = 0))
  resnames <- c("GOI", "barcode", "expression", "dist")
  colnames(res) <- resnames
  
  
  for (i in 1:(length(dist_list)-1)) {
    a = dist_list[i]
    b = dist_list[i+1]
    
    spot_list <- names(dist_to_pt[dist_to_pt>= a & dist_to_pt<=b])
    #fetch data in expression matrix
    expr = as.data.frame(rowMeans(as.data.frame(FetchData(Region, GOI, cells = spot_list, layer = "data")))) 
    colnames(expr) = "expression"
    
    # Fill results dataframe
    generep <- rep(GOI, nrow(expr))
    distrep <- rep(i, nrow(expr))
    tmp_df <- data.frame(matrix(nrow = nrow(expr)))
    tmp_df$GOI <- generep
    tmp_df$barcode <- rownames(expr)
    tmp_df$expression <- expr$expression
    tmp_df$dist <- distrep
    res <- rbind(res, tmp_df)
  }
  
  rstart = names(dist_to_pt[dist_to_pt>= 0.01 & dist_to_pt<=1.01])
  rmax = spot_list

  
  res <- res[, 2:5]
  
  # Create a distance column with the actual distances in microns
  res$distance = as.double((res$dist - 1)*100)
  
  #Create plot of the distance around the spot
  spat <- SpatialDimPlot(Region, 
                         crop=TRUE,
                         cells.highlight = c(rstart, rmax),
                         pt.size.factor = 2.5) + labs(fill = "Starting Point")
 
  
  
  p <- ggplot(res, aes(x = distance, y = expression)) +
    geom_smooth(method = "loess") +
    theme_minimal()
  
  
  if (trig ==FALSE) {
    return(p)
  }
  else {
    return(res)
  }
}


# Function that shows the expression profile for several regions -------------

SeveralRegions <- function(sample,
                           left_starting_point=NULL, # typically portal
                           right_starting_point =NULL, #typically central
                           dist = 400,
                           GOI) {
  
  total.frame = data.frame()
  total.frame.2 = data.frame()
  
  #Computes the layering regions of the left starting point for each starting point
  if(!is.null(left_starting_point)){
    for (i in 1:length(left_starting_point)) {
      
      start = left_starting_point[i]
      df_temp = SingleRegion(sample, start, dist, GOI, trig = TRUE)
      total.frame = rbind(total.frame, df_temp)
      
    }
  }
  
  #Computes the layering regions of the right starting point for each starting point
  if(!is.null(right_starting_point)){
    for (j in 1:length(right_starting_point)) {
      
      start.2 = right_starting_point[j]
      df_temp.2 = SingleRegion(sample, start.2, dist, GOI, trig = TRUE)
      total.frame.2 = rbind(total.frame.2, df_temp.2)
      
    }
  }
  
  #Prints the plot if only right region of origin
  if(is.null(left_starting_point)){
    for (gene in GOI) {
      df_temp.2 <- total.frame.2[GOI==gene,]
      print(ggplot(df_temp.2,aes(x=distance, y=expr)) 
             + geom_smooth(aes(fill="blue"),level = .7,  colour="red",formula = y ~ x, method = "loess")   
             + theme_minimal()  + scale_fill_discrete(labels = c("Central")) 
             +  labs(fill="Starting Point") + ggtitle(label=as.character(gene)) 
             + theme(plot.title = element_text(hjust = 0.5)) ) 
    }
  } 
  
  
  # Prints the plot if only left region of origin
  else if(is.null(right_starting_point)) {
    for (gene in GOI) {
      df_temp <- total.frame[GOI==gene,]
      print(ggplot(df_temp,aes(x=distance, y=expr)) 
             + geom_smooth(aes(fill="blue"),level = .7,  colour="red",formula = y ~ x, method = "loess")   
             + theme_minimal()  + scale_fill_discrete(labels = c("Portal")) 
             +  labs(fill="Starting Point") + ggtitle(label=as.character(gene)) 
             + theme(plot.title = element_text(hjust = 0.5)) ) 
    }
  } 
  
  #prints the plot if both regions are selected 
  else {
    for (gene in GOI) {
      df_temp <- total.frame[GOI==gene,]
      df_temp.2 <- total.frame.2[GOI==gene,]
      print(ggplot(df_temp,aes(x=distance, y=expression)) 
             + geom_smooth(aes(fill="red"),
                           level = .7,
                           colour="blue",
                           formula = y ~ x, method = "loess")   
             + geom_smooth(data = df_temp.2,
                           aes(x = distance,
                               y = expression, 
                               fill="blue"),
                           colour="red",
                           level = .7,  formula = y ~ x, method = "loess") 
             + theme_minimal() 
             + scale_fill_discrete(labels = c("Central", "Portal")) 
             +labs(fill="Starting Point", 
                   x= "Distance from starting point",
                   y = "Expression") 
             + ggtitle(label=as.character(gene)) 
             + theme(plot.title = element_text(hjust = 0.5, size=16, face="bold")) ) 
    }
  }
}

# Calculates euclidian distance between two points -----------------------------
euc.dist <- function(x1, x2) {return(sqrt(sum((x1 - x2)^2)))}

# Sampling of spots ------------------------------------------------------------
start.sampling <- function(sample,
                           nsample = 10,
                           zone,
                           spacing = 0,
                           seed = 42) {
  
  starter.list = c()
  list.of.spots = colnames(sample[,which(sample$Bio.Order==zone)])
  first.candidate = sample(list.of.spots, 1)
  set.seed(seed)
  starter.list = append(starter.list, first.candidate)
  it = 1 
  while(length(starter.list)<nsample) {
    
    candidate = sample(list.of.spots, 1)
    candidate.coords = as.matrix(sample@images$slice1@coordinates[candidate,2:3])
    flag = FALSE
    
    #Si il y a un point précédent proche on déclenche un flag
    for (i in starter.list) {
      previous.coords = as.matrix(sample@images$slice1@coordinates[i,2:3])
      if (euc.dist(previous.coords, candidate.coords) < spacing) {
        flag = TRUE
      }  
    }
    
    if (!flag) {
      starter.list = append(starter.list, candidate)
      it = it+1
    }
  }
  return(starter.list)
}


# Function that goes from one region to another linearily --------------------

SingleLinearRegion <- function(sample,
                               start = NULL,
                               target_region,
                               GOI,
                               ID =1,
                               show_region = F) {
  
  pt.start = start
  sample@meta.data$imagerow <- sample@images$slice1@coordinates[ , "imagerow" ]
  sample@meta.data$imagecol <- sample@images$slice1@coordinates[ , "imagecol" ]
  
  #fetch the coordinates of the starting point
  pt_coords =  sample@images$slice1@coordinates[pt.start,]
  row = pt_coords$row/(1/sqrt(0.75)) #rescaling so 1 dimension unit == 100 micron
  col = pt_coords$col/2
  imrow = pt_coords$imagerow
  imcol = pt_coords$imagecol
  
  #fetch all of teh coordinates
  total_coords = sample@images$slice1@coordinates
  total_coords$row = total_coords$row/(1/sqrt(0.75)) #rescaling so 1 dimension unit == 100 micron
  total_coords$col = total_coords$col/2
  
  #Keep the coordinates in memory
  memory_coords = total_coords
  
  # Chenge the origin of the coordinate system
  total_coords$row = total_coords$row - total_coords[pt.start,]$row
  total_coords$col = total_coords$col - total_coords[pt.start,]$col
  
  #Switch the base ?
  x = total_coords$row
  y = total_coords$col
  R = sqrt(x^2 + y^2)
  theta = 2*atan(y/(x + sqrt(x^2 + y^2))) # in radians
  theta = theta*180/pi
  polar_df = cbind(total_coords, "Radius" = R, "theta" = theta)
  
  #Find among the spots of the target zone the closest one to the starting point
  sub = polar_df[colnames(sample[,which(sample$Bio.Order == target_region)]),]
  sub <-  sub[order(sub$Radius), ]
  
  closest.arrival = rownames(sub[which(sub$Radius == min(sub$Radius)),])
  
  
  #Create a line from the starting point to the arrival point
  target.theta = polar_df[closest.arrival,]$theta
  target.radius = polar_df[closest.arrival,]$Radius
  
  spot.in.that.line = rownames(polar_df[which( polar_df$theta <= (target.theta+20) &  polar_df$theta >= (target.theta-20) & polar_df$Radius <= (target.radius+1)) ,])
  
  no_spot = length(spot.in.that.line) == 0
  
  #In case the closest spot is in fact too close, take the 2nd closest spot
  if (no_spot){
    closest.arrival = rownames(sub[2, ])
    target.theta = polar_df[closest.arrival,]$theta
    target.radius = polar_df[closest.arrival,]$Radius
    
    spot.in.that.line = rownames(polar_df[which( polar_df$theta <= (target.theta+20) &  polar_df$theta >= (target.theta-20) & polar_df$Radius <= (target.radius+1)) ,])
  
    print(SpatialDimPlot(sample, 
                   cells.highlight = c(closest.arrival, pt.start, spot.in.that.line),
                   pt.size.factor = 1) + labs(fill = "Starting Point"))
    
    }
  
  # Make the line larger to get adjacent spots and create a rectangle
  spot.in.that.rectangle = NULL
  for (i in 1:length(spot.in.that.line)) {
    
    mem.row = memory_coords[spot.in.that.line[i],c("row")]
    mem.col = memory_coords[spot.in.that.line[i],c("col")]
    
    spot.in.that.rectangle = append(spot.in.that.rectangle,
                                    rownames(memory_coords[which(memory_coords$row <= (mem.row+1.1) &  memory_coords$row >= (mem.row-1.1) & memory_coords$col <= (mem.col+1.1) & memory_coords$col >= (mem.col-1.1)  ),])
    )
  }
  
  #keep only spots within the rectangle
  studied.rectangle = polar_df[which(rownames(polar_df)%in%spot.in.that.rectangle),]
  
  #CCalculation of the mean expression level in each of the spots
  max.radius = max(studied.rectangle$Radius)
  radius.list = seq(0.01, (max.radius+1), 1)
  
  res <- data.frame(matrix(ncol = 4, nrow = 0))
  resnames <- c("GOI", "barcode", "expression", "dist")
  
  sample@meta.data$visual = rep(0, ncol(sample))
  
  for (i in 1:length(radius.list)-1){
    a = radius.list[i]
    b = radius.list[i+1]
    
    if (i == 0){
      spot_list <-start
    }
    
    else { 
      spot_list <- rownames(studied.rectangle[which(studied.rectangle$Radius>= a & studied.rectangle$Radius<=b),])
      }
    
    
    expr = as.data.frame(rowMeans(as.data.frame(FetchData(sample, GOI, cells = spot_list, layer = "data", assay = "SCT"))))
    colnames(expr) = "expression"
    
    if (nrow(expr) == 0){
      next
    }
    
    generep <- rep(GOI, nrow(expr))
    distrep <- rep(i, nrow(expr))
    tmp_df <- data.frame(matrix(nrow = nrow(expr)))
    tmp_df$GOI <- generep
    tmp_df$barcode <- rownames(expr)
    tmp_df$expression <- expr$expression
    tmp_df$dist <- distrep
    res <- rbind(res, tmp_df)
    
    sample@meta.data[spot_list,]$visual = i
  }
  
  res <- res[, 2:5]
  res$ID <- rep(ID, length(res$GOI))
  
  my_colors = c("black","#ff4040" ,"#ff7f00", "#ffd700","#00cd00","#89cdff","#3399ff","#0000ff", "#6614b8", "white")
  
  #Show the rectangle region
  if (show_region) {
    print(SpatialDimPlot(sample, group.by = "visual")  + 
            scale_fill_manual(values = my_colors))
  }
  

  
  return(res)
  
  }




#function that repeats the linear regionning -----------------------------------

CombineLinearRegioning <- function(sample,
                                   starts,
                                   target_region,
                                   GOI){
  
  total_frame = SingleLinearRegion(sample,
                                   starts[1],
                                   target_region = target_region,
                                   GOI = GOI,
                                   ID = 1)

  
  if(!is.null(starts)){
    for (i in 2:length(starts)) {
      start = starts[i]
      df_temp = SingleLinearRegion(sample,
                                   start,
                                   target_region = target_region,
                                   GOI = GOI,
                                   ID = i)
      total_frame = rbind(total_frame, df_temp)
    }
  }
  
  MaxLayer = max(total_frame$dist)
  interpolate_df = data.frame()
  return(total_frame)
  for(id in 1:max(total_frame$ID)) {
    sub = total_frame[total_frame$ID == id,]
    #Add NAs to intermedierary missing layers
    while(length(sub$expression)<MaxLayer) {
      sub <- sub %>% add_row(expr = NA, ID = id, GOI = GOI, .before = (floor(length(sub$expression)/2)+1))
    }
    sub$dist = seq(1,MaxLayer)
    # Interpolate the missing lines
    sub$expression = na.spline(sub$expression)
    interpolate_df = rbind(interpolate_df, sub)
  }
  
  
  return(interpolate_df)
  
  
}


# Function that runs linear regionning -----------------------------------------

RunLinearRegionning <- function(sample,
                                nsample,
                                start_region,
                                target_region,
                                spacing = 1,
                                GOI,
                                seed = 42,
                                return_df = F) {
  
  color_grad = rev(c("#ce1212","#f57a00","#e6d11d","#6FD577","#048FEB","#240CE1"))
  
  start_list <- start.sampling(sample, nsample, start_region, spacing = spacing, seed = seed)
  
  plot_df<- CombineLinearRegioning(sample = sample,
                                   starts = start_list,
                                   target_region = target_region,
                                   GOI = GOI)
  plot_df <- plot_df %>%
    group_by(ID)%>%
    mutate(max_dist = max(dist))
  
  plot_df$scaled_dist <- plot_df$dist/plot_df$max_dist

  max_dist_counts <- plot_df %>%
    group_by(ID, max_dist) %>%
    summarise(count = n()) %>%
    ungroup() 
  
  
  p3 <- ggplot(max_dist_counts, aes(x=max_dist)) +
    geom_histogram()+ 
    xlab("Minimal between portal and central")
  
  print(p3)
  
  
  
  p1 <- ggplot(plot_df, aes(x = scaled_dist, y = expression, fill = GOI)) +
    geom_smooth(method = "loess", colour="black") +
    theme_minimal() + ggtitle(label = plot_df$GOI) + 
    xlab("Minimal distance between portal and central") + 
    theme(plot.title = element_text(hjust = 0.5, size=16, face="bold")) 
  

    
  
  # p2 <- VlnPlot(sample,
  #              features = GOI, 
  #              layer = "data",
  #              group.by = "Bio.Order", 
  #              cols = color_grad) +
  #   geom_smooth(data = plot_df,
  #               aes(x = scaled_dist,
  #                  y = expression,
  #                  fill="blue"),
  #               colour="black",
  #               level = .7,  formula = y ~ x, method = "loess")+
  #   scale_x_continuous(sec.axis = sec_axis(~./6))+
  #   theme(legend.position = "none")
  # 
  print(p1)
  #print(p2)
  
  if (return_df){
    return(plot_df)
  }
}


# Function that runs hexagonal regionning --------------------------------------

RunHexagonalRegionning <- function(sample,
                                   nsample,
                                   left_region = NULL,
                                   right_region = NULL,
                                   distance = 400,
                                   spacing = 1,
                                   GOI,
                                   seed = 42) {
  
  if (!is.null(left_region)){
    left_start = start.sampling(sample, nsample, left_region, spacing = spacing, seed = seed)
  }
  
  else {left_start = NULL}
  
  if (!is.null(right_region)){
    right_start = start.sampling(sample, nsample, right_region, spacing = spacing, seed = seed)
  }
  else {right_start = NULL}
  
  
  p1 <- SeveralRegions(sample,
                       left_starting_point  = left_start,
                       right_starting_point  = right_start,
                       dist=distance, GOI=GOI) + 
    ggtitle(label = GOI) 
  print(p1)
  
}


# testing zone. ---------------------------------------------------------------

testing = F

if (testing) {
  
  sample <- readRDS("data/clean/samp_1340_bio_ordered.Rds")
  sample_1352 <- readRDS("data/clean/samp_1352_bio_ordered.Rds")
  sample_1353 <- readRDS("data/clean/samp_1353_bio_ordered.Rds")
  sample_1354 <- readRDS("data/clean/samp_1354_bio_ordered.Rds")
  
  levels(Idents(sample_1352)) <- c("Portal", "PeriPortal", "MidPortal", "MidCentral", "PeriCentral","Central")
  VlnPlot(sample_1352, "Glul", cols = color_grad)
  VlnPlot(sample_1352, "Spp1", cols = color_grad)
  VlnPlot(sample_1352, "Cyp2f2", cols = color_grad)
  VlnPlot(sample_1352, "Mup11", cols = color_grad)
  
  levels(Idents(sample_1353)) <- c("Portal", "PeriPortal", "MidPortal", "MidCentral", "PeriCentral","Central")
  VlnPlot(sample_1353, "Glul", cols = color_grad)
  VlnPlot(sample_1353, "Spp1", cols = color_grad)
  VlnPlot(sample_1353, "Cyp2f2", cols = color_grad)
  VlnPlot(sample_1353, "Mup11", cols = color_grad)
  

  VlnPlot(sample_1354, "Glul", cols = color_grad)
  VlnPlot(sample_1354, "Spp1", cols = color_grad)
  VlnPlot(sample_1354, "Cyp2f2", cols = color_grad)
  VlnPlot(sample_1354, "Mup11", cols = color_grad)
  
  
  
  
  gene = "Mup11"
  samp1 <- RunLinearRegionning(sample_1354, 40, start_region = "1", target_region = "6", GOI = gene, seed =20, return_df = T)
  
  ggplot(samp1, aes(x = scaled_dist, y = expression, fill = GOI)) +
    geom_smooth(method = "loess", colour="black") +
    theme_minimal() + ggtitle(label = samp1$GOI) + 
    xlab("Minimal distance between portal and central") + 
    theme(plot.title = element_text(hjust = 0.5, size=16, face="bold")) 
  
  
  
  samp1$sample <- "CTRL"
  samp2$sample <- "STIM"
  samp3$sample <- "DaKO"
  samp4$sample <- "DaKO_STIM"
  
  big_df <- rbind(hehe, samp2, samp3, samp4)
  
  p1 <- ggplot(big_df, aes(x = scaled_dist, y = expression, fill = sample, group = sample)) +
    geom_smooth(method = "loess", colour="black") +
    theme_minimal() + ggtitle(label = big_df$GOI) + 
    xlab("Minimal distance between portal and central")+
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  
  p1
  
  
  
  
  
  
  RunHexagonalRegionning(sample_1354, 20, left_region = "1", right_region = "6", GOI = "Spp1", distance = 400)
  RunHexagonalRegionning(sample_1354, 20, left_region = "1", right_region = "6", GOI = "Cyp2f2", distance = 400)
  RunHexagonalRegionning(sample_1353, 20, left_region = "1", right_region = "6", GOI = "Mup11", distance = 400)
  RunHexagonalRegionning(sample_1354, 20, left_region = "1", right_region = "6", GOI = "Glul", distance = 400)
  
  SpatialDimPlot(sample)+ scale_fill_manual(values = color_grad)
  DimPlot(sample, cols = color_grad)
  SpatialDimPlot(sample_1354) + scale_fill_manual(values = color_grad)
  DimPlot(sample_135, cols = color_grad)
  
  nsample = 10
  spacing = 1
  C1.portal.start = start.sampling(sample, nsample, "1", spacing=spacing)
  C1.central.start = start.sampling(sample, nsample, "6", spacing=spacing)
  
  gene = "Cyp2f2"
  p1 <- SeveralRegions(sample,
                       left_starting_point  = C1.portal.start,
                       right_starting_point  = C1.central.start,
                       dist=400, GOI=gene) + ggtitle("Contrôle", subtitle = gene) + theme(legend.position = "none")
  p1
  
  
  res = SingleLinearRegion(sample, C1.portal.start[10], target_region = "6", GOI = "Glul")
  
  plot_spp1 = CombineLinearRegioning(sample = sample,
                                     starts = C1.portal.start,
                                     target_region = "6",
                                     GOI = "Spp1")
  plot_spp1 <- plot_spp1 %>%
    group_by(ID)%>%
    mutate(max_dist = max(dist))
  
  
  
  plot_apoc <- CombineLinearRegioning(sample = sample,
                                      starts = C1.portal.start,
                                      target_region = "6",
                                      GOI = "Apoc2")
  plot_cyp2 <- CombineLinearRegioning(sample = sample,
                                      starts = C1.portal.start,
                                      target_region = "6",
                                      GOI = "Cyp2f2")
  plot_mup <- CombineLinearRegioning(sample = sample,
                                     starts = C1.portal.start,
                                     target_region = "6",
                                     GOI = "Mup11")
  plot_scd <- CombineLinearRegioning(sample = sample,
                                     starts = C1.portal.start,
                                     target_region = "6",
                                     GOI = "Rnase4")
  plot_glul <- CombineLinearRegioning(sample = sample,
                                      starts = C1.central.start,
                                      target_region = "1",
                                      GOI = "Glul")
  
  
  odd = c("GCACAACCTCGGGCGT-1" ,"GACTCCTTCCAATTGA-1")
  spat <- SpatialDimPlot(sample, 
                         cells.highlight = odd,
                         pt.size.factor = 1) + labs(fill = "Starting Point")
  spat
  
  
  Spatial
  
  
  cols = c("Spp1" = "#240CE1",
           "Apoc2" = "#048FEB",
           "Cyp2f2" = "#6FD577",
           "Mup11" = "#e6d11d",
           "Rnase4" = "#f57a00",
           "Glul" = "#ce1212")
  
  ggplot(plot_spp1, aes(x = dist*100, y = expression, fill = GOI)) +
    geom_smooth(method = "loess", colour="black") +
    theme_minimal() + ggtitle(label = plot_spp1$GOI) +
    geom_smooth(data = plot_apoc,
                aes(x = dist*100,
                    y = expression,
                    fill=GOI),
                colour="black",
                level = .7,  formula = y ~ x, method = "loess")+
    geom_smooth(data = plot_cyp2,
                aes(x = dist*100,
                    y = expression,
                    fill=GOI),
                colour="black",
                level = .7,  formula = y ~ x, method = "loess")+
    geom_smooth(data = plot_mup,
                aes(x = dist*100,
                    y = expression,
                    fill=GOI),
                colour="black",
                level = .7,  formula = y ~ x, method = "loess")+
    geom_smooth(data = plot_scd,
                aes(x = dist*100,
                    y = expression,
                    fill=GOI),
                colour="black",
                level = .7,  formula = y ~ x, method = "loess")+
    geom_smooth(data = plot_glul,
                aes(x = dist*100,
                    y = expression,
                    fill=GOI),
                colour="black",
                level = .7,  formula = y ~ x, method = "loess")+
    scale_fill_manual(values = cols) +
    xlab("Trajet minimal entre portal vers central") 
  
  
  plot_df<- CombineLinearRegioning(sample = sample_1352,
                                   starts = C1.portal.start,
                                   target_region = "6",
                                   GOI = "Fgf21")
  
  color_grad = rev(c("red", "#ce1212","#f57a00","#e6d11d","#6FD577","#048FEB","#240CE1"))
  
  VlnPlot(sample_1352,
          features = "Fgf21", 
          layer = "data",
          group.by = "Bio.Order", 
          cols = color_grad) +
    geom_smooth(data = plot_df,
                aes(x = dist+0.5,
                    y = expression,
                    fill="blue"),
                colour="black",
                level = .7,  formula = y ~ x, method = "loess")+ theme(legend.position = "none")
  
}

p1 <- ggplot(hehe, aes(x = - scaled_dist, y = expression, fill = GOI)) +
  geom_smooth(method = "loess", colour="black") +
  theme_minimal() + ggtitle(label = hehe$GOI) + 
  xlab("Minimal distance between portal and central") 
p1

  

