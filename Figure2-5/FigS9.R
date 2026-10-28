# =====================================================================================================
# Fig. S
# =====================================================================================================
# INPUTS: 
#    chromosome2.tsv, distances_chr2L100-3000.txt
# =====================================================================================================

library(plotly)
library(ggridges)
library(circlize)
library(ggExtra)
library(tidyr)
library(dbscan)

IMR_2 <- read.delim("chromosome2.tsv")
print(colnames(IMR_2))
print(dim(IMR_2))

cells = unique(IMR_2$Chromosome.copy.number)
spots = unique(IMR_2$Genomic.coordinate)
print(paste(length(cells),length(spots)))

Co = IMR_2$Genomic.coordinate
Coordinate = rep(0,nrow = length(Co))
for (ii in 1:length(spots)) {
  Co_tmp = which(Co==spots[ii])
  tmp = strsplit(spots[ii],split =':')[[1]]
  tmp = strsplit(tmp[2],split ='-')[[1]]
  Coordinate[Co_tmp] =(as.numeric(tmp[1]))
}
IMR_2$Coordinate = Coordinate

spots2 = unique(Coordinate)

st = 0
ed = 0
for (jj in 1:length((spots2))) {
  st = c(st,rep(spots2[jj],length(spots2)-jj))
  ed = c(ed, spots2[-c(1:jj)])
}
st = st[-1]
ed = ed[-1]
lner = as.numeric(ed)-as.numeric(st)
tt1 = which(lner>=100000000)

distances_2m100 = read.table(paste0('distances_chr2L100-3000.txt'),skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

chr2 = 242200000
rsl = 50
d3li = 299

conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

rmv = which(colSums(conformation_matrix)>10000)
conformation_matrix = conformation_matrix[,-rmv]
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))
colnames(distances_stat2) = c('start','end','genom.dist.','numblks')
distances_stat2 = data.frame(distances_stat2)
print(dim(distances_stat2))
print(dim(conformation_matrix))
print(head(distances_stat2))
print(head(IMR_2))

bin_thresh = 1000000
tmp = which(distances_stat2[,4]>23)
distances_stat22 = distances_stat2[tmp,]
conformation_matrix22 = conformation_matrix[tmp,]
print(dim(distances_stat22))
print(dim(conformation_matrix22))

####################################################################################################
# Fig. S9d, S9e ####
####################################################################################################
# position of loci on chr2 250k ####

# Perform Topological Data Analysis
analyze_point_distribution <- function(points_df) {
  points_df2 = data.frame(x = points_df$x, y = points_df$y, z = points_df$z)
  # Convert dataframe to matrix
  points_matrix <- as.matrix(points_df2)/1000
  
  # Compute distance matrix
  dist_matrix <- as.matrix(dist(points_matrix))
  
  # Persistent Homology Analysis
  # We'll look at 0-dimensional (connected components) and 1-dimensional (holes) features
  max_distance <- max(dist_matrix)
  
  # # Compute persistent homology
  # ph_diagram <- gridDiag(
  #   X = points_matrix,
  #   FUN = kde,  # Kernel Density Estimation
  #   lim = t(cbind(
  #     apply(points_matrix, 2, min),
  #     apply(points_matrix, 2, max)
  #   )),
  #   by = 0.1,
  #   maxdimension = 1,
  #   sublevel = FALSE,
  #   epsilon = 0.1,
  #   sigma = 0.1,
  #   n = 40
  # )
  # 
  # # Persistent Homology Diagram
  # plot(ph_diagram$diagram, 
  #      main = "Persistent Homology Diagram",
  #      xlab = "Birth", 
  #      ylab = "Death")
  # 
  # Identify peripheral points using DBSCAN
  # This helps identify points that are less densely surrounded
  dbscan_result <- dbscan(points_matrix, eps = 0.3, minPts = 3)
  
  # Compute local point density
  knn_dist <- kNNdist(points_matrix, k = 5)
  
  # Create results dataframe
  results <- data.frame(
    points_df2/1000,
    cluster = dbscan_result$cluster,
    local_density = 1 / knn_dist,
    coord = points_df$cord
  )
  
  # Return comprehensive results
  return(list(
    # persistent_homology = ph_diagram,
    point_analysis = results,
    dbscan_clusters = dbscan_result
  ))
}
# Function to filter points with different 3d distancne from targeted points
filter_nearby_points <- function(point_matrix,lc) {
  point_matrix2 = as.matrix(point_matrix[,2:4])
  # Extract the third point's coordinates
  reference_point <- point_matrix2[lc,]
  
  # Calculate Euclidean distances for all points from the reference point
  distances <- apply(point_matrix2, 1, function(point) {
    sqrt(sum((point - reference_point)^2))
  })
  
  # Filter points within xx nm
  nearby_points <- point_matrix[distances <= 299,]
  nearby_points = abs(nearby_points$cord-point_matrix$cord[lc])
  neigbhor = which((nearby_points>0)&(nearby_points<=1))
  near = which((nearby_points>1)&(nearby_points<=10))
  medi = which((nearby_points>10)&(nearby_points<=50))
  long = which((nearby_points>50)&(nearby_points<=100))
  superl = which((nearby_points>100))
  return(c(length(neigbhor),length(near),length(medi),length(long),length(superl)))
}
# Function to calculate angular distribution
calculate_angular_distribution <- function(central_point, surrounding_points, radius) {
  # Convert to spherical coordinates relative to central point
  calculate_spherical_coords <- function(point) {
    # Relative vector
    rel_vector <- point - central_point
    
    # Distance from central point
    r <- sqrt(sum(rel_vector^2))
    
    # Check if point is within the shell radius
    if (r > radius) return(NULL)
    
    # Azimuthal angle (theta) in x-y plane
    theta <- atan2(rel_vector[2], rel_vector[1])
    
    # Polar angle (phi) from z-axis
    phi <- acos(rel_vector[3] / r)
    
    # Additional check to determine which side of the diameter plane
    # Using z-coordinate to split each sector
    z_sign <- sign(rel_vector[3])
    
    return(c(r = r, theta = theta, phi = phi, z_sign = z_sign))
  }
  
  # Apply to all surrounding points
  spherical_coords <- lapply(
    1:nrow(surrounding_points), 
    function(i) calculate_spherical_coords(as.numeric(surrounding_points[i,]))
  )
  
  # Remove NULL values (points outside radius)
  spherical_coords <- spherical_coords[!sapply(spherical_coords, is.null)]
  
  # Convert to matrix
  spherical_matrix <- do.call(rbind, spherical_coords)
  spherical_matrix = spherical_matrix[which(abs(spherical_matrix[,4])==1),]
  
  # Divide sphere into sectors with additional subdivision
  detailed_sector_analysis <- function(coords, num_primary_sectors = 8) {
    # Divide theta (azimuthal) into primary sectors
    theta_breaks <- seq(-pi, pi, length.out = num_primary_sectors + 1)
    theta_sectors <- cut(coords[,'theta'], breaks = theta_breaks)
    
    # Combine primary sector with z-sign for detailed sectoring
    detailed_sectors <- interaction(theta_sectors, coords[,'z_sign'])
    
    # Count points in each detailed sector
    sector_distribution <- table(detailed_sectors)
    
    # Calculate variance of sector point counts
    sector_variance <- var(sector_distribution)
    
    return(list(
      sector_distribution = sector_distribution,
      sector_variance = sector_variance,
      detailed_sectors = detailed_sectors
    ))
  }
  
  # Perform detailed sector analysis
  distribution_analysis <- detailed_sector_analysis(spherical_matrix)
  
  return(distribution_analysis)
}

# Main execution
locus_density = matrix(0, nrow = 1, ncol = 10)
locus_distribution = matrix(0, nrow = 1, ncol = 16)
for (ss in (1:3000)) {
  if (ss%%10==0) {
    print(ss)
    print(dim(locus_density))
    print(dim(locus_distribution))
  }
  tmp = which(IMR_2$Chromosome.copy.number==ss)
  data <- data.frame(
    index = 1:length(tmp),
    x = IMR_2$X.nm.[tmp],
    y = IMR_2$Y.nm.[tmp],
    z = IMR_2$Z.nm.[tmp],
    cord = (IMR_2$Coordinate[tmp]-1)/1000000
  )
  data = na.omit(data)

  angular_results = rep(0,nrow(data))
  locus_distribution_tmp = matrix(0,nrow(data),16)
  for (pp in 1:nrow(data)) {
    surrounding_points = data[,2:4]
    # Calculate Euclidean distances for all points from the reference point
    distances <- apply(surrounding_points, 1, function(point) {
      sqrt(sum((point - as.numeric(data[pp,2:4]))^2))
    })
    distances = as.numeric(distances)
    tmp = which(distances<299)
    if (length(tmp)<7) {
      next
    }else{
      surrounding_points = surrounding_points[tmp,]
      tmp <- calculate_angular_distribution(as.numeric(data[pp,2:4]), surrounding_points, 299)
      angular_results[pp] = round(tmp$sector_variance,3)
      locus_distribution_tmp[pp,] = as.numeric(tmp$sector_distribution)
    }
  }
  locus_distribution = rbind(locus_distribution,locus_distribution_tmp)
  dis_tmp = matrix(0,nrow = nrow(data), ncol = 6)
  for (ll in 1:nrow(data)) {
    dis_tmp[ll,] = c(data$cord[ll],filter_nearby_points(data,ll))
  }
  analysis_results <- analyze_point_distribution(data)
  dis_tmp = cbind(rep(ss,nrow(data)),dis_tmp,analysis_results$point_analysis$local_density,analysis_results$point_analysis$coord,angular_results)
  tmp = diff(c(length(rep(ss,nrow(data))),nrow(dis_tmp),length(analysis_results$point_analysis$local_density),length(analysis_results$point_analysis$coord),length(angular_results)))
  if (length(which(diff(tmp)!=0))>0) {
    print(ss)
    break
  }
  locus_density = rbind(locus_density,dis_tmp)
}
locus_distribution = locus_distribution[-1,]
locus_density = locus_density[-1,]
colnames(locus_density) = c('chr_copy','locus','neigbhor','near','median','long','superlong','dbscan','locus2','sector_var')
locus_density = data.frame(locus_density)
locus_density = locus_density[,-9]

locus_pos = matrix(0,nrow = length(unique(locus_density$locus)), ncol = 6)
locus_pos[,1] = sort(unique(locus_density$locus))
for (ii in 1:nrow(locus_pos)) {
  tmp = locus_pos[ii,1]
  neigbhor_tmp = seq(max(0,tmp-1),min(242,tmp+1),0.25)
  near_tmp = seq(max(0,tmp-10),min(242,tmp+10),0.25)
  median_tmp = seq(max(0,tmp-50),min(242,tmp+50),0.25)
  long_tmp = seq(max(0,tmp-100),min(242,tmp+100),0.25)
  superlong_tmp = seq(0,242,0.25)
  superlong_tmp = superlong_tmp[!(superlong_tmp%in%long_tmp)]
  long_tmp = long_tmp[!(long_tmp%in%median_tmp)]
  median_tmp = median_tmp[!(median_tmp%in%near_tmp)]
  near_tmp = near_tmp[!(near_tmp%in%neigbhor_tmp)]
  locus_pos[ii,2:6] = c(length(neigbhor_tmp),length(near_tmp),length(median_tmp),length(long_tmp),length(superlong_tmp))
}
colnames(locus_pos) = c('locus','neibhor','near','median','long','superlong')


####################################################################################################
# save the result
# SOURCE: Figure5.R, lines 1051-1051   | PANELS: S9d, S9e
####################################################################################################
save.image("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/locus_density2.RData")

####################################################################################################
# load / remove the copies with > 10000 LRC
# SOURCE: Figure5.R, lines 1052-1058   | PANELS: S9d, S9e
####################################################################################################
load("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/locus_density2.RData")
print(dim(locus_density))
print(head(locus_density))
locus_pos = data.frame(locus_pos)
print(dim(locus_pos))
print(head(locus_pos))
locus_density = locus_density[!locus_density$chr_copy%in%rmv,]

####################################################################################################
# local density vs sLRC: RD-SN / RDH at both ends of every LRC (3D distance < 299 nm) of every chr2 copy (XRC_density) and an equally sized random non-LRC set (XRC_density_notLRC)
# SOURCE: Figure5.R, lines 1580-1674   | PANELS: S9d, S9e
# NOTE: 'XRC' = LRC. sLRC are picked out afterwards by matching the locus pair to distances_stat22.
####################################################################################################
## local density v.s. sLRC ####
load("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/locus_density2.RData")
print(dim(locus_density))
print(head(locus_density))
print(head(locus_pos))
sLRC = which(distances_stat2[,4]>23)
library(data.table)
# Create a list to store matrices
XRC_density <- vector("list", 3000)
for (ss in 1:ncol(distances_matrix)) {

  # print(c(ss,log10(rcr[ss,1]),log10(rcr[ss,2])))
  distances_stat_tmp = cbind(distances_stat2,as.numeric(as.character(distances_matrix[,ss])))
  colnames(distances_stat_tmp) = c(colnames(distances_stat2),'D3dist')
  distances_stat_tmp$start = (distances_stat_tmp$start-1)/1000000
  distances_stat_tmp$end = (distances_stat_tmp$end-1)/1000000
  
  locus_density_tmp = locus_density[which(locus_density$chr_copy==ss),]
  
  st_tmp = match(distances_stat_tmp$start,locus_density_tmp$locus)
  ed_tmp = match(distances_stat_tmp$end,locus_density_tmp$locus)
  
  LRC_dens = data.frame(
    locus_st = locus_density_tmp$locus[st_tmp],
    locus_ed = locus_density_tmp$locus[ed_tmp],
    dbscan_st = locus_density_tmp$dbscan[st_tmp],
    dbscan_ed = locus_density_tmp$dbscan[ed_tmp],
    near_st = locus_density_tmp$neigbhor[st_tmp]+locus_density_tmp$near[st_tmp],
    near_ed = locus_density_tmp$neigbhor[ed_tmp]+locus_density_tmp$near[ed_tmp],
    median_st = locus_density_tmp$median[st_tmp],
    median_ed = locus_density_tmp$median[ed_tmp],
    long_st = locus_density_tmp$long[st_tmp],
    long_ed = locus_density_tmp$long[ed_tmp],
    superlong_st = locus_density_tmp$superlong[st_tmp],
    superlong_ed = locus_density_tmp$superlong[ed_tmp],
    distribution_st = locus_density_tmp$sector_var[st_tmp],
    distribution_ed = locus_density_tmp$sector_var[ed_tmp]
    )

  LRC = which(distances_stat_tmp$D3dist<299)
  LRC_dens = LRC_dens[LRC,]
  if(ss%%500==0){
    print(c('ss',ss))
    print(dim(distances_stat_tmp))
    print(dim(locus_density_tmp))
    print(dim(LRC_dens))
    }
  current_dt <- as.data.table(LRC_dens)
  current_dt[, chr_copy := ss]
  XRC_density[[ss]] <- current_dt
}

set.seed(123)
XRC_density_notLRC <- vector("list", 3000)
for (ss in 1:ncol(distances_matrix)) {
  # print(c(ss,log10(rcr[ss,1]),log10(rcr[ss,2])))
  distances_stat_tmp = cbind(distances_stat2,as.numeric(as.character(distances_matrix[,ss])))
  colnames(distances_stat_tmp) = c(colnames(distances_stat2),'D3dist')
  distances_stat_tmp$start = (distances_stat_tmp$start-1)/1000000
  distances_stat_tmp$end = (distances_stat_tmp$end-1)/1000000
  
  locus_density_tmp = locus_density[which(locus_density$chr_copy==ss),]
  
  st_tmp = match(distances_stat_tmp$start,locus_density_tmp$locus)
  ed_tmp = match(distances_stat_tmp$end,locus_density_tmp$locus)
  
  LRC_dens = data.frame(
    locus_st = locus_density_tmp$locus[st_tmp],
    locus_ed = locus_density_tmp$locus[ed_tmp],
    dbscan_st = locus_density_tmp$dbscan[st_tmp],
    dbscan_ed = locus_density_tmp$dbscan[ed_tmp],
    near_st = locus_density_tmp$neigbhor[st_tmp]+locus_density_tmp$near[st_tmp],
    near_ed = locus_density_tmp$neigbhor[ed_tmp]+locus_density_tmp$near[ed_tmp],
    median_st = locus_density_tmp$median[st_tmp],
    median_ed = locus_density_tmp$median[ed_tmp],
    long_st = locus_density_tmp$long[st_tmp],
    long_ed = locus_density_tmp$long[ed_tmp],
    superlong_st = locus_density_tmp$superlong[st_tmp],
    superlong_ed = locus_density_tmp$superlong[ed_tmp],
    distribution_st = locus_density_tmp$sector_var[st_tmp],
    distribution_ed = locus_density_tmp$sector_var[ed_tmp]
  )
  
  LRC = which(distances_stat_tmp$D3dist<299)
  LRC_dens = LRC_dens[sample((1:nrow(LRC_dens))[-LRC],length(LRC)),]
  if(ss%%500==0){
    print(c('ss',ss))
    print(dim(distances_stat_tmp))
    print(dim(locus_density_tmp))
    print(dim(LRC_dens))
  }
  current_dt <- as.data.table(LRC_dens)
  current_dt[, chr_copy := ss]
  XRC_density_notLRC[[ss]] <- current_dt
}

####################################################################################################
# save
# SOURCE: Figure5.R, lines 1675-1676   | PANELS: S9d, S9e
####################################################################################################
save.image("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/locus_density2.RData")


####################################################################################################
# load, remove removed copies
# SOURCE: Figure5.R, lines 1677-1699   | PANELS: S9d, S9e
####################################################################################################
library(data.table)
load("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/locus_density2.RData")
# ADDED:
# ADDED: XRC_density / XRC_density_notLRC are lists of per-copy data.tables at this point; data.frame(<list>) below needs ONE table
# (the original does not show this step - the RData it loads must already have contained the combined tables).
if (is.list(XRC_density) && !is.data.frame(XRC_density)) XRC_density = rbindlist(XRC_density)
if (is.list(XRC_density_notLRC) && !is.data.frame(XRC_density_notLRC)) XRC_density_notLRC = rbindlist(XRC_density_notLRC)
locus_pos = data.frame(locus_pos)
print(dim(locus_pos))
print(head(locus_pos))
print(dim(locus_density))
print(head(locus_density))
print(head(XRC_density))
XRC_density = data.frame(XRC_density)
XRC_density <- XRC_density[!XRC_density$chr_copy %in% rmv,]
print(dim(XRC_density))
locus_density = locus_density[!locus_density$chr_copy%in%rmv,]
print(dim(XRC_density))
plot(XRC_density$dbscan_ed,XRC_density$dbscan_st,pch = '.')
hist(XRC_density$dbscan_st)

XRC_density_notLRC = data.frame(XRC_density_notLRC)
XRC_density_notLRC <- XRC_density_notLRC[!XRC_density_notLRC$chr_copy %in% rmv,]
print(dim(XRC_density))
locus_density = locus_density[!locus_density$chr_copy%in%rmv,]
print(dim(XRC_density_notLRC))
plot(XRC_density$dbscan_ed,XRC_density$dbscan_st,pch = '.')
hist(XRC_density$dbscan_st)

####################################################################################################
# normalise by the number of genomic positions per window; RD-SN classes of the two ends (> 6 = high); keep only the sLRC (sRC_density, sRC_density2)
# SOURCE: Figure5.R, lines 1819-1850   | PANELS: S9d, S9e
# NOTE: near_gp: hh = both ends RD-SN > 6, ll = both <= 6, hl/lh = one of each.  sRC_density = RAW counts (used for the plots), sRC_density2 = counts divided by the number of windows (locus_pos).
####################################################################################################
tmp_st = match(XRC_density$locus_st,locus_pos$locus)
tmp_ed = match(XRC_density$locus_ed,locus_pos$locus)

XRC_density2 = XRC_density
XRC_density2$near_st = XRC_density$near_st/(locus_pos$neibhor[tmp_st]+locus_pos$near[tmp_st])
XRC_density2$near_ed = XRC_density$near_ed/(locus_pos$neibhor[tmp_ed]+locus_pos$near[tmp_ed])
XRC_density2$median_st = XRC_density$median_st/(locus_pos$median[tmp_st])
XRC_density2$median_ed = XRC_density$median_ed/(locus_pos$median[tmp_ed])
XRC_density2$long_st = XRC_density$long_st/(locus_pos$long[tmp_st])
XRC_density2$long_ed = XRC_density$long_ed/(locus_pos$long[tmp_ed])
XRC_density2$superlong_st = XRC_density$superlong_st/(locus_pos$superlong[tmp_st])
XRC_density2$superlong_ed = XRC_density$superlong_ed/(locus_pos$superlong[tmp_ed])
XRC_density2$distribution_st = XRC_density$distribution_st
XRC_density2$distribution_ed = XRC_density$distribution_ed

near_gp = rep(0, nrow(XRC_density))
near_hh = which(XRC_density$near_st>6&XRC_density$near_ed>6)
near_gp[near_hh] ='hh'
near_ll = which(XRC_density$near_st<=6&XRC_density$near_ed<=6)
near_gp[near_ll] ='ll'
near_hl = which(XRC_density$near_st>6&XRC_density$near_ed<=6)
near_gp[near_hl] ='hl'
near_lh = which(XRC_density$near_st<=6&XRC_density$near_ed>6)
near_gp[near_lh] ='lh'
XRC_density2$near_gp = near_gp
XRC_density$near_gp = near_gp

LRC = paste(XRC_density$locus_st, XRC_density$locus_ed)
sLRC = paste((distances_stat22$start-1)/1000000,(distances_stat22$end-1)/1000000)
tmp = match(LRC, sLRC)
sRC_density = XRC_density[which(!is.na(tmp)),]
sRC_density2 = XRC_density2[which(!is.na(tmp)),]

####################################################################################################
# Fig. S9e  ####
####################################################################################################
data = XRC_density[sample(1:nrow(XRC_density),100000),]
df <- cbind(data$near_st,data$near_ed,data$near_gp)
colnames(df) <- c("start", "end",'near_gp')
print(dim(df))

data = sRC_density[sample(1:nrow(sRC_density),100000),]
df <- cbind(data$near_st,data$near_ed,data$near_gp)
print(dim(df))

df = as.matrix(df)
df2<-t(apply(t(df), 2, sort, decreasing = F))
df <- as.data.frame(cbind(df,data$near_gp))
colnames(df) <- c("start", "end",'near_gp')
# df = df[which(df$pca2_gp!='grey'),]

ggplot(df, aes(x = as.numeric(start), y = as.numeric(end),colour = near_gp)) +
  geom_point(size = 0.5, alpha = 1) +
  scale_color_manual(values = c("hh" = "#ff9999",
                                "ll" = "#b3c38a",
                                "lh" = "#f7dd83",
                                "hl" = "#f7dd83")) +
  geom_hline(yintercept = 6, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 6, color = "blue", linetype = "dashed") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")  # Add axis lines
  ) +
  labs(title = 'near')

####################################################################################################
# Fig. S9d  ####
####################################################################################################
data <- cbind(XRC_density2$distribution_st,XRC_density2$distribution_ed,XRC_density2$near_gp)
data <- cbind(sRC_density2$distribution_st,sRC_density2$distribution_ed,sRC_density2$near_gp)
colnames(data) <- c("start", "end",'near_gp')
df <- as.data.frame(data)
df = df[sample(1:length(sRC_density2$dbscan_st),100000),]
df$near_gp[which(df$near_gp=='lh')] = 'hl'

df$near_gp <- factor(df$near_gp, levels = c("hh", "ll", "hl"))

ggplot(df, aes(x = as.numeric(start), y = as.numeric(end))) +
  geom_point(size = 0.1, alpha = 0.5, color = case_when(
    df$near_gp == "hh" ~ "#f5bdc8",
    df$near_gp == "ll" ~ "#b3c38a",
    df$near_gp == "hl" ~ "#f7dd83"
  )) +
  facet_wrap(~near_gp, ncol = 3) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    strip.placement = "outside",
    axis.text.x = element_text(angle = 0, hjust = 0.5)  # Ensures x-axis text is shown
  ) +
  labs(title = 'distribution')

ggplot(df, aes(x = as.numeric(start), y = as.numeric(end))) +
  geom_point(aes(color = near_gp), size = 0.1, alpha = 0.5) +
  scale_color_manual(values = c(
    "hh" = "#f5bdc8",
    "ll" = "#b3c38a",
    "hl" = "#f7dd83"
  )) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"  # You can adjust legend position or remove it with "none"
  ) +
  xlim(0,210)+
  ylim(0,210)+
  labs(
    title = 'Combined Distribution',
    color = 'Group'  # This changes the legend title
  )

####################################################################################################
# Fig. S9f ####
####################################################################################################
## 3D plot ####
d3plot <- function(ss, sLRC_pair){
  tmp = which(IMR_2$Chromosome.copy.number==ss)
  set.seed(123)
  data <- data.frame(
    index = 1:length(tmp),
    x = IMR_2$X.nm.[tmp],
    y = IMR_2$Y.nm.[tmp],
    z = IMR_2$Z.nm.[tmp],
    cord = (IMR_2$Coordinate[tmp]-1)/1000000
  )
  print(dim(data))
  
  point_matrix2 = as.matrix(data[,2:4])
  st = which(data$cord==sLRC_pair[1])
  reference_point <- point_matrix2[st,]
  distances <- apply(point_matrix2, 1, function(point) {
    sqrt(sum((point - reference_point)^2))
  })
  nearby_points_st <- data[distances <= 299,]
  nearby_points_st = na.omit(nearby_points_st)
  nearby_points_cord = abs(nearby_points_st$cord-data$cord[st])
  near_st = which((nearby_points_cord>0)&(nearby_points_cord<=10))
  medi_st = which((nearby_points_cord>10)&(nearby_points_cord<=100))
  superl_st = which((nearby_points_cord>100))
  
  ed = which(data$cord==sLRC_pair[2])
  reference_point <- point_matrix2[ed,]
  distances <- apply(point_matrix2, 1, function(point) {
    sqrt(sum((point - reference_point)^2))
  })
  nearby_points_ed <- data[distances <= 299,]
  nearby_points_ed = na.omit(nearby_points_ed)
  nearby_points_cord = abs(nearby_points_ed$cord-data$cord[ed])
  near_ed = which((nearby_points_cord>0)&(nearby_points_cord<=10))
  medi_ed = which((nearby_points_cord>10)&(nearby_points_cord<=100))
  superl_ed = which((nearby_points_cord>100))
  print(c(length(near_st),length(near_ed)))
  print(c(length(superl_st),length(superl_ed)))
  
  # Create a group column in the data
  all_colors = colorspace::rainbow_hcl(length(tmp)+length(tmp))
  data$group <- "Other"
  plot <- plot_ly()
  # Add all points
  plot <- add_trace(plot,
                    data = subset(data), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers+lines",
                    marker = list(size = 3, color = all_colors, symbol = "circle",opacity = 0.3),
                    line = list(color = "grey", width = .5),
                    name = "Other",
                    visible = TRUE,
                    text = data$cord,  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  XRC_indices_st = which(data$cord %in% (nearby_points_st$cord))
  XRC_indices_ed = which(data$cord %in% (nearby_points_ed$cord))
  XRC_indices = c(XRC_indices_st, XRC_indices_ed)
  data$group[XRC_indices_st] <- "XRC_st"
  all_colors[XRC_indices_st] = '#734889'
  data$group[XRC_indices_ed] <- "XRC_ed"
  all_colors[XRC_indices_ed] = '#00B9EF'
  plot <- add_trace(plot,
                    data = subset(data, group == "XRC_st"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 5, color = all_colors[XRC_indices_st], symbol = "circle",opacity = 0.5),
                    # line = list(color = "grey", width = 1),
                    name = "XRC_st",
                    legendgroup = "XRC_st",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="XRC_st")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  plot <- add_trace(plot,
                    data = subset(data, group == "XRC_ed"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 5, color = all_colors[XRC_indices_ed], symbol = "circle",opacity = 0.5),
                    # line = list(color = "grey", width = 1),
                    name = "XRC_ed",
                    legendgroup = "XRC_ed",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="XRC_ed")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  XRC_near_st = which(data$cord %in% (nearby_points_st$cord[near_st]))
  data$group[XRC_near_st] <- "XRC_near_st"
  all_colors[XRC_near_st] = '#9b3f5c'
  plot <- add_trace(plot,
                    data = subset(data, group == "XRC_near_st"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 7, color = '#9b3f5c', symbol = "circle"),
                    # line = list(color = "grey", width = 1),
                    name = "XRC_near_st",
                    legendgroup = "XRC_near_st",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="XRC_near_st")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  XRC_superl_st = which(data$cord %in% (nearby_points_st$cord[superl_st]))
  data$group[XRC_superl_st] <- "XRC_superl_st"
  all_colors[XRC_superl_st] = 'pink'
  plot <- add_trace(plot,
                    data = subset(data, group == "XRC_superl_st"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 7, color = 'pink', symbol = "circle"),
                    # line = list(color = "grey", width = 1),
                    name = "XRC_superl_st",
                    legendgroup = "XRC_superl_st",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="XRC_superl_st")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  XRC_near_ed = which(data$cord %in% (nearby_points_ed$cord[near_ed]))
  data$group[XRC_near_ed] <- "XRC_near_ed"
  all_colors[XRC_near_ed] = '#4d7e54'
  plot <- add_trace(plot,
                    data = subset(data, group == "XRC_near_ed"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 7, color = '#4d7e54', symbol = "circle"),
                    # line = list(color = "grey", width = 1),
                    name = "XRC_near_ed",
                    legendgroup = "XRC_near_ed",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="XRC_near_ed")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  XRC_superl_ed = which(data$cord %in% (nearby_points_ed$cord[superl_ed]))
  data$group[XRC_superl_ed] <- "XRC_superl_ed"
  all_colors[XRC_superl_ed] = '#00B9EF'
  plot <- add_trace(plot,
                    data = subset(data, group == "XRC_superl_ed"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 7, color = '#00B9EF', symbol = "circle"),
                    # line = list(color = "grey", width = 1),
                    name = "XRC_superl_ed",
                    legendgroup = "XRC_superl_ed",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="XRC_superl_ed")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  
  ol1 = intersect(nearby_points_st[near_st,]$cord,nearby_points_ed[superl_ed,]$cord)
  ol2 = intersect(nearby_points_st[superl_st,]$cord,nearby_points_ed[near_ed,]$cord)
  # overlap = intersect(nearby_points_st$cord,nearby_points_ed$cord)
  overlap = c(ol1, ol2)
  if (length(overlap)!=0) {
    XRC_ol = which(data$cord %in% (overlap))
    data$group[XRC_ol] <- "XRC_ol"
    all_colors[XRC_ol] = '#60281e'
    plot <- add_trace(plot,
                      data = subset(data, group == "XRC_ol"), 
                      x = ~x, y = ~y, z = ~z, 
                      type = "scatter3d", mode = "markers",
                      marker = list(size = 7, color = '#60281e', symbol = "circle"),
                      # line = list(color = "grey", width = 1),
                      name = "XRC_ol",
                      legendgroup = "XRC_ol",
                      showlegend = T,
                      visible = TRUE,
                      text = data$cord[which(data$group=="XRC_ol")],  # Use the index as the text
                      hoverinfo = "text",    # Show only the text (index) on hover
                      hoverlabel = list(namelength = 0)) 
  }
  
  tg_indices_st = which(data$cord  ==(sLRC_pair[1]))
  tg_indices_ed = which(data$cord  ==(sLRC_pair[2]))
  tg_indices = c(tg_indices_st, tg_indices_ed)
  data$group[tg_indices] <- "sLRC"
  all_colors[tg_indices] = c('red','blue')
  plot <- add_trace(plot,
                    data = subset(data, group == "sLRC"), 
                    x = ~x, y = ~y, z = ~z, 
                    type = "scatter3d", mode = "markers",
                    marker = list(size =10, color = all_colors[tg_indices], symbol = "circle"),
                    # line = list(color = "grey", width = 1),
                    name = "sLRC",
                    legendgroup = "sLRC",
                    showlegend = T,
                    visible = TRUE,
                    text = data$cord[which(data$group=="sLRC")],  # Use the index as the text
                    hoverinfo = "text",    # Show only the text (index) on hover
                    hoverlabel = list(namelength = 0))
  # Add DistanceGroup links
  pp1 <- which(data$cord  == (sLRC_pair[1]))
  pp2 <- which(data$cord  == (sLRC_pair[2]))
  plot <- add_trace(plot,
                    x = c(data$x[pp1], data$x[pp2]),
                    y = c(data$y[pp1], data$y[pp2]),
                    z = c(data$z[pp1], data$z[pp2]),
                    type = "scatter3d",
                    mode = "lines",
                    line = list(color = "black", width = 10),
                    name = "sLRC",
                    legendgroup = "sLRC",
                    showlegend = F,
                    visible = TRUE)
  
  # Add DistanceGroup links
  pp1 <- which(data$cord   == (sLRC_pair[1]))
  XRC_near_st = which(data$cord %in% (nearby_points_st$cord[near_st]))
  for (pp2 in XRC_near_st) {
    plot <- add_trace(plot,
                      x = c(data$x[pp1], data$x[pp2]),
                      y = c(data$y[pp1], data$y[pp2]),
                      z = c(data$z[pp1], data$z[pp2]),
                      type = "scatter3d",
                      mode = "lines",
                      line = list(color = "blue", width = 3),
                      name = "sLRC",
                      legendgroup = "sLRC",
                      showlegend = F,
                      visible = TRUE)
  }
  pp2 <- which(data$cord   == (sLRC_pair[2]))
  XRC_near_ed = which(data$cord %in% (nearby_points_ed$cord[near_ed]))
  for (pp1 in XRC_near_ed) {
    plot <- add_trace(plot,
                      x = c(data$x[pp1], data$x[pp2]),
                      y = c(data$y[pp1], data$y[pp2]),
                      z = c(data$z[pp1], data$z[pp2]),
                      type = "scatter3d",
                      mode = "lines",
                      line = list(color = "blue", width = 3),
                      name = "sLRC",
                      legendgroup = "sLRC",
                      showlegend = F,
                      visible = TRUE)
  }
  pp1 <- which(data$cord == (sLRC_pair[1]))
  XRC_superl_st = which(data$cord %in% (nearby_points_st$cord[superl_st]))
  for (pp2 in XRC_superl_st) {
    plot <- add_trace(plot,
                      x = c(data$x[pp1], data$x[pp2]),
                      y = c(data$y[pp1], data$y[pp2]),
                      z = c(data$z[pp1], data$z[pp2]),
                      type = "scatter3d",
                      mode = "lines",
                      line = list(color = "red", width = 3),
                      name = "sLRC",
                      legendgroup = "sLRC",
                      showlegend = F,
                      visible = TRUE)
  }
  pp2 <- which(data$cord == (sLRC_pair[2]))
  XRC_superl_ed = which(data$cord %in% (nearby_points_ed$cord[superl_ed]))
  for (pp1 in XRC_superl_ed) {
    plot <- add_trace(plot,
                      x = c(data$x[pp1], data$x[pp2]),
                      y = c(data$y[pp1], data$y[pp2]),
                      z = c(data$z[pp1], data$z[pp2]),
                      type = "scatter3d",
                      mode = "lines",
                      line = list(color = "red", width = 3),
                      name = "sLRC",
                      legendgroup = "sLRC",
                      showlegend = F,
                      visible = TRUE)
  }
  # Add scale bar
  scale_length <- 1000
  min_x <- min(data$x,na.rm = T)
  min_y <- min(data$y,na.rm = T)
  min_z <- min(data$z,na.rm = T)
  
  plot <- add_trace(plot,
                    x = c(min_x, min_x + scale_length),
                    y = c(min_y, min_y),
                    z = c(min_z, min_z),
                    type = "scatter3d",
                    mode = "lines",
                    line = list(color = "red", width = 5),
                    name = "Scale (100 units)",
                    showlegend = F)
  plot <- add_trace(plot,
                    x = c(min_x, min_x),
                    y = c(min_y, min_y + scale_length),
                    z = c(min_z, min_z),
                    type = "scatter3d",
                    mode = "lines",
                    line = list(color = "red", width = 5),
                    name = "Scale (100 units)",
                    showlegend = F)
  plot <- add_trace(plot,
                    x = c(min_x, min_x),
                    y = c(min_y, min_y),
                    z = c(min_z, min_z + scale_length),
                    type = "scatter3d",
                    mode = "lines",
                    line = list(color = "red", width = 5),
                    name = "Scale (100 units)",
                    showlegend = F)
  # Set up layout
  plot <- plot %>% layout(
    scene = list(
      aspectmode = "cube",
      xaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = F, title = ""),
      yaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = F, title = ""),
      zaxis = list(showgrid = FALSE, zeroline = T, showline = T, showticklabels = F, title = "")
    ),
    updatemenus = list(
      list(
        type = "buttons",
        direction = "right",
        xanchor = "center",
        yanchor = "top",
        pad = list(r = 0, t = 10, b = 10),
        x = 0.5,
        y = 1.1,
        buttons = list(
          list(method = "update",
               args = list(visible = list(TRUE, TRUE)),
               label = "All"),
          list(method = "update",
               args = list(visible = list(TRUE, FALSE)),
               label = "Distance Group"),
          list(method = "update",
               args = list(visible = list(FALSE, TRUE)),
               label = "Other Group")
        )
      )
    )
  )
  print(plot)
}


####################################################################################################
# Fig. S9f  ####
####################################################################################################
d3plot(2410, c(33,     189.75))   # example
