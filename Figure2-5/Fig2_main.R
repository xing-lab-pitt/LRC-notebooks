library(ggplot2)
library(plotly)
library(dplyr)
library(circlize)
library(htmlwidgets)

IMR_2 <- read.delim("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/Verify/chromosome2.tsv")
print(colnames(IMR_2))
print(dim(IMR_2))
tmp = c(1,2,3,198,546,556,1259,2368,2862)
tmp = which(IMR_2$Chromosome.copy.number%in%tmp)
IMR_2222 = IMR_2[tmp,]
write.csv(IMR_2222, paste0("DNAFISH_IMR90_chr2_lite.csv"), row.names = FALSE)

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

fils = seq(100,3000,100)
ff = 1
distances_2 = read.table(paste0("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/Verify/distances_2-",fils[ff],".txt"),skip = 1,header = FALSE)
dim(distances_2)
distances_2 = round(distances_2[,1],0)[-1]
print(length(distances_2))

spots2 = unique(Coordinate)
distances_tmp = matrix(0,ncol = length(spots2),nrow = length(spots2))
colnames(distances_tmp) = spots2
rownames(distances_tmp) = spots2
distances_tmp[lower.tri(distances_tmp)] <- distances_2[1:(length(distances_2)/100)]

st = 0
ed = 0
dst = 0
for (jj in 1:length(colnames(distances_tmp))) {
  st = c(st,rep(colnames(distances_tmp)[jj],nrow(distances_tmp)-jj))
  ed = c(ed, rownames(distances_tmp)[-c(1:jj)])
  # dst = c(dst,distances_tmp[(jj+1):nrow(distances_tmp),jj])
}
st = st[-1]
ed = ed[-1]
# dst = dst[-1]
lner = as.numeric(ed)-as.numeric(st)

tt1 = which(lner>=100000000)
tt2 = which((lner>50000000)&(lner<100000000))

# genomic distances (l) between pair of bins####
## Fig. S1B ####
IMR_21 <- read.delim("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/Verify/chromosome21.tsv")
print(head(IMR_21))
print(colnames(IMR_21))
print(dim(IMR_21))

Co = IMR_21$Genomic.coordinate
Coordinate = rep(0,nrow = length(Co))
for (ii in 1:length(spots)) {
  Co_tmp = which(Co==spots[ii])
  tmp = strsplit(spots[ii],split =':')[[1]]
  tmp = strsplit(tmp[2],split ='-')[[1]]
  Coordinate[Co_tmp] =(as.numeric(tmp[1]))
}
IMR_21$Coordinate = Coordinate

# chr21 32.45 MB - 33.35 MB, Fig. S1B
tmp = which(Coordinate>32450000&Coordinate<33350000)
IMR_211 = IMR_21[tmp,]

neighbor50k2 = 0
for (ii in 2:nrow(IMR_211)) {
  if (ii%%100000==0) {
    print(ii)
  }
  if ((IMR_211$Coordinate[ii]-IMR_211$Coordinate[ii-1])==50000) {
    tmp = ((IMR_211$Z.nm.[ii-1]-IMR_211$Z.nm.[ii])^2+(IMR_211$X.nm.[ii-1]-IMR_211$X.nm.[ii])^2+(IMR_211$Y.nm.[ii-1]-IMR_211$Y.nm.[ii])^2)^0.5
    neighbor50k2 = c(neighbor50k2,as.integer(tmp))
  }
}
neighbor50k2 = neighbor50k2[-1]
print(median(neighbor50k2,na.rm = T))

A = as.vector(neighbor50k2)
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
print(quartiles)
# Create the histogram with density curve
hist(A, probability = TRUE, breaks = 400, xlim = c(1,600), 
     main = "chr21:32.45 MB-33.35 MB", 
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c( "Median"),
       col = c( "brown3"),
       lty = c( 2), lwd = 2)

neighbor50k = 0
for (ii in 2:nrow(IMR_21)) {
  if (ii%%100000==0) {
    print(ii)
  }
  if (IMR_21$Coordinate[ii]-IMR_21$Coordinate[ii-1]==50000) {
    tmp = ((IMR_21$Z.nm.[ii-1]-IMR_21$Z.nm.[ii])^2+(IMR_21$X.nm.[ii-1]-IMR_21$X.nm.[ii])^2+(IMR_21$Y.nm.[ii-1]-IMR_21$Y.nm.[ii])^2)^0.5
    # if (is.na(tmp)) {
    #   print(ii)
    # }
    neighbor50k = c(neighbor50k,as.integer(tmp))
  }
}
neighbor50k = neighbor50k[-1]
print(median(neighbor50k,na.rm = T))

A = as.vector(neighbor50k)
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
print(quartiles)
# Create the histogram with density curve
hist(A, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "chr21", 
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

IMR_2 <- read.delim("/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/Verify/chromosome2.tsv")
print(colnames(IMR_2))
print(dim(IMR_2))

cells = unique(IMR_21$Chromosome.copy.number)
spots = unique(IMR_21$Genomic.coordinate)
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

## the 3D distances between pairs of neibhor 250k ####
neighbor250k = matrix(0,nrow = 1000, ncol = 3000)
aa = 1
for (ll in 1:5) {
  tt = which((spots2>(ll-1)*50000000)&(spots2<(ll*50000000+1)))
  distances_2m_tt = read.table(paste0('/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/distances_chr2R',(ll-1)*50,'-',(ll*50),'M_3000.txt'))
  print(dim(distances_2m_tt))
  colnames(distances_2m_tt) = distances_2m_tt[1,]
  distances_2m_tt = t(distances_2m_tt[-1,])
  print(dim(distances_2m_tt))
  if((dim(distances_2m_tt)[2]/3000)!=dim(distances_2m_tt)[1]){
    print('WRONG!!!')
    break
  }
  bb = aa+dim(distances_2m_tt)[1]-2
  for (cells in 1:3000) {
    start_col <- (cells - 1) * dim(distances_2m_tt)[1] + 1
    end_col <- (cells) * dim(distances_2m_tt)[1]
    xx<- distances_2m_tt[, start_col:end_col]
    neighbor250k[aa:bb,cells] = xx[row(xx) == col(xx) - 1]
  }
  aa = bb+1
  print(c(bb,aa))
}
print(aa)
neighbor250k = neighbor250k[1:(bb),]
# # double check
# neighbor250k2 = 0
# for (ii in 2:nrow(IMR_2)) {
#   if (ii%%100000==0) {
#     print(ii)
#   }
#   if (IMR_2$Coordinate[ii]-IMR_2$Coordinate[ii-1]==250000) {
#     tmp = ((IMR_2$Z.nm.[ii-1]-IMR_2$Z.nm.[ii])^2+(IMR_2$X.nm.[ii-1]-IMR_2$X.nm.[ii])^2+(IMR_2$Y.nm.[ii-1]-IMR_2$Y.nm.[ii])^2)^0.5
#     neighbor250k2 = c(neighbor250k2,as.integer(tmp))
#   }
# }
# neighbor250k2 = neighbor250k2[-1]

print(median(neighbor250k,na.rm = T))
A = as.vector(neighbor250k)
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
# Create the histogram with density curve
hist(A, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "chr2", 
     xlab = "3D dist. (nm)", col = "lightblue", border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

## compared with distances 50-100 and >100 MB #####
distances_2m100 = read.table(paste0('/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/distances_chr2L100-3000.txt'),skip = 1)
distances_2m1000 = data.frame(distances_2m100)
colnames(distances_2m1000) = c('distance','homolog')

distances_2m50 = read.table(paste0('/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/distances_chr2L50-3000.txt'),skip = 1)
distances_2m500 = data.frame(distances_2m50)
colnames(distances_2m500) = c('distance','homolog')


# Combine them
color_A <- "orange"  # Soft blue
color_B <- "#9b59b6"  # Coral red
color_C <- "darkolivegreen3"  # Coral red
A = as.numeric(distances_2m1000$distance)
# calculate l between 100 MB bins
quartiles <- quantile(A, probs = c(0.25, 0.5, 0.75),na.rm = T)
hist(A, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "3D dist. between pair of loci >100 MB,chr2", 
     xlab = "3D dist. (nm)", col = color_A, border = "white")
box(lwd = 3)
lines(density(A,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

B = as.numeric(distances_2m500$distance)
# calculate l between 50 MB bins
quartiles <- quantile(B, probs = c(0.25, 0.5, 0.75),na.rm = T)
# Create the histogram with density curve
hist(B, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "3D dist. between pair of loci 50-100 MB,chr2", 
     xlab = "3D dist. (nm)", col = color_B, border = "white")
box(lwd = 3)
lines(density(B,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)

C = as.vector(neighbor250k)
# calculate l between neibhor bins
quartiles <- quantile(C, probs = c(0.25, 0.5, 0.75),na.rm = T)
# Create the histogram with density curve
hist(C, probability = TRUE, breaks = 200, xlim = c(1,2500), 
     main = "3D dist. between neighbor loci,chr2", 
     xlab = "3D dist. (nm)", col = color_C, border = "white")
box(lwd = 3)
lines(density(C,na.rm = T), col ="darkblue", lwd = 2)
abline(v = quartiles[2], col ="brown3", lty = 2, lwd = 2)  
legend("topright", 
       legend = c("Median"),
       col = c("brown3"),
       lty = c(2), lwd = 2)
# Calculate the range for both A and B
# x_range <- range(c(A, B),na.rm = T)
# Calculate densities
dens_A <- density(A,na.rm = T)
dens_B <- density(B,na.rm = T)
dens_C <- density(C,na.rm = T)
# Calculate the y-axis limit
y_max <- max(c(dens_A$y, dens_B$y,dens_C$y))
# Create an empty plot
plot(0, 0, type = "n", xlim = c(0,6000), ylim = c(0, y_max * 1.1),
     main = "Overlaid Histograms with Density Curves",
     xlab = "3D dist. (nm)", ylab = "Density")
box(lwd = 3)
# Add subtle grid lines for better readability
grid(col = "gray90", lty = "dotted")
hist(A, freq = FALSE, col = adjustcolor(color_A, alpha.f = 0.5), add = TRUE, border = NA,breaks = 200)
hist(B, freq = FALSE, col = adjustcolor(color_B, alpha.f = 0.5), add = TRUE, border = NA,breaks = 200)
hist(C, freq = FALSE, col = adjustcolor(color_C, alpha.f = 0.5), add = TRUE, border = NA,breaks = 200)
lines(dens_A, col = color_A, lwd = 2)
lines(dens_B, col = color_B, lwd = 2)
lines(dens_C, col = color_C, lwd = 2)
abline(v = quantile(C,c(0.25,0.5,0.75),na.rm=T)[1], col = "darkblue", lty = 2, lwd = 2) 
# Add a legend
legend("topright", legend = c("bin pairs >100M","bin pairs 50-100M", "bin pairs = 250K","1st QTL of the neibhor"),
       col = c(adjustcolor(color_A, alpha.f = 0.5), adjustcolor(color_B, alpha.f = 0.5), adjustcolor(color_C, alpha.f = 0.5),adjustcolor("darkblue", alpha.f = 0.5)),
       lty = c(1,1,1,2), lwd = 2,
       border = NA, bty = "n")

# Figure 2C, relationship between convex hull volume v.s. LRC counts ####
convex_hull_volume <- function(points) {
  library(geometry)
  hull <- convhulln(points, output.options = "FA")
  return(hull$vol)
}

distances_2m100 = read.table(paste0('/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/distances_chr2L100-3000.txt'),skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

chr2 = 242200000
rsl = 50
d3li = 299

## distribution of distances of loci pairs ####
conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

chv = rep(NA,length(cells))
for (ii in 1:length(cells)) {
  if (ii%%500==0) {
    print(ii)
  }
  tmp = which(IMR_2$Chromosome.copy.number ==ii)
  set = as.matrix(IMR_2[tmp,1:3])
  set = na.omit(set)
  
  vol <- convex_hull_volume(set)
  chv[ii] = vol
}
rcr = cbind(chv, colSums(conformation_matrix))
colnames(rcr) = c('chv','numbLRC')
plot(rcr[,1],rcr[,2],col = 'darkgrey',pch = '.',cex = 2)
plot(log10(rcr[,1]),log10(rcr[,2]),col = 'darkgrey',pch = '.')
abline(v = 2.5)

# contours plot
data <- data.frame(x = log10(rcr[,1]), y = log10(rcr[,2]))
data %>% ggplot(aes(x, y))+
  geom_point(colour = "brown", size = 0.1)+
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "darkblue",alpha =0.3,bins = 20)+
  scale_fill_distiller(palette = "YlGn", direction = 1) +
  theme_classic()+
  labs(title = "relationship between # of LRC and compaction",
       x = "convex hull volume",
       y = "# LRC/chr.")

## number of count of LRC in individual chromosome ####
hist(log10(colSums(conformation_matrix)), probability = TRUE, breaks = 200,
     main = "count of LRC in individual chr.", 
     xlab = "count of LRC", col = "lightgreen", border = "white",xaxt = "n")  # Suppress default x-axis

##  distribution of individual LRC ####
rmv = which(colSums(conformation_matrix)>10000)
conformation_matrix = conformation_matrix[,-rmv]
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))
hist(distances_stat2[,4])

tgt = which(rowSums(conformation_matrix)>24)
print(length(tgt))
data_df =distances_matrix[tgt,]
data_df2 = data.frame(value=as.vector(data_df), pairs =rep(1:nrow(data_df)+length(tgt),3000))
rdm = 1:nrow(distances_stat2)
rdm = rdm[-tgt]
rdm = sample(rdm,length(tgt))
data_rdm =distances_matrix[rdm,]
data_rdm2 = data.frame(value=as.vector(data_rdm), pairs =rep(1:nrow(data_rdm),3000))
colors = c(rep('random',nrow(data_rdm2)),rep('targets',nrow(data_df2)))
data_df3 = rbind(data_rdm2,data_df2)
data_df3$colors = factor(colors,levels = c('targets','random'))
print(length(tgt))

custom_colors <- c("targets" = "purple", "random" = "lightblue")

# ggplot(data_df3, aes(x = log2(value), group = pairs,color = colors)) +
#   geom_line(stat = "density", size = 0.1, alpha = 0.5) +
#   scale_color_manual(values = custom_colors) +
#   theme_minimal() +
#   labs(title = "Distribution of 1000 Vectors",
#        x = "Value",
#        y = "Density") +
#   xlim(5,13) +  # Limit x-axis from 1 to 100
#   geom_vline(xintercept = log2(299), linetype = "dashed", color = "red")

data_df3t = data_df3[which(data_df3$colors=='random'),]
ggplot(data_df3t, aes(x = (value), group = pairs,color = colors)) +
  geom_line(stat = "density", size = 0.1, alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(title = "3D dist. in each chr.",
       x = "3D dist. (nm)",
       y = "Density") +
  xlim(-10,6000) +  # Limit x-axis from 1 to 100
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red")+
  geom_vline(xintercept = 299, linetype = "dashed", color = "red")
ggsave("random.png",width = 947/72, height = 405/72, units = "in", dpi = 300,  bg = "transparent")

data_df3t = data_df3[which(data_df3$colors=='targets'),]
ggplot(data_df3t, aes(x = log2(value), group = pairs,color = colors)) +
  geom_line(stat = "density", size = 0.1, alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(title = "3D dist. in each chr.",
       x = "3D dist. (nm)",
       y = "Density") +
  xlim(5,13) +  # Limit x-axis from 1 to 100
  geom_vline(xintercept = 5.1, linetype = "dashed", color = "red")
ggsave("targets.png",width = 947/72, height = 405/72, units = "in", dpi = 72,  bg = "transparent")

# Figure 2E, distribution of number of cells that a LRC exist in #####
distances_2m100 = read.table(paste0('/Volumes/InMoLab/Chr14_paper/BioInfo_Analysis/distances_chr2L100-3000.txt'),skip = 1)
distances_2m = data.frame(distances_2m100)
colnames(distances_2m) = c('distance','homolog')
cells = unique(distances_2m[,2])
distances_matrix = matrix(as.numeric(distances_2m$distance),nrow = length(tt1),ncol = length(cells))

conformation_matrix = matrix(0,nrow = length(tt1),ncol = length(cells))
tmp = which((distances_matrix<d3li),arr.ind = T)
conformation_matrix[tmp] = 1
print(dim(conformation_matrix))

distances_stat2 = cbind(as.numeric(as.character(st[tt1])),as.numeric(as.character(ed[tt1])),as.numeric(as.character(lner[tt1])),as.numeric(as.character(rowSums(conformation_matrix))))

data = rowSums(conformation_matrix)
a = boxplot(data,notch = T,col = 'orange',outpch = 19,outcex = 0.5,outcol = 'darkorange3',horizontal = T)
hist(data, breaks = 200, prob = TRUE, col = "orange", xlab = "Values", main = " ")
box(lwd = 3)
abline(v = a$stats[5],lty = 3, col = 'orchid', lwd = 2)

# interactive plots, linkage maps, 3D, networks ####
