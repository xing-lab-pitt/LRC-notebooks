library(ggplot2)
library(plotly)
library(circlize)

# create_coordinate_pairs: create_coordinate_pairs
## chr_length is the length of targeted chromosome, chr2 = 242200000, chr14 = 108000000, chr21 = 46700000; res: resolution, unit of res is kb.
create_coordinate_pairs <- function(chr_length, res) {
  # Generate sequence of coordinates
  positions <- seq(0, chr_length, res * 1000)
  positions_str <- as.character(format(positions, scientific = FALSE))
  n <- length(positions_str)
  
  # Pre-allocate vectors with the exact size needed
  total_pairs <- sum(n - (2:n))
  loc1 <- character(total_pairs)
  loc2 <- character(total_pairs)
  
  # Fill vectors
  idx <- 1
  for (i in 1:(n-2)) {
    # Add progress indicator for long chromosomes
    if(i %% 100 == 0) {
      cat(sprintf("Processing position %d of %d\n", i, n-2))
    }
    
    # Calculate number of pairs for this position
    num_pairs <- n - i - 1
    
    # Assign values directly to pre-allocated vectors
    loc1[idx:(idx+num_pairs-1)] <- positions_str[i]
    loc2[idx:(idx+num_pairs-1)] <- positions_str[(i+2):n]
    
    # Update index
    idx <- idx + num_pairs
  }
  
  # Return as data frame
  return(data.frame(loc1 = loc1, loc2 = loc2))
}
# Hi-C data statistic analysis function
stat_analyze_hic <- function(filename, res, chr_length) {
  # Load and preprocess data
  hic <- read.delim(filename, header = TRUE)
  print(paste("File name:",filename))
  message(sprintf("Initial dimensions: %d x %d", nrow(hic), ncol(hic)))
  
  # Set column names based on file format
  colnames(hic) <- c('chr1','start1','end1','chr2','start2','end2','count','balanced')
  message("Data sample:")
  print(head(hic))
  
  # Filter data 
  valid_indices <- which(abs(hic$start2 - hic$start1) > (res*1000))
  hic2 <- hic[valid_indices, ]
  message(sprintf("After filtering: %d x %d", nrow(hic2), ncol(hic2)))
  
  # Create a data frame directly instead of cbind operations
  formatted_data <- data.frame(
    start1 = format(hic2$start1, scientific = FALSE),
    start2 = format(hic2$start2, scientific = FALSE),
    balanced = hic2$balanced,
    stringsAsFactors = FALSE
  )
  
  # Create coordinate pairs using vectorized operations
  st <- 0
  ed <- max(as.numeric(hic2$end1), as.numeric(hic2$end2))
  positions_str <- format(seq(st, ed, res*1000), scientific = FALSE)
  
  # Create coordinate pairs more efficiently
  n <- length(positions_str)
  pair_indices <- expand.grid(i = 1:(n-2), j = 3:n)
  pair_indices <- subset(pair_indices, j > i)
  
  coord_pairs <- data.frame(
    loc1 = positions_str[pair_indices$i],
    loc2 = positions_str[pair_indices$j],
    stringsAsFactors = FALSE
  )
  
  # Match data using merge instead of match for better performance
  merged_data <- merge(
    coord_pairs,
    formatted_data,
    by.x = c("loc1", "loc2"),
    by.y = c("start1", "start2"),
    all.x = TRUE
  )
  
  # Convert to numeric in one step
  merged_data <- transform(
    merged_data,
    loc1 = as.numeric(loc1),
    loc2 = as.numeric(loc2),
    balanced = as.numeric(balanced)
  )
  
  # Calculate distances and transform data
  dis <- abs(merged_data$loc2 - merged_data$loc1)
  hic_transformed <- data.frame(
    dis = log10(dis),
    interactions = log10(merged_data$balanced),
    counts = merged_data$balanced
  )
  
  # Remove NA values once
  hic_transformed <- na.omit(hic_transformed)
  message(sprintf("Final transformed data: %d x %d", nrow(hic_transformed), ncol(hic_transformed)))
  
  # Group data more efficiently
  tmp_seq <- round(seq(0.1, max(hic_transformed$dis) + 0.2, 0.2), 1)
  gp <- findInterval(hic_transformed$dis, tmp_seq)
  grouped_data <- cbind(hic_transformed, grp = tmp_seq[gp])
  
  
  # Return results as a list
  return(data = grouped_data)
}
# process_hic_data: pre-process Hi-C data
process_hic_data <- function(filename, res, coord_pair) {
  message(paste("Processing file:", filename))
  
  # Read and preprocess data
  hic <- read.delim(filename, header = TRUE)
  message(sprintf("Initial dimensions: %d x %d", nrow(hic), ncol(hic)))
  # Set column names
  colnames(hic) <- c('chr1','start1','end1','chr2','start2','end2','count','balanced')
  
  # Filter data more efficiently using logical operations
  valid_indices <- (abs(hic$start2 - hic$start1) > (res*1000))
  hic_filtered <- hic[valid_indices, ]
  message(sprintf("After filtering: %d x %d", nrow(hic_filtered), ncol(hic_filtered)))
  
  # Create formatted data
  formatted_data <- data.frame(
    start1 = format(hic_filtered$start1, scientific = FALSE),
    start2 = format(hic_filtered$start2, scientific = FALSE),
    count = hic_filtered$count,
    balanced = hic_filtered$balanced,
    stringsAsFactors = FALSE
  )
  
  # Match coordinates more efficiently
  match_key1 <- paste(coord_pair[,1], coord_pair[,2])
  match_key2 <- paste(formatted_data$start1, formatted_data$start2)
  match_indices <- match(match_key1, match_key2)
  
  # Create result data frame
  result <- data.frame(
    loc1 = coord_pair[,1],
    loc2 = coord_pair[,2],
    counts = formatted_data$count[match_indices],
    balance = formatted_data$balanced[match_indices],
    stringsAsFactors = FALSE
  )
  
  # Convert all columns to numeric in one step
  result[] <- lapply(result, function(x) as.numeric(as.character(x)))
  
  return(result)
}
# create_log10_matrix: convert the data frame to log10 matrix
# hic_balanced is a data.frame with at least the loci of paired regions and the balanced interaction strength (from cooler). 
# st: the smallest number of all involved loci; ed: the biggest number of all involved loci; res: resolution of hic. 
# out put is the log10 hic matrix
create_log10_matrix <- function(hic_balanced, st, ed, res) {
  # Create position sequence
  locs <- seq(st, ed, res*1000)
  n_locs <- length(locs)
  
  # Pre-allocate matrix
  hic_matrix <- matrix(NA, ncol = n_locs, nrow = n_locs)
  rownames(hic_matrix) <- locs
  colnames(hic_matrix) <- locs
  
  # Create lookup table for locations (much faster than repeated matches)
  loc_indices <- match(hic_balanced$loc1, locs)
  loc2_indices <- match(hic_balanced$loc2, locs)
  
  # For each unique loc1 value, find all corresponding entries
  unique_loc1 <- unique(hic_balanced$loc1)
  
  for (loc in unique_loc1) {
    # Find the row index for this location
    i <- match(loc, locs)
    if (is.na(i)) next
    
    # Find all entries with this loc1
    entries <- hic_balanced[hic_balanced$loc1 == loc, ]
    
    # Find the column indices for these entries
    j_indices <- match(entries$loc2, locs)
    
    # Assign values to the matrix
    hic_matrix[i, j_indices] <- log10(entries$balance)
  }
  
  # Convert rows to numeric (though this should already be handled by the matrix creation)
  for (i in 1:nrow(hic_matrix)) {
    hic_matrix[i, ] <- as.numeric(hic_matrix[i, ])
  }
  
  return(hic_matrix)
}
# calculate_neighborhood_matrix: Convert log10 matrix to 7×7-neighborhood mean matrix
# hic_matrix is the log10 converted hiC matrix, but other matrix can also be used. nnb is the same nnb used as in calculate_neighborhood_mean(). 
calculate_neighborhood_matrix <- function(hic_matrix, st, ed, res, nnb) {
  # Create position sequence
  locs <- seq(st, ed, res*1000)
  nn <- nrow(hic_matrix)
  
  # Pre-allocate result matrix
  hic_matrix2 <- matrix(NA, nrow = nn, ncol = nn)
  
  # Add progress reporting
  total_iterations <- (nn-4) * (nn-5) / 2
  iterations_done <- 0
  message("Starting neighborhood averaging...")
  progress_interval <- max(1, floor(total_iterations / 10))
  
  # Calculate neighborhood means
  for (ii in 1:(nn-4)) {
    for (jj in (ii+3):(nn-3)) {
      hic_matrix2[ii, jj] <- calculate_neighborhood_mean(hic_matrix, ii, jj, nnb)
      
      # Report progress
      iterations_done <- iterations_done + 1
      if (iterations_done %% progress_interval == 0) {
        percent_done <- round(iterations_done / total_iterations * 100)
        message(sprintf("Processing: %d%% complete", percent_done))
      }
    }
  }
  
  # Set row and column names
  rownames(hic_matrix2) <- locs
  colnames(hic_matrix2) <- locs
  
  message("Neighborhood averaging complete")
  return(hic_matrix2)
}
# calculate_neighborhood_mean: calculate the mean of one 7×7-neighborhood block and substitute the original log10 data. 
# for res = 100 kb, #NA<15, for res = 50 kb, #NA<10; since biologically, higher resolution 
calculate_neighborhood_mean <- function(matrix, i, j, nnb) {
  n <- nrow(matrix)
  neighborhood <- matrix[
    max(1, i-3):min(n, i+3),
    max(1, j-3):min(n, j+3)
  ]
  if (length(which(!is.na(neighborhood)))<nnb) {
    return(NA)
  }
  else{
    return(mean(neighborhood,na.rm = T))
  }
}
# calculate_distance_iqrs: calculate the iqr
calculate_distance_iqrs <- function(matrix, st, ed, res) {
  # Create positions
  locs <- seq(st, ed, res*1000)
  n_locs <- length(locs)
  
  # Pre-allocate vectors for better performance
  total_pairs <- sum((n_locs-1):(n_locs-(n_locs-2)))
  lg_balance <- numeric(total_pairs)
  loc1 <- numeric(total_pairs)
  loc2 <- numeric(total_pairs)
  
  # Fill vectors more efficiently
  idx <- 1
  for (ii in 1:(n_locs-2)) {
    # Number of pairs for this position
    num_pairs <- n_locs - ii - 1
    
    # Extract values directly
    lg_balance[idx:(idx+num_pairs-1)] <- matrix[ii, (ii+2):n_locs]
    loc1[idx:(idx+num_pairs-1)] <- locs[ii]
    loc2[idx:(idx+num_pairs-1)] <- locs[(ii+2):n_locs]
    
    # Update index
    idx <- idx + num_pairs
  }
  
  # Calculate distances
  dis <- loc2 - loc1
  dis1 <- sort(unique(dis))
  
  # Pre-allocate IQRs vector
  IQRs <- numeric(length(lg_balance))
  
  # Process first distance separately to avoid boundary issues
  tmp <- dis <= dis1[2]
  tmpp <- lg_balance[tmp]
  tmp1 <- dis == dis1[1]
  tmpp1 <- lg_balance[tmp1]
  q3 <- quantile(tmpp, 0.75, na.rm = TRUE)
  iqr_value <- IQR(tmpp, na.rm = TRUE)
  # if (iqr_value == 0) iqr_value <- 0# Avoid division by zero
  IQRs[tmp1] <- (tmpp1 - q3) / iqr_value
  
  # Process remaining distances with progress reporting
  n_dis <- length(dis1)
  message("Processing distances:")
  report_interval <- max(1, floor(n_dis/20))
  
  for (ii in 2:(n_dis-1)) {
    # Report progress
    if(ii %% report_interval == 0) {
      message(sprintf("  Processing %d of %d (%.1f%%)", ii, n_dis-1, ii/(n_dis-1)*100))
    }
    
    # Get values within distance window
    tmp <- (dis >= dis1[ii-1]) & (dis <= dis1[ii+1])
    tmpp <- lg_balance[tmp]
    
    # Get values at exact distance
    tmp1 <- dis == dis1[ii]
    tmpp1 <- lg_balance[tmp1]
    
    # Calculate IQR scores
    q3 <- quantile(tmpp, 0.75, na.rm = TRUE)
    iqr_value <- IQR(tmpp, na.rm = TRUE)
    # if (iqr_value == 0) iqr_value <- 0  # Avoid division by zero
    IQRs[tmp1] <- (tmpp1 - q3) / iqr_value
  }
  
  # Create data frame with results
  hic_IQRs <- data.frame(
    loc1 = loc1,
    loc2 = loc2, 
    lg_balance = lg_balance,
    IQRs = IQRs
  )
  
  # Clean up IQRs - replace Inf and NA with 0
  hic_IQRs$IQRs[is.infinite(hic_IQRs$IQRs) | is.na(hic_IQRs$IQRs)] <- 0
  
  return(hic_IQRs)
}
# generate_circle_plot：generate circle plot for chromosome interactions
generate_circle_plot <- function(chr, chr_length, res, filename, color) {
  
  # Calculate chromosome segments
  segments_count <- round((chr_length + 2000000) / 2000000) + 1
  chrom <- (1:segments_count) * 2 - 2
  start <- rep(0, segments_count)
  end <- rep(2000000/(res*1000), segments_count)
  chr_sizes_df <- data.frame(chrom = chrom, start = start, end = end)
  
  # Create location sequence
  locs <- seq(0, (chr_length + 2000000), 2000000)
  
  # Process each sample
  plot_color <- color
  transparent_color <- adjustcolor(color, alpha.f = 0.4)
  
  # Read data file
  message("Processing file: ", filename)
  
  hic_IQRs <- tryCatch({
    read.csv(filename, header = TRUE)
  }, error = function(e) {
    stop(paste("Error reading file:", file_path, "-", e$message))
  })
  
  # Clear and initialize circos plot
  circos.clear()
  circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0, 0.05),
    start.degree = 90,
    gap.degree = 0.1,
    clock.wise = TRUE
  )
  
  circos.initialize(
    factors = chr_sizes_df$chrom,
    xlim = cbind(chr_sizes_df$start, chr_sizes_df$end)
  )
  
  # Add track with labels
  circos.track(
    ylim = c(0, 1), 
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      xlim = CELL_META$xlim
      circos.text(mean(xlim), 0.5, chr, cex = 0.5)
    },
    bg.col = transparent_color, 
    bg.border = TRUE, 
    track.height = 0.1
  )
  
  # Add title
  title(paste(flss, '- chr', chr), cex.main = 1.5)
  
  # Filter interaction data
  distances_stat <-hic_IQRs[,2:5]
  colnames(distances_stat) <- c("loc1", "loc2", "balance", "IQR") # Ensure column names
  
  # Apply distance filter
  long_distance <- distances_stat$loc2 - distances_stat$loc1 > 50000000
  strong_signal <- distances_stat[, 4] > 2
  filtered_data <- distances_stat[long_distance & strong_signal, ]
  
  message("Filtered data dimensions: ", nrow(filtered_data), "x", ncol(filtered_data))
  
  # Draw links between interacting regions
  for (segment_i in 1:(segments_count - 1)) {
    loc_start <- locs[segment_i]
    loc_end <- locs[segment_i + 1]
    loc_seq_start <- seq(loc_start, loc_end, res * 1000)
    
    # Find interactions starting in this segment
    segment_interactions <- filtered_data[
      filtered_data$loc1 >= loc_start & filtered_data$loc1 < loc_end, 
    ]
    
    if (nrow(segment_interactions) > 0) {
      # For each possible target segment
      for (segment_j in (segment_i + 1):(segments_count - 1)) {
        target_start <- locs[segment_j]
        target_end <- locs[segment_j + 1]
        loc_seq_end <- seq(target_start, target_end, res * 1000)
        
        # Find interactions ending in the target segment
        target_interactions <- segment_interactions[
          segment_interactions$loc2 >= target_start & 
            segment_interactions$loc2 < target_end,
        ]
        
        if (nrow(target_interactions) > 0) {
          # Create matrix for vectorized operations
          for (pos_start_idx in 1:length(loc_seq_start)) {
            pos_start <- loc_seq_start[pos_start_idx]
            start_interactions <- target_interactions[target_interactions$loc1 == pos_start, ]
            
            if (nrow(start_interactions) > 0) {
              for (pos_end_idx in 1:length(loc_seq_end)) {
                pos_end <- loc_seq_end[pos_end_idx]
                interaction <- start_interactions[start_interactions$loc2 == pos_end, ]
                
                if (nrow(interaction) > 0 && interaction$IQR != 0) {
                  # Calculate link thickness based on interaction strength
                  link_thickness <- min(interaction$IQR / 20, 2) # Cap thickness
                  
                  # Draw the link
                  circos.link(
                    segment_i * 2 - 2, pos_start_idx, 
                    segment_j * 2 - 2, pos_end_idx, 
                    col = color, 
                    lwd = link_thickness
                  )
                }
              }
            }
          }
        }
      }
    }
  }
}
is_outlier <- function(x,coef) {
  return(x > quantile(x, 0.75) + coef * IQR(x))
}
iqr =2

chr2 = 242200000
chr14 = 108000000
chr21 = 46700000
# input, use IMR90 chr14 as an example ####
filenms = c(’IMR90‘)
chr_length = chr14
chr = 14
res = 100

for (filename in filenms) {
  print(res)
  
  coord_pair <- create_coordinate_pairs(chr_length, res) # combination of region pairs
  print(dim(coord_pair))
  
  filename = paste0(filenms,'.chr',chr,'.output.',res,'k.txt')
  print(filename)
  grouped_data <- stat_analyze_hic(filename, res, chr_length)
  
  ## Create plot ####
  # scatter boxplot 
  scatter_plot <- ggplot(grouped_data, aes(x = dis, y = interactions, color = interactions)) +
    geom_point(shape = 20, size = 0.5, alpha = 0.7) +
    scale_color_gradientn(colors = c("goldenrod3", "slateblue")) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_line(colour = "black"),
      axis.ticks.length = unit(0.2, "cm"),
    ) +
    labs(
      title = paste('chr', chr, 'resolution:', res, 'kb'),
      x = "Genomic distances (log10)",
      y = "Interactions (log10)"
    )
  print(scatter_plot)
  # Calculate medians for boxplot
  medians <- tapply(grouped_data$interactions, grouped_data$grp, median, na.rm = TRUE)
  # Create boxplot 
  box_plot <- ggplot(grouped_data, aes(x = factor(grp), y = interactions, group = grp)) +
    stat_boxplot(geom = 'errorbar', coef = iqr, width = 0.5) +
    geom_boxplot(
      fill = "burlywood1", 
      outlier.colour = "black",
      outlier.size = 0.5,
      coef = iqr
    ) +
    stat_summary(
      fun = median, 
      geom = "line", 
      aes(group = 1), 
      color = "blue", 
      size = 1
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    labs(
      title = basename(filename),
      x = "Genomic distance (MB)",
      y = "Interactions (log10)"
    )
  print(box_plot)
  
  ## CLOD pipeline ####
  ##### preprocess hi-c data  #####
  hic3 = process_hic_data(filename, res, coord_pair)
  
  ##### convert the data frame to log10 matrix  #####
  st <- 0
  ed <- max(as.numeric(coord_pair$loc1), as.numeric(coord_pair$loc2))
  locs <- seq(st, ed, res*1000)

  hic_matrix = create_log10_matrix(hic3,st,ed,res)
  # 3D plot the log10 matrix
  tmp1 = which(locs>20000000&locs<30000000)
  tmp2 = which(locs>100000000&locs<107000000)
  hic_matrix_tmp= hic_matrix[tmp1,tmp2]
  hic_matrix_tmp[which(is.na(hic_matrix_tmp),arr.ind = T)] = min(hic_matrix,na.rm = T)
  plot_ly() %>% add_surface(x = ~locs[tmp2], y = ~locs[tmp1], z = ~hic_matrix_tmp)%>%
    layout(
      scene = list(
        xaxis = list(
          showgrid = FALSE,
          zeroline = FALSE,
          showticklabels = FALSE,
          title = ""
        ),
        yaxis = list(
          showgrid = FALSE,
          zeroline = FALSE,
          showticklabels = FALSE,
          title = ""
        ),
        zaxis = list(
          showgrid = FALSE,
          zeroline = FALSE,
          showticklabels = FALSE,
          title = ""
        )
      )
    )
  
  #### convert log10 matrix to 7×7-neighborhood mean matrix #####
  nnb = 15
  hic_matrix2 = calculate_neighborhood_matrix(hic_matrix, st, ed, res, nnb)
  # plot
  hic_matrix_tmp = hic_matrix2[tmp1,tmp2]
  rownames(hic_matrix_tmp) = locs[tmp1]/1000000
  colnames(hic_matrix_tmp) = locs[tmp2]/1000000
  hic_matrix_tmp[which(is.na(hic_matrix_tmp),arr.ind = T)] = min(hic_matrix,na.rm = T)
  plot_ly() %>% add_surface(x = ~locs[tmp2], y = ~locs[tmp1], z = ~hic_matrix_tmp)%>% colorbar(title = "log10")
  
  #### convert neighborhood mean matrix to IQR matrix by smoothed data #####
  hic_IQRs  = calculate_distance_iqrs(hic_matrix2, st, ed, res)
  # quick check
  tmp = which(hic_IQRs[,1]>20000000&hic_IQRs[,1]<30000000&hic_IQRs[,2]>103000000&hic_IQRs[,2]<107000000)
  View(hic_IQRs[tmp,])
  
 # ## save IQR>2 #####
 #     tmp = which(hic_IQRs$IQRs>2)
 #     print(length(tmp))
 #     hic_IQRs2 = hic_IQRs[tmp,]
 #     hic_IQRs2 = cbind(rep(paste0('chr',chr),length(tmp)),hic_IQRs2)
 #     hic_IQRs2[,4] = round(hic_IQRs2[,4],3)
 #     hic_IQRs2[,5] = round(hic_IQRs2[,5],3)
 #     colnames(hic_IQRs2)[1] = 'chr'
 #     print(dim(hic_IQRs2))
 #     head(hic_IQRs2)
 #     write.csv(hic_IQRs2, paste0("outliers_",flss,'-',res,"k-chr",chr,".csv"), row.names = FALSE)
  }

# plot the linkage maps                     
filename = paste0("outliers_",flss,'-',res,"k-chr",chr,".csv")
color = '#FDB46F'
generate_circle_plot(chr, chr_length, res, filename, color)                 
