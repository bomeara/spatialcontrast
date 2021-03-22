#' @param latlon a data.frame with taxon, taxon_full, lon, lat
#' @param minlat all points will be greater than or equal to this latitude
#' @param maxlat all points will be less than this latitude
#' @param minlon all points will be greater than or equal to this longitude
#' @param maxlon all points will be less than this longitude
subset_cell <- function(latlon, minlat=-90, maxlat=90, minlon=-180, maxlon=180) {
	return(subset(latlon, latlon$lat>=minlat & latlon$lat<maxlat & latlon$lon>=minlon & latlon$lon < maxlon))
}



present_in_cell <- function(phy, subsetted_cell) {
	trait <- rep(0, ape::Ntip(phy))
	names(trait) <- phy$tip.label
	trait[names(trait) %in% subsetted_cell$taxon] <- 1
	return(trait)
}

present_in_cluster <- function(phy, latlon, cluster_id) {
	trait <- rep(0, ape::Ntip(phy))
	names(trait) <- phy$tip.label
	latlon_sub <- subset(latlon, cluster==cluster_id)
	trait[names(trait) %in% latlon_sub$taxon] <- 1
	return(trait)
}



cluster_latlon <- function(latlon, mean_points_per_cluster=10000, nclusters=NULL) {
	latlon_only <- latlon[,c("lat", "lon")]
	if(is.null(nclusters)){
		 nclusters<-floor(nrow(latlon)/mean_points_per_cluster)
	}
	clusters <- stats::kmeans(latlon_only, centers=nclusters, nstart=20, iter.max=100)
	latlon$cluster <- clusters$cluster
	latlon$cluster_lat <- clusters$centers[latlon$cluster,1]
	latlon$cluster_lon <- clusters$centers[latlon$cluster,2]
	return(latlon)
}


compare_clusters <- function(phy, latlon, tiprates, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
	coarse_rates <- list()
	fine_rates <- list()
	all_results <- as.list(rep(NA, max(latlon$cluster)*max(latlon$cluster) ))
	dim(all_results) <- c(max(latlon$cluster), max(latlon$cluster))
	#rate_matrix <- matrix(NA, nrow=max(latlon$cluster), ncol=max(latlon$cluster))
	for(focal in sequence(max(latlon$cluster))) {
		for (target in sequence(max(latlon$cluster))) {
			focal_presence <- present_in_cluster(phy, latlon, focal)
			target_presence <- present_in_cluster(phy, latlon, target)
			total_presence <- focal_presence + target_presence
			taxa_to_keep <- names(total_presence)[which(total_presence==1)]
			if(length(taxa_to_keep)>=2) {
				phy_pruned <- ape::keep.tip(phy, taxa_to_keep)
				target_pruned <- target_presence[taxa_to_keep]
				sink("/dev/null")
				comparisons <- sanitize_comparisons(sisters::sis_format_comparison(sisters=sisters::sis_get_sisters(phy_pruned, ncores=1), trait=target_pruned, phy=phy_pruned))
				sink()
				global_df <- data.frame()
				if(length(comparisons$node)>=mincomparisons) {
					heights <- phytools::nodeHeights(phy_pruned)
					depths <- max(heights) - heights
					for (sister_index in seq_along(comparisons$node)) {
						focal.taxa <- comparisons$taxa.0[sister_index][[1]]
						target.taxa <- comparisons$taxa.1[sister_index][[1]]
						local_df <- data.frame(node=comparisons$node[sister_index], depth=depths[phy_pruned$edge==comparisons$node[sister_index]][1], ntax.focal=comparisons$ntax.0[sister_index], ntax.target=comparisons$ntax.1[sister_index])
						if(local_df$depth<=maxdepth) {
							pairwise_result <- t(compare_two_clades(focal.taxa, target.taxa, tiprates=tiprates, phy_pruned=phy_pruned, rates=rates, ndraws=ndraws))
							if(!is.na(pairwise_result[1])) {
								local_df <- cbind(local_df,pairwise_result)
								if(sister_index==1) {
									global_df <-local_df
								} else {
									global_df <- rbind(global_df, local_df)
								}
							}
						}
					}
					coarse_df <- apply(subset(global_df, select=c(paste0(rates, "_mean_target_minus_focal"))), 2, mean)
					#coarse_df["npairs"] <- nrow(coarse_df)
					print(c(focal, target))
					print(coarse_df)

				}
				#print(coarse_df)
				all_results[[focal, target]] <- global_df

			}
		}
	}
	return(all_results)
}


summarize_cluster_results <- function(results, latlon, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
	rate_list <- list(rep(NA, length(rates)))
	sign_test_p <- list(rep(NA, length(rates)))
	t_test_p<- list(rep(NA, length(rates)))
	t_test_estimate<- list(rep(NA, length(rates)))
	t_test_lower_ci <- list(rep(NA, length(rates)))
	t_test_upper_ci <- list(rep(NA, length(rates)))
	n_comparisons_list <- list(rep(NA, length(rates)))


	for (rate_index in seq_along(rates)) {
		rate_matrix <- matrix(NA, nrow=nrow(results), ncol=ncol(results))
		sign_test_matrix <- matrix(NA, nrow=nrow(results), ncol=ncol(results))
		t_test_matrix <- matrix(NA, nrow=nrow(results), ncol=ncol(results))
		t_test_estimate_matrix <- matrix(NA, nrow=nrow(results), ncol=ncol(results))
		t_test_lower_ci_matrix <- matrix(NA, nrow=nrow(results), ncol=ncol(results))
		t_test_upper_ci_matrix <- matrix(NA, nrow=nrow(results), ncol=ncol(results))
		n_comparisons_matrix <- matrix(NA, nrow=nrow(results), ncol=ncol(results))
		for (row_index in sequence(nrow(results))) {
			for (col_index in sequence(ncol(results))) {
				global_df <- NULL
				try(global_df <- results[[row_index,col_index]], silent=TRUE)
				if(!is.null(global_df)) {
					count_positive <- count_total <- raw_values <- t_test_result <- NA
					coarse_df <- NULL
					
					try(coarse_df <- apply(subset(global_df, select=c(paste0(rates[rate_index], "_mean_target_minus_focal"))), 2, mean), silent=TRUE)
					try(count_positive <- sum(sign(results[[row_index,col_index]][,paste0(rates[rate_index], "_mean_target_minus_focal")])>0), silent=TRUE)
					try(count_total <- nrow(global_df), silent=TRUE)
					try(raw_values <- results[[row_index,col_index]][,paste0(rates[rate_index], "_mean_target_minus_focal")], silent=TRUE)
					try(sign_test_matrix[row_index, col_index] <- stats::binom.test(x=count_positive, n=count_total, p=0.5, alternative="two.sided")$p.value, silent=TRUE)
					try(t_test_result <- stats::t.test(raw_values, alternative="two.sided"), silent=TRUE)
					try(t_test_matrix[row_index, col_index] <- t_test_result$p.value, silent=TRUE)
					try(t_test_estimate_matrix[row_index, col_index] <- t_test_result$estimate, silent=TRUE)
					try(t_test_lower_ci_matrix[row_index, col_index] <- t_test_result$conf.int[1], silent=TRUE)
					try(t_test_upper_ci_matrix[row_index, col_index] <- t_test_result$conf.int[2], silent=TRUE)
					try(rate_matrix[row_index, col_index]<- unname(coarse_df[paste0(rates[rate_index], "_mean_target_minus_focal")]), silent=TRUE)
					try(n_comparisons_matrix[row_index, col_index] <- nrow(global_df), silent=TRUE)
				}
			}
		}
		rate_list[[rate_index]] <- rate_matrix
		sign_test_p[[rate_index]] <- sign_test_matrix
		t_test_p[[rate_index]] <- t_test_matrix
		t_test_estimate[[rate_index]] <- t_test_estimate_matrix
		t_test_lower_ci[[rate_index]] <- t_test_lower_ci_matrix
		t_test_upper_ci[[rate_index]] <- t_test_upper_ci_matrix
		n_comparisons_list[[rate_index]] <- n_comparisons_matrix


	}
	names(rate_list) <- rates
	names(sign_test_p) <- rates
	names(t_test_p) <- rates
	names(t_test_estimate) <- rates
	names(t_test_lower_ci) <- rates
	names(t_test_upper_ci) <- rates
	names(n_comparisons_list) <- rates

	return(list(rate_difference=rate_list, sign_test_p=sign_test_p, t_test_p=t_test_p, t_test_estimate=t_test_estimate, t_test_lower_ci=t_test_lower_ci, t_test_upper_ci=t_test_upper_ci,n_comparisons=n_comparisons_list))
}

sanitize_comparisons <- function(comparisons) {
	left.unique <- unlist(comparisons$left.unique)	
	right.unique <- unlist(comparisons$right.unique)	
	traitsums <- apply(cbind(left.unique, right.unique), 1, sum) # want those with just state 1 in one clade, state 0 in its sister
	traitsums[is.na(traitsums)] <- 1e10 # some have a mixture on one side, so their unique is NA, so make it a number that won't match
	comparisons <- comparisons[traitsums==1,] # 0 + 1 = 1. Having zero on both sides leads to 0, 1's on both leads to 2, and a mixture on at least one side leads to NA so 1e10 above
	if(nrow(comparisons)>0) {
		comparisons$left.trait <- unlist(comparisons$left.unique)	
		comparisons$right.trait <- unlist(comparisons$right.unique)	
		comparisons$taxa.0 <- NA
		comparisons$taxa.1 <- NA
		comparisons$ntax.0 <- NA
		comparisons$ntax.1 <- NA
		comparisons$taxa.0[comparisons$left.trait==0] <- comparisons$left[comparisons$left.trait==0]
		comparisons$taxa.0[comparisons$right.trait==0] <- comparisons$right[comparisons$right.trait==0]
		comparisons$taxa.1[comparisons$left.trait==1] <- comparisons$left[comparisons$left.trait==1]
		comparisons$taxa.1[comparisons$right.trait==1] <- comparisons$right[comparisons$right.trait==1]
		comparisons$ntax.0 <- sapply(comparisons$taxa.0, length)
		comparisons$ntax.1 <- sapply(comparisons$taxa.1, length)
	}
	return(comparisons)
}


#' @param maxdepth How far back in time to make comparisons
#' @examples 
#' data(magnoliiadae)
#' compare_cells_fine(phy_clean, latlon, tiprates, minlat_focal=-10, maxlat_focal=10, minlon_focal=-100, maxlon_focal=-50, minlat_target=10, maxlat_target=20, minlon_target=-100, maxlon_target=-50)
compare_cells_fine <- function(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
	focal_presence <- present_in_cell(phy, subset_cell(latlon, minlat=minlat_focal, maxlat=maxlat_focal, minlon=minlon_focal, maxlon=maxlon_focal))
	target_presence <- present_in_cell(phy, subset_cell(latlon, minlat=minlat_target, maxlat=maxlat_target, minlon=minlon_target, maxlon=maxlon_target))
	total_presence <- focal_presence + target_presence
	taxa_to_keep <- names(total_presence)[which(total_presence==1)]
	if(length(taxa_to_keep)<2) {
		warning("Not enough taxa to compare these cells")
		means <- as.data.frame(t(rep(NA, length(rates))))
		colnames(means) <- c(paste0(rates, "_mean_target_minus_focal"))
		sds <- as.data.frame(t(rep(NA, length(rates))))
		colnames(sds) <- c(paste0(rates, "_mean_target_minus_focal"))		
		return(cbind(data.frame(node=NA, depth=NA, ntax.left=0, ntax.right=0), means, sds))
	}
	global_df <- data.frame()
	phy_pruned <- ape::keep.tip(phy, taxa_to_keep)
	target_pruned <- target_presence[taxa_to_keep]
	sink("/dev/null")
	comparisons <- sanitize_comparisons(sisters::sis_format_comparison(sisters=sisters::sis_get_sisters(phy_pruned, ncores=1), trait=target_pruned, phy=phy_pruned))
	sink()
	if(length(comparisons$node)<mincomparisons) {
		warning("Not enough sister pairs to compare these cells")
		means <- as.data.frame(t(rep(NA, length(rates))))
		colnames(means) <- c(paste0(rates, "_mean_target_minus_focal"))
		sds <- as.data.frame(t(rep(NA, length(rates))))
		colnames(sds) <- c(paste0(rates, "_mean_target_minus_focal"))		
		return(cbind(data.frame(node=NA, depth=NA, ntax.left=0, ntax.right=0), means, sds))
	}
	heights <- phytools::nodeHeights(phy_pruned)
	depths <- max(heights) - heights
	#cat("\n")
	#pb = utils::txtProgressBar(min = 0, max = length(comparisons$node), initial = 0) 
	for (sister_index in seq_along(comparisons$node)) {
		focal.taxa <- comparisons$taxa.0[sister_index][[1]]
		target.taxa <- comparisons$taxa.1[sister_index][[1]]
		local_df <- data.frame(node=comparisons$node[sister_index], depth=depths[phy_pruned$edge==comparisons$node[sister_index]][1], ntax.focal=comparisons$ntax.0[sister_index], ntax.target=comparisons$ntax.1[sister_index])
		if(local_df$depth<=maxdepth) {
			pairwise_result <- t(compare_two_clades(focal.taxa, target.taxa, tiprates=tiprates, phy_pruned=phy_pruned, rates=rates, ndraws=ndraws))
			if(!is.na(pairwise_result[1])) {
				local_df <- cbind(local_df,pairwise_result)
				if(sister_index==1) {
					global_df <-local_df
				} else {
					global_df <- rbind(global_df, local_df)
				}
			}
		}
	}
	return(global_df)
}


compare_cells_coarse <- function(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
	fine_results <- compare_cells_fine(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates, ndraws, maxdepth)
	coarse_results <- apply(subset(fine_results, select=c(paste0(rates, "_mean_target_minus_focal"))), 2, mean)
	coarse_results["npairs"] <- nrow(fine_results)
	return(coarse_results)
}

# Get mean and sd of difference across all pairs in a sister comparison of target (state 1, right) minus focal (state 0, left)
#' @param ndraws how many random pairs of taxa to sample
compare_two_clades <- function(focal.taxa, target.taxa, tiprates, phy_pruned, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100) {
	tiprates_focal <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[focal.taxa],]
	tiprates_target <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[target.taxa],]
	if(nrow(tiprates_focal)<1 | nrow(tiprates_target)<1) {
		means <- as.data.frame(t(rep(NA, length(rates))))
		colnames(means) <- c(paste0(rates, "_mean_target_minus_focal"))
		sds <- as.data.frame(t(rep(NA, length(rates))))
		return(cbind(means, sds))
	}
	local_results <- data.frame(matrix(NA, nrow=ndraws, ncol=length(rates)))
	for (sample_index in sequence(ndraws)) {
		focal_sample <- sample.int(nrow(tiprates_focal),1)
		target_sample <- sample.int(nrow(tiprates_target),1)
		for(rate_index in seq_along(rates)) {
			local_results[sample_index,rate_index] <- tiprates_target[target_sample, rates[rate_index]] - tiprates_focal[focal_sample, rates[rate_index]]
		}
	}
	summary_results <- c(apply(local_results, 2, mean), apply(local_results, 2, sd))
	names(summary_results) <- c(paste0(rates, "_mean_target_minus_focal"), paste0(rates, "_sd_target_minus_focal"))
	return(summary_results)
}

#' compare_rook(cell_bounds = c(minlat_focal=-10, maxlat_focal=10, minlon_focal=-100, maxlon_focal=-50), phy=phy_clean, latlon=latlon, tiprates=tiprates)
compare_rook <- function(cell_bounds, phy, latlon, tiprates, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
	minlat_focal <- unname(cell_bounds["minlat_focal"])
	maxlat_focal <- unname(cell_bounds["maxlat_focal"])
	minlon_focal <- unname(cell_bounds["minlon_focal"])
	maxlon_focal <- unname(cell_bounds["maxlon_focal"])

	lat_mid <- mean(c(maxlat_focal,minlat_focal))
	lon_mid <- mean(c(maxlon_focal, minlon_focal))
	lat_step <- maxlat_focal-minlat_focal
	lon_step <- maxlon_focal-minlon_focal
	north <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal+lat_step, maxlat_target=maxlat_focal+lat_step, minlon_target=minlon_focal, maxlon_target=maxlon_focal, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
	south <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal-lat_step, maxlat_target=maxlat_focal-lat_step, minlon_target=minlon_focal, maxlon_target=maxlon_focal, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
	east <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal, maxlat_target=maxlat_focal, minlon_target=minlon_focal+lon_step, maxlon_target=maxlon_focal+lon_step, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
	west <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal, maxlat_target=maxlat_focal, minlon_target=minlon_focal-lon_step, maxlon_target=maxlon_focal-lon_step, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
	direction_east_minus_west <- east-west
	direction_north_minus_south <- north-south
	tidy_df <- as.data.frame(rbind(north, south, east, west, direction_east_minus_west, direction_north_minus_south))
	tidy_df$direction <- c("north", "south", "east", "west", "east_minus_west", "north_minus_south")
	tidy_df$lat_mid=lat_mid
	tidy_df$lon_mid=lon_mid
	tidy_df$lat_step=lat_step
	tidy_df$lon_step=lon_step
	return(tidy_df)
	#return(list(tidy_df=tidy_df, lat_mid=lat_mid, lon_mid=lon_mid, lat_step=lat_step, lon_step=lon_step, north=north, south=south, east=east, west=west, direction_east_minus_west=direction_east_minus_west, direction_north_minus_south=direction_north_minus_south))
}

compute_sampling_grid <- function(latbins=25, lonbins=25, minlat=-70, maxlat=70, minlon=-170, maxlon=170) {
	lat_breaks <- seq(from=minlat, to=maxlat, length.out=1+latbins)
	lon_breaks <- seq(from=minlon, to=maxlon, length.out=1+lonbins)
	indices_matrix <- expand.grid(lat_index=sequence(latbins), lon_index=sequence(lonbins))
	focal_bounds <- data.frame(
		minlat_focal=lat_breaks[indices_matrix$lat_index], 
		maxlat_focal=lat_breaks[1+indices_matrix$lat_index],
		minlon_focal=lon_breaks[indices_matrix$lon_index], 
		maxlon_focal=lon_breaks[1+indices_matrix$lon_index]
	)
	return(focal_bounds)
}

#' results <- compare_all_grid(phy=phy_clean, latlon=latlon, tiprates=tiprates, sample_grid=compute_sampling_grid(), ncores=parallel::detectCores())
compare_all_grid <- function(phy, latlon, tiprates, sample_grid=compute_sampling_grid(), ncores=2, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
	results <- pbapply::pbapply(sample_grid, MARGIN=1, FUN=compare_rook, phy=phy, latlon=latlon, tiprates=tiprates, rates=rates, ndraws=ndraws, maxdepth=Inf, mincomparisons=mincomparisons, cl=ncores)
	return(results)
}

#from compare_all_grid
synthesize_results <- function(results, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
	for (i in seq_along(results)) {
		results[[i]]$cell_number <- i
	}
	is.bad <- function(x) {
		return(any(is.na(x[,1])))
	}
	# 
	results_all<- dplyr::bind_rows(results)
	rownames(results_all) <- NULL
	for (i in seq_along(rates)) {
		results_all[,paste0(rates[i], "_normalized")] <- results_all[,paste0(rates[i], "_mean_target_minus_focal")]/max(abs(results_all[,paste0(rates[i], "_mean_target_minus_focal")]))
	}
	
	results_good <- subset(results_all, !is.na(results_all[,1]))

	return(list(good=results_good, all=results_all))
}

plot_results_good <- function(results_good, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
	scale <- 1.4*min(c(results_good$lat_step, results_good$lon_step))
	eastwest <- subset(results_good,direction=="east_minus_west")
	northsouth <- subset(results_good,direction=="north_minus_south")
	valid_cells <- unique(eastwest$cell_number)[unique(eastwest$cell_number) %in% unique(northsouth$cell_number)]
	for(rate_index in seq_along(rates)) {
		plot(results_good$lon_mid, results_good$lat_mid,  asp=1, main=rates[rate_index], type="n")

		for (cell_index in seq_along(valid_cells)) {
			cell_df <- subset(results_good, cell_number==valid_cells[cell_index])
			eastwest_index <- which(cell_df$direction=="east_minus_west")
			northsouth_index <- which(cell_df$direction=="north_minus_south")

			#segments(x0=cell_df$lon_mid[1], y0=cell_df$lat_mid[1], x1=cell_df$lon_mid[1]+scale*cell_df[eastwest_index,paste0(rates[rate_index], "_normalized")], y1=cell_df$lat_mid[1]+scale*cell_df[northsouth_index,paste0(rates[rate_index], "_normalized")])
			points(eastwest$lon_mid, eastwest$lat_mid, pch=21, col="red")

		}
	}
}

krig_results <- function(results, focal_rate="turnover", npoints_latitude=500, npoints_longitude=1000) {
	results_zero <- results$all
	results_zero[is.na(results_zero)] <- 0
	results_zero <- results_zero[order(results_zero$lat_mid, results_zero$lon_mid),]
	results_zero_eastwest <- subset(results_zero, direction=="east_minus_west")
	results_zero_northsouth <- subset(results_zero, direction=="north_minus_south")
	longitudes <- sort(unique(results_zero$lon_mid))
	latitudes <- sort(unique(results_zero$lat_mid))
	latitudes_fine <- seq(from=min(longitudes), to=max(longitudes), length.out=npoints_latitude)
	longitudes_fine <- seq(from=min(latitudes), to=max(latitudes), length.out=npoints_longitude)
	z_eastwest <- matrix(nrow=length(latitudes), ncol=length(longitudes),0)
	z_northsouth <- z_eastwest 
	for (lat_index in seq_along(latitudes)) {
		for (lon_index in seq_along(longitudes)) {
			z_eastwest[lon_index, lat_index] <- subset(results_zero_eastwest, lat_mid==latitudes[lat_index] & lon_mid==longitudes[lon_index])[,paste0(focal_rate, "_mean_target_minus_focal")]
			z_northsouth[lon_index, lat_index] <- subset(results_zero_northsouth, lat_mid==latitudes[lat_index] & lon_mid==longitudes[lon_index])[,paste0(focal_rate, "_mean_target_minus_focal")]

		}
	}
	xpyp_grid <- expand.grid(x=longitudes_fine, y=latitudes_fine)
	krig_eastwest <- pracma::interp2(x=longitudes, y=latitudes, Z=z_eastwest, xp=xpyp_grid$x, yp=xpyp_grid$y)
	krig_eastwest_grid <- matrix(krig_eastwest, nrow=npoints_latitude, byrow=TRUE)

	krig_northsouth <- pracma::interp2(x=longitudes, y=latitudes, Z=z_northsouth, xp=xpyp_grid$x, yp=xpyp_grid$y)
	krig_northsouth_grid <- matrix(krig_eastwest, nrow=npoints_latitude, byrow=TRUE)
	return(list(latitudes=latitudes_fine, longitudes=longitudes_fine, eastwest=krig_eastwest_grid, northsouth=krig_northsouth_grid))
}

# Uses ideas from Vincent Zoonekynd, https://stackoverflow.com/a/14939043
dynamic_plot <- function(results, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
	results_zero <- results$all
	results_zero[is.na(results_zero)] <- 0
	results_zero_eastwest <- subset(results_zero, direction=="east_minus_west")
	results_zero_northsouth <- subset(results_zero, direction=="north_minus_south")
	longitudes <- sort(unique(results_zero$lon_mid))
	latitudes <- sort(unique(results_zero$lat_mid))
	y <- latitudes
	ymin <- min(y)
	ymax <- max(y)

	x <- longitudes
	xmin <- min(x)
	xmax <- max(x)

	i_to_x <- function(i) xmin + i / width  * (xmax - xmin)
	j_to_y <- function(j) ymin + j / height * (ymax - ymin)
	x_to_i <- function(x) pmin( width,  pmax( 1, floor( (x-xmin)/(xmax-xmin) * width  ) ) )
	y_to_j <- function(y) pmin( height, pmax( 1, floor( (y-ymin)/(ymax-ymin) * height ) ) )
	i <- col(z)
	j <- row(z)
	x <- i_to_x(i)
	y <- j_to_y(j)


	z_original <- matrix(runif(length(longitudes)*length(latitudes)),nr=length(latitudes))
	z <- z_original
	for (rate_index in seq_along(rates)) {
		z <- z_original
		res <- z_original
		for (k in 1:50) {

		}
	}



}