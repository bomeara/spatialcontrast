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
	comparisons <- sisters::sis_format_comparison(sisters=sisters::sis_get_sisters(phy_pruned, ncores=1), trait=target_pruned, phy=phy_pruned)
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
		left <- comparisons$left[sister_index][[1]]
		right <- comparisons$right[sister_index][[1]]
		local_df <- data.frame(node=comparisons$node[sister_index], depth=depths[phy_pruned$edge==comparisons$node[sister_index]][1], ntax.left=comparisons$ntax.left[sister_index], ntax.right=comparisons$ntax.right[sister_index])
		if(local_df$depth<=maxdepth) {
			pairwise_result <- t(compare_two_clades(left, right, tiprates=tiprates, phy_pruned=phy_pruned, rates=rates, ndraws=ndraws))
			if(!is.na(pairwise_result[1])) {
				local_df <- cbind(local_df,pairwise_result)
				if(sister_index==1) {
					global_df <-local_df
				} else {
					global_df <- rbind(global_df, local_df)
				}
			}
		}
		#setTxtProgressBar(pb,sister_index)
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
compare_two_clades <- function(left, right, tiprates, phy_pruned, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100) {
	tiprates_left <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[left],]
	tiprates_right <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[right],]
	if(nrow(tiprates_right)<1 | nrow(tiprates_left)<1) {
		means <- as.data.frame(t(rep(NA, length(rates))))
		colnames(means) <- c(paste0(rates, "_mean_target_minus_focal"))
		sds <- as.data.frame(t(rep(NA, length(rates))))
		return(cbind(means, sds))
	}
	local_results <- data.frame(matrix(NA, nrow=ndraws, ncol=length(rates)))
	for (sample_index in sequence(ndraws)) {
		left_sample <- sample.int(nrow(tiprates_left),1)
		right_sample <- sample.int(nrow(tiprates_right),1)
		for(rate_index in seq_along(rates)) {
			local_results[sample_index,rate_index] <- tiprates_right[right_sample, rates[rate_index]] - tiprates_left[left_sample, rates[rate_index]]
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

compute_sampling_grid <- function(latbins=10, lonbins=10, minlat=-70, maxlat=70, minlon=-170, maxlon=170) {
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
synthesize_results <- function(results) {

}