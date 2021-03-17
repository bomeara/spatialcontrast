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
	comparisons <- sisters::sis_format_comparison(sisters=sisters::sis_get_sisters(phy_pruned, ncores=1), trait=target_pruned, phy=phy_pruned)
	if(length(comparisons$nodes)<mincomparisons) {
		warning("Not enough sister pairs to compare these cells")
		means <- as.data.frame(t(rep(NA, length(rates))))
		colnames(means) <- c(paste0(rates, "_mean_target_minus_focal"))
		sds <- as.data.frame(t(rep(NA, length(rates))))
		colnames(sds) <- c(paste0(rates, "_mean_target_minus_focal"))		
		return(cbind(data.frame(node=NA, depth=NA, ntax.left=0, ntax.right=0), means, sds))
	}
	heights <- phytools::nodeHeights(phy_pruned)
	depths <- max(heights) - heights
	for (sister_index in seq_along(comparisons$node)) {
		left <- comparisons$left[sister_index][[1]]
		right <- comparisons$right[sister_index][[1]]
		local_df <- data.frame(node=comparisons$node[sister_index], depth=depths[phy_pruned$edge==comparisons$node[sister_index]][1], ntax.left=comparisons$ntax.left[sister_index], ntax.right=comparisons$ntax.right[sister_index])
		if(local_df$depth<=maxdepth) {
			local_df <- cbind(local_df, t(compare_two_clades(left, right, tiprates=tiprates, phy_pruned=phy_pruned, rates=rates, ndraws=ndraws)))
			if(sister_index==1) {
				global_df <-local_df
			} else {
				global_df <- rbind(global_df, local_df)
			}
		}
	}

	return(global_df)
}

compare_cells_coarse <- function(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
	fine_results <- compare_cells_fine(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates, ndraws, maxdepth)
	coarse_results <- apply(subset(fine_results, select=c(paste0(rates, "_mean_target_minus_focal"))), 2, mean)
	coarse_results["npairs"] <- nrow(fine_results)
}

# Get mean and sd of difference across all pairs in a sister comparison of target (state 1, right) minus focal (state 0, left)
#' @param ndraws how many random pairs of taxa to sample
compare_two_clades <- function(left, right, tiprates, phy_pruned, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100) {
	tiprates_left <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[left],]
	tiprates_right <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[right],]
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

