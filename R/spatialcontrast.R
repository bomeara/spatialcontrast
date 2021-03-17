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

#' @examples 
#' data(magnoliiadae)
#' compare_cells(phy_clean, latlon, minlat_focal=-10, maxlat_focal=10, minlon_focal=-100, maxlon_focal=-50, minlat_target=10, maxlat_target=20, minlon_target=-100, maxlon_target=-50)
compare_cells <- function(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target) {
	focal_presence <- present_in_cell(phy, subset_cell(latlon, minlat=minlat_focal, maxlat=maxlat_focal, minlon=minlon_focal, maxlon=maxlon_focal))
	target_presence <- present_in_cell(phy, subset_cell(latlon, minlat=minlat_target, maxlat=maxlat_target, minlon=minlon_target, maxlon=maxlon_target))
	total_presence <- focal_presence + target_presence
	taxa_to_keep <- names(total_presence)[which(total_presence==1)]
	if(length(taxa_to_keep)<2) {
		warning("Not enough taxa to compare these cells")
		return(NA)
	}
	phy_pruned <- ape::keep.tip(phy, taxa_to_keep)
	target_pruned <- target_presence[taxa_to_keep]
	comparisons <- sisters::sis_format_comparison(sisters=sisters::sis_get_sisters(phy_pruned, ncores=1), trait=target_pruned, phy=phy_pruned)
	heights <- phytools::nodeHeights(phy_pruned)
	depths <- max(heights) - heights
	for (sister_index in seq_along(comparisons$node)) {
		left <- comparisons$left[sister_index][[1]]
		right <- comparisons$right[sister_index][[1]]
	}

	return(comparisons)
}

compare_two_clades <- function(left, right, tiprates, phy_pruned, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
	tiprates_left <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[left],]
	tiprates_right <- tiprates[tiprates$taxon %in% phy_pruned$tip.label[right],]

}