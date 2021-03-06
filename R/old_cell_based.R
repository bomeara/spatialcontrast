# #' @param latlon a data.frame with taxon, taxon_full, lon, lat
# #' @param minlat all points will be greater than or equal to this latitude
# #' @param maxlat all points will be less than this latitude
# #' @param minlon all points will be greater than or equal to this longitude
# #' @param maxlon all points will be less than this longitude
# subset_cell <- function(latlon, minlat=-90, maxlat=90, minlon=-180, maxlon=180) {
# 	return(subset(latlon, latlon$lat>=minlat & latlon$lat<maxlat & latlon$lon>=minlon & latlon$lon < maxlon))
# }



# present_in_cell <- function(phy, subsetted_cell) {
# 	trait <- rep(0, ape::Ntip(phy))
# 	names(trait) <- phy$tip.label
# 	trait[names(trait) %in% subsetted_cell$taxon] <- 1
# 	return(trait)
# }




# #' @param maxdepth How far back in time to make comparisons
# #' @examples 
# #' data(magnoliiadae)
# #' compare_cells_fine(phy_clean, latlon, tiprates, minlat_focal=-10, maxlat_focal=10, minlon_focal=-100, maxlon_focal=-50, minlat_target=10, maxlat_target=20, minlon_target=-100, maxlon_target=-50)
# compare_cells_fine <- function(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
# 	focal_presence <- present_in_cell(phy, subset_cell(latlon, minlat=minlat_focal, maxlat=maxlat_focal, minlon=minlon_focal, maxlon=maxlon_focal))
# 	target_presence <- present_in_cell(phy, subset_cell(latlon, minlat=minlat_target, maxlat=maxlat_target, minlon=minlon_target, maxlon=maxlon_target))
# 	total_presence <- focal_presence + target_presence
# 	taxa_to_keep <- names(total_presence)[which(total_presence==1)]
# 	if(length(taxa_to_keep)<2) {
# 		warning("Not enough taxa to compare these cells")
# 		means <- as.data.frame(t(rep(NA, length(rates))))
# 		colnames(means) <- c(paste0(rates, "_mean_target_minus_focal"))
# 		sds <- as.data.frame(t(rep(NA, length(rates))))
# 		colnames(sds) <- c(paste0(rates, "_mean_target_minus_focal"))		
# 		return(cbind(data.frame(node=NA, depth=NA, ntax.left=0, ntax.right=0), means, sds))
# 	}
# 	global_df <- data.frame()
# 	phy_pruned <- ape::keep.tip(phy, taxa_to_keep)
# 	target_pruned <- target_presence[taxa_to_keep]
# 	sink("/dev/null")
# 	comparisons <- sanitize_comparisons(sisters::sis_format_comparison(sisters=sisters::sis_get_sisters(phy_pruned, ncores=1), trait=target_pruned, phy=phy_pruned))
# 	sink()
# 	if(length(comparisons$node)<mincomparisons) {
# 		warning("Not enough sister pairs to compare these cells")
# 		means <- as.data.frame(t(rep(NA, length(rates))))
# 		colnames(means) <- c(paste0(rates, "_mean_target_minus_focal"))
# 		sds <- as.data.frame(t(rep(NA, length(rates))))
# 		colnames(sds) <- c(paste0(rates, "_mean_target_minus_focal"))		
# 		return(cbind(data.frame(node=NA, depth=NA, ntax.left=0, ntax.right=0), means, sds))
# 	}
# 	heights <- phytools::nodeHeights(phy_pruned)
# 	depths <- max(heights) - heights
# 	#cat("\n")
# 	#pb = utils::txtProgressBar(min = 0, max = length(comparisons$node), initial = 0) 
# 	for (sister_index in seq_along(comparisons$node)) {
# 		focal.taxa <- comparisons$taxa.0[sister_index][[1]]
# 		target.taxa <- comparisons$taxa.1[sister_index][[1]]
# 		local_df <- data.frame(node=comparisons$node[sister_index], depth=depths[phy_pruned$edge==comparisons$node[sister_index]][1], ntax.focal=comparisons$ntax.0[sister_index], ntax.target=comparisons$ntax.1[sister_index])
# 		if(local_df$depth<=maxdepth) {
# 			pairwise_result <- t(compare_two_clades(focal.taxa, target.taxa, tiprates=tiprates, phy_pruned=phy_pruned, rates=rates, ndraws=ndraws))
# 			if(!is.na(pairwise_result[1])) {
# 				local_df <- cbind(local_df,pairwise_result)
# 				if(sister_index==1) {
# 					global_df <-local_df
# 				} else {
# 					global_df <- rbind(global_df, local_df)
# 				}
# 			}
# 		}
# 	}
# 	return(global_df)
# }


# compare_cells_coarse <- function(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
# 	fine_results <- compare_cells_fine(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target, maxlat_target, minlon_target, maxlon_target, rates, ndraws, maxdepth)
# 	coarse_results <- apply(subset(fine_results, select=c(paste0(rates, "_mean_target_minus_focal"))), 2, mean)
# 	coarse_results["npairs"] <- nrow(fine_results)
# 	return(coarse_results)
# }


# #' compare_rook(cell_bounds = c(minlat_focal=-10, maxlat_focal=10, minlon_focal=-100, maxlon_focal=-50), phy=phy_clean, latlon=latlon, tiprates=tiprates)
# compare_rook <- function(cell_bounds, phy, latlon, tiprates, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
# 	minlat_focal <- unname(cell_bounds["minlat_focal"])
# 	maxlat_focal <- unname(cell_bounds["maxlat_focal"])
# 	minlon_focal <- unname(cell_bounds["minlon_focal"])
# 	maxlon_focal <- unname(cell_bounds["maxlon_focal"])

# 	lat_mid <- mean(c(maxlat_focal,minlat_focal))
# 	lon_mid <- mean(c(maxlon_focal, minlon_focal))
# 	lat_step <- maxlat_focal-minlat_focal
# 	lon_step <- maxlon_focal-minlon_focal
# 	north <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal+lat_step, maxlat_target=maxlat_focal+lat_step, minlon_target=minlon_focal, maxlon_target=maxlon_focal, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
# 	south <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal-lat_step, maxlat_target=maxlat_focal-lat_step, minlon_target=minlon_focal, maxlon_target=maxlon_focal, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
# 	east <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal, maxlat_target=maxlat_focal, minlon_target=minlon_focal+lon_step, maxlon_target=maxlon_focal+lon_step, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
# 	west <- compare_cells_coarse(phy, latlon, tiprates, minlat_focal, maxlat_focal, minlon_focal, maxlon_focal, minlat_target=minlat_focal, maxlat_target=maxlat_focal, minlon_target=minlon_focal-lon_step, maxlon_target=maxlon_focal-lon_step, rates=rates, ndraws=ndraws, maxdepth=maxdepth, mincomparisons=mincomparisons)
# 	direction_east_minus_west <- east-west
# 	direction_north_minus_south <- north-south
# 	tidy_df <- as.data.frame(rbind(north, south, east, west, direction_east_minus_west, direction_north_minus_south))
# 	tidy_df$direction <- c("north", "south", "east", "west", "east_minus_west", "north_minus_south")
# 	tidy_df$lat_mid=lat_mid
# 	tidy_df$lon_mid=lon_mid
# 	tidy_df$lat_step=lat_step
# 	tidy_df$lon_step=lon_step
# 	return(tidy_df)
# 	#return(list(tidy_df=tidy_df, lat_mid=lat_mid, lon_mid=lon_mid, lat_step=lat_step, lon_step=lon_step, north=north, south=south, east=east, west=west, direction_east_minus_west=direction_east_minus_west, direction_north_minus_south=direction_north_minus_south))
# }

# compute_sampling_grid <- function(latbins=25, lonbins=25, minlat=-70, maxlat=70, minlon=-170, maxlon=170) {
# 	lat_breaks <- seq(from=minlat, to=maxlat, length.out=1+latbins)
# 	lon_breaks <- seq(from=minlon, to=maxlon, length.out=1+lonbins)
# 	indices_matrix <- expand.grid(lat_index=sequence(latbins), lon_index=sequence(lonbins))
# 	focal_bounds <- data.frame(
# 		minlat_focal=lat_breaks[indices_matrix$lat_index], 
# 		maxlat_focal=lat_breaks[1+indices_matrix$lat_index],
# 		minlon_focal=lon_breaks[indices_matrix$lon_index], 
# 		maxlon_focal=lon_breaks[1+indices_matrix$lon_index]
# 	)
# 	return(focal_bounds)
# }

# #' results <- compare_all_grid(phy=phy_clean, latlon=latlon, tiprates=tiprates, sample_grid=compute_sampling_grid(), ncores=parallel::detectCores())
# compare_all_grid <- function(phy, latlon, tiprates, sample_grid=compute_sampling_grid(), ncores=2, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction"), ndraws=100, maxdepth=Inf, mincomparisons=3) {
# 	results <- pbapply::pbapply(sample_grid, MARGIN=1, FUN=compare_rook, phy=phy, latlon=latlon, tiprates=tiprates, rates=rates, ndraws=ndraws, maxdepth=Inf, mincomparisons=mincomparisons, cl=ncores)
# 	return(results)
# }

# #from compare_all_grid
# synthesize_results <- function(results, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
# 	for (i in seq_along(results)) {
# 		results[[i]]$cell_number <- i
# 	}
# 	is.bad <- function(x) {
# 		return(any(is.na(x[,1])))
# 	}
# 	# 
# 	results_all<- dplyr::bind_rows(results)
# 	rownames(results_all) <- NULL
# 	for (i in seq_along(rates)) {
# 		results_all[,paste0(rates[i], "_normalized")] <- results_all[,paste0(rates[i], "_mean_target_minus_focal")]/max(abs(results_all[,paste0(rates[i], "_mean_target_minus_focal")]))
# 	}
	
# 	results_good <- subset(results_all, !is.na(results_all[,1]))

# 	return(list(good=results_good, all=results_all))
# }

# plot_results_good <- function(results_good, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
# 	scale <- 1.4*min(c(results_good$lat_step, results_good$lon_step))
# 	eastwest <- subset(results_good,direction=="east_minus_west")
# 	northsouth <- subset(results_good,direction=="north_minus_south")
# 	valid_cells <- unique(eastwest$cell_number)[unique(eastwest$cell_number) %in% unique(northsouth$cell_number)]
# 	for(rate_index in seq_along(rates)) {
# 		plot(results_good$lon_mid, results_good$lat_mid,  asp=1, main=rates[rate_index], type="n")

# 		for (cell_index in seq_along(valid_cells)) {
# 			cell_df <- subset(results_good, cell_number==valid_cells[cell_index])
# 			eastwest_index <- which(cell_df$direction=="east_minus_west")
# 			northsouth_index <- which(cell_df$direction=="north_minus_south")

# 			#segments(x0=cell_df$lon_mid[1], y0=cell_df$lat_mid[1], x1=cell_df$lon_mid[1]+scale*cell_df[eastwest_index,paste0(rates[rate_index], "_normalized")], y1=cell_df$lat_mid[1]+scale*cell_df[northsouth_index,paste0(rates[rate_index], "_normalized")])
# 			points(eastwest$lon_mid, eastwest$lat_mid, pch=21, col="red")

# 		}
# 	}
# }

# krig_results <- function(results, focal_rate="turnover", npoints_latitude=500, npoints_longitude=1000) {
# 	results_zero <- results$all
# 	results_zero[is.na(results_zero)] <- 0
# 	results_zero <- results_zero[order(results_zero$lat_mid, results_zero$lon_mid),]
# 	results_zero_eastwest <- subset(results_zero, direction=="east_minus_west")
# 	results_zero_northsouth <- subset(results_zero, direction=="north_minus_south")
# 	longitudes <- sort(unique(results_zero$lon_mid))
# 	latitudes <- sort(unique(results_zero$lat_mid))
# 	latitudes_fine <- seq(from=min(longitudes), to=max(longitudes), length.out=npoints_latitude)
# 	longitudes_fine <- seq(from=min(latitudes), to=max(latitudes), length.out=npoints_longitude)
# 	z_eastwest <- matrix(nrow=length(latitudes), ncol=length(longitudes),0)
# 	z_northsouth <- z_eastwest 
# 	for (lat_index in seq_along(latitudes)) {
# 		for (lon_index in seq_along(longitudes)) {
# 			z_eastwest[lon_index, lat_index] <- subset(results_zero_eastwest, lat_mid==latitudes[lat_index] & lon_mid==longitudes[lon_index])[,paste0(focal_rate, "_mean_target_minus_focal")]
# 			z_northsouth[lon_index, lat_index] <- subset(results_zero_northsouth, lat_mid==latitudes[lat_index] & lon_mid==longitudes[lon_index])[,paste0(focal_rate, "_mean_target_minus_focal")]

# 		}
# 	}
# 	xpyp_grid <- expand.grid(x=longitudes_fine, y=latitudes_fine)
# 	krig_eastwest <- pracma::interp2(x=longitudes, y=latitudes, Z=z_eastwest, xp=xpyp_grid$x, yp=xpyp_grid$y)
# 	krig_eastwest_grid <- matrix(krig_eastwest, nrow=npoints_latitude, byrow=TRUE)

# 	krig_northsouth <- pracma::interp2(x=longitudes, y=latitudes, Z=z_northsouth, xp=xpyp_grid$x, yp=xpyp_grid$y)
# 	krig_northsouth_grid <- matrix(krig_eastwest, nrow=npoints_latitude, byrow=TRUE)
# 	return(list(latitudes=latitudes_fine, longitudes=longitudes_fine, eastwest=krig_eastwest_grid, northsouth=krig_northsouth_grid))
# }

# # Uses ideas from Vincent Zoonekynd, https://stackoverflow.com/a/14939043
# dynamic_plot <- function(results, rates=c("turnover", "net.div", "speciation", "extinct.frac", "extinction")) {
# 	results_zero <- results$all
# 	results_zero[is.na(results_zero)] <- 0
# 	results_zero_eastwest <- subset(results_zero, direction=="east_minus_west")
# 	results_zero_northsouth <- subset(results_zero, direction=="north_minus_south")
# 	longitudes <- sort(unique(results_zero$lon_mid))
# 	latitudes <- sort(unique(results_zero$lat_mid))
# 	y <- latitudes
# 	ymin <- min(y)
# 	ymax <- max(y)

# 	x <- longitudes
# 	xmin <- min(x)
# 	xmax <- max(x)

# 	i_to_x <- function(i) xmin + i / width  * (xmax - xmin)
# 	j_to_y <- function(j) ymin + j / height * (ymax - ymin)
# 	x_to_i <- function(x) pmin( width,  pmax( 1, floor( (x-xmin)/(xmax-xmin) * width  ) ) )
# 	y_to_j <- function(y) pmin( height, pmax( 1, floor( (y-ymin)/(ymax-ymin) * height ) ) )
# 	i <- col(z)
# 	j <- row(z)
# 	x <- i_to_x(i)
# 	y <- j_to_y(j)


# 	z_original <- matrix(runif(length(longitudes)*length(latitudes)),nr=length(latitudes))
# 	z <- z_original
# 	for (rate_index in seq_along(rates)) {
# 		z <- z_original
# 		res <- z_original
# 		for (k in 1:50) {

# 		}
# 	}
# }