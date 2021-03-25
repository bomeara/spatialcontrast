

present_in_cluster <- function(phy, latlon, cluster_id) {
	trait <- rep(0, ape::Ntip(phy))
	names(trait) <- phy$tip.label
	latlon_sub <- subset(latlon, cluster==cluster_id)
	trait[names(trait) %in% latlon_sub$taxon] <- 1
	return(trait)
}


#' Cluster points based on latitude, longitude, and potentially other traits
#' 
#' To compare between regions, we must first have assigned points to regions. This function uses k-means to do this. By default you can use latitude and longitude, but you can optionally use other factors instead or in addition to these factors. For example, if you want the algorithm to attempt to not split species into separate clusters, you could make a column that converts the taxon to a number and use that number. If nclusters is specified, the function will split into that many clusters; if not, it will use whatever mean_points_per_cluster is set to in order to divide into clusters of approximately that size. If data are missing for any of the split_traits columns, that row will be deleted.
#' 
#' @param latlon Data.frame with fields for taxon, lon, lat and perhaps other fields
#' @param mean_points_per_cluster How many points *on average* to assign to each cluster -- some will have more or fewer
#' @param nclusters If specified, how many clusters to divide points into.
#' @param split_traits The columns to use for clustering. 
#' @param do_scale If TRUE, normalizes each column before running
#' @return The latlon data.frame with an additional column for cluster id and cluster_* columns having the raw mean of all points in that cluster for that trait (i.e., cluster_lat for the mean latitude). 
#' @export
cluster_latlon <- function(latlon, mean_points_per_cluster=10000, nclusters=NULL, split_traits = c("lat", "lon"), do_scale=TRUE) {
	for (trait_index in seq_along(split_traits)) {
		latlon <- latlon[!is.na(latlon[,split_traits[trait_index]]),]	

	}
	latlon_focal <- latlon[,split_traits]
	if(do_scale) {
		latlon_focal[,split_traits[trait_index]] <- scale(latlon_focal[,split_traits[trait_index]])
	}
	if(is.null(nclusters)){
		 nclusters<-floor(nrow(latlon)/mean_points_per_cluster)
	}
	clusters <- stats::kmeans(latlon_focal, centers=nclusters, nstart=20, iter.max=100)
	latlon$cluster <- clusters$cluster
	for(trait_index in seq_along(split_traits)) {
		latlon[,paste0("cluster_", split_traits[trait_index])] <- NA
		for(cluster_index in sequence(max(latlon$cluster))) {
			relevant_rows <- which(latlon$cluster==cluster_index)
			latlon[relevant_rows,paste0("cluster_", split_traits[trait_index])] <- mean(latlon[relevant_rows, split_traits[trait_index]])
		}
	}
	return(latlon)
}

#' Compare spatial clusters
#' 
#' This is the main workhorse function: it compares pairs of clusters by subsampling a tree to only taxa present in at least one of those areas. Then it uses the sisters package to find all bifurcations where all taxa in one branch are in one area and all taxa in the other branch are in the other area (if one taxon is in both areas, it doesn't use this comparison). It takes the rate estimate for a taxon in the "target" cluster branch and subtracts the rate estimate for a taxon in the "focal" cluster branch (the columns and rows of the output matrix, respectively). If at least one of the clades being compared has more than one taxon, taxa are repeatedly sampled at random.
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
					print(c(focal, target, max(latlon$cluster)))
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


plot_map <- function(latlon, summaries, focal_rate="net.div", focal_cluster=1) {
	library(ggplot2)
	mp <- NULL
	mapWorld <- ggplot2::borders("world", colour="gray50", fill="gray50") # create a layer of borders
	mp <- ggplot2::ggplot() + mapWorld
	values <- summaries$t_test_estimate[focal_rate][[1]][focal_cluster,]
	latlon$values_to_plot <- values[latlon$cluster]
	mp <- mp + geom_point(data=latlon, mapping=aes(x=lon, y=lat, col=values_to_plot)) + scale_colour_gradient2(low="blue", high="red", na.value="black", name=paste0("Diff in ",focal_rate)) + ggtitle(paste0(focal_rate, " difference for cluster ", focal_cluster))
	local_latlon <- latlon[!duplicated(latlon$cluster),]
	mp <- mp + geom_label(aes(x=cluster_lon, y=cluster_lat, label=cluster), data=local_latlon)
	print(mp)
}


plot_map_signif <- function(latlon, summaries, focal_rate="net.div") {
	library(ggplot2)
	cluster_info <- apply(summaries$t_test_estimate[focal_rate][[1]], 2, stats::t.test)
	limits <- as.data.frame(lapply(cluster_info, "[[", "conf.int"))
	centers <- simplify2array(lapply(cluster_info, "[[", "estimate"))
	mp <- NULL
	mapWorld <- ggplot2::borders("world", colour="gray50", fill="gray50") # create a layer of borders
	mp <- ggplot2::ggplot() + mapWorld
	#values <- summaries$t_test_estimate[focal_rate][[1]][focal_cluster,]
	latlon$values_to_plot <- unlist(unname(limits[1,latlon$cluster]))
	mp <- mp + geom_point(data=latlon, mapping=aes(x=lon, y=lat, col=values_to_plot, size=1.5)) + scale_colour_gradient2(low="blue", high="red", limits=range(limits), na.value="black", name=paste0("Diff in ",focal_rate)) 
	
	latlon$values_to_plot <- unlist(unname(limits[2,latlon$cluster]))
	mp <- mp + geom_point(data=latlon, mapping=aes(x=lon, y=lat, col=values_to_plot, size=1.1)) + scale_colour_gradient2(low="blue", high="red", limits=range(limits), na.value="black", name=paste0("Diff in ",focal_rate)) 

	latlon$values_to_plot <- unlist(unname(centers[latlon$cluster]))
	mp <- mp + geom_point(data=latlon, mapping=aes(x=lon, y=lat, col=values_to_plot, size=1)) + scale_colour_gradient2(low="blue", high="red", limits=range(limits), na.value="black", name=paste0("Diff in ",focal_rate)) 


	local_latlon <- latlon[!duplicated(latlon$cluster),]
	mp <- mp + geom_label(aes(x=cluster_lon, y=cluster_lat, label=cluster), data=local_latlon)
	print(mp)
}
