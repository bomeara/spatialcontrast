### setwd("~/Desktop/spatialcontrast")
# load("inst/data/magnoliidae.rda")
# source("R/spatialcontrast.R")

#--------------
# nclusters_vector <- c(2,3,4,5,10,15,20,25,30,35,40)
# nclusters_vector <- c(3,4,5) # testing with fewer clusters for speed
# nclusters_results <- list()
# for (cluster_index in seq_along(nclusters_vector)) {
#     latlon_new <- cluster_latlon(latlon, nclusters=nclusters_vector[cluster_index], do_scale=F, split_traits="temp") # testing with temp
#     results <- compare_clusters(phy=phy_clean, latlon_new, tiprates)
#     summaries <- summarize_cluster_results(results)
#     nclusters_results[[cluster_index]] <- list(latlon=latlon_new, results=results, summaries=summaries)
#     pdf(file=paste0("inst/data/cluster_", nclusters_vector[cluster_index], ".pdf"))
#    for(individual_cluster in sequence(nclusters_vector[cluster_index])) {
#         plot_map(latlon=nclusters_results[[cluster_index]]$latlon, summaries=nclusters_results[[cluster_index]]$summaries, focal_rate="net.div", focal_cluster=individual_cluster)
#     }
#    #plot_map_signif(latlon=nclusters_results[[cluster_index]]$latlon, summaries=nclusters_results[[cluster_index]]$summaries, focal_rate="turnover") 
#    dev.off()
# }
# save(nclusters_results, file="inst/data/nclusters_results.rda")

#--------------
# Test with pre-defined bioregions
 latlon_new <-  cluster_bioregion(bioregion_from_points(latlon))
 results <- compare_clusters(phy=phy_clean, latlon_new, tiprates)
 summaries <- summarize_cluster_results(results)

 pdf(file=paste0("inst/data/cluster_test_biomes.pdf"))
 for(individual_cluster in sequence(length(unique(latlon_new$cluster)))) {
     plot_map(latlon=latlon_new, summaries=summaries, focal_rate="turnover", focal_cluster=individual_cluster)
 }
 dev.off()

 pdf(file=paste0("inst/data/cluster_test_biomes_signif.pdf"))
 plot_map_signif(latlon=latlon_new, summaries=summaries, focal_rate="turnover") 
 dev.off()

 
 