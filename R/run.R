load("inst/data/magnoliidae.rda")
source("R/spatialcontrast.R")
nclusters_vector <- c(2,3,4,5,10,15,20,25,30,35,40)

#nclusters_vector <- c(2,3,4)
nclusters_results <- list()
for (cluster_index in seq_along(nclusters_vector)) {
    latlon_new <- cluster_latlon(latlon, nclusters=nclusters_vector[cluster_index])
    results <- compare_clusters(phy=phy_clean, latlon_new, tiprates)
    summaries <- summarize_cluster_results(results)
    nclusters_results[[cluster_index]] <- list(latlon=latlon_new, results=results, summaries=summaries)
    pdf(file=paste0("inst/data/cluster_", nclusters_vector[cluster_index], ".pdf"))
    plot_map(latlon=nclusters_results[[cluster_index]]$latlon, summaries=nclusters_results[[cluster_index]]$summaries, "net.div")
    dev.off()
}
save(nclusters_results, file="inst/data/nclusters_results.rda")
