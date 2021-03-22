library(ape)
gbif_dir <- "/Users/bomeara/Documents/MyDocuments/GitClones/MiSSEGradient/full_run/Data/2_gbif"
setwd(gbif_dir)
phy_taxized <- read.tree("../1_trees/magnoliidae.taxized.tre") 
phy_clean <- read.tree("../1_trees/magnoliidae.tre")

latlon <- data.frame()
for (i in seq_along(phy_clean$tip.label)) {
	local.csv <- NULL
	relevant.files <- list.files(pattern=paste0("magnoliidae__", phy_clean$tip.label[i], ".*"))
	for (j in seq_along(relevant.files)) {
		try(local.csv <- read.delim(relevant.files[j], header=FALSE))		
		if(!is.null(local.csv)) {
			local.csv$taxon <- phy_clean$tip.label[i]
			latlon <- rbind(latlon, local.csv)
		}
	}
	print(paste0(i, " of ", ape::Ntip(phy_clean)))
}
latlon <- data.frame(taxon=latlon$taxon, taxon_full=latlon$V1, lon=latlon$V4, lat=latlon$V3)
tiprates <- read.csv("../8_rates/magnoliidae_tip_rates.csv")
tiprates$taxon <- gsub(" ", "_", tiprates$taxon)
save(phy_taxized, phy_clean, latlon, tiprates, file="/Users/bomeara/Documents/MyDocuments/GitClones/spatialcontrast/inst/data/magnoliidae.rda")
