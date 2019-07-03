# Figure 3A Marker tSNE
# Matthew Buckley

rm(list=ls())
library(Seurat)
library(ggplot2)

# Change directory for testing purposes. Subsequent paths are relative.
setwd("~/Desktop/Dropbox/DBN2019/10x/tsne")

# Load Data
load(file = "../seurat/data/svz_All6_Filtered_2019-01-31.rda")
d <- svz_alldata
features <- colnames(d)

#===================================================================================================
# Ifng response genes
# Figure 3A

load(file = "../data/mouse.gene.hallmark.kegg.reactome.rda")
genelist <- tolower(hall.kegg.react.list$HALLMARK_INTERFERON_GAMMA_RESPONSE)

colnames(d) <- tolower(colnames(d))
ifn_data <- d[, colnames(d) %in% genelist] # dim: 14685   310
ifn_data$ifng_response <- rowSums(ifn_data)
meta <- d[, colnames(d) %in% c("age", "replicate", "celltype", "tsne_1", "tsne_2")]
d2 <- cbind(meta, ifn_data)
d2$age <- factor(d2$age, levels = c("y", "o"), ordered=TRUE)

low_col <- "grey"
high_col <- "darkred"
size_range <- c(.8, 1.2)
alpha_range <- c(0.10, .90)

# Saturate color at 99.5% percentile of response
threshold <- quantile(d2$ifng_response, .995)
ifng_sig_cap <- d2$ifng_response
ifng_sig_cap[ifng_sig_cap > threshold] <- threshold
d2$ifng_sig_cap <- ifng_sig_cap

# Plot IFN gamma response
q <- ggplot(data = d2, aes(x = tsne_1, y = tsne_2,
			size = ifng_sig_cap, alpha = ifng_sig_cap, color = ifng_sig_cap))
q <- q + geom_point() + scale_size(range = size_range) + scale_alpha(range = alpha_range)
q <- q + scale_colour_gradient(low = low_col, high = high_col)
q <- q + theme(legend.position="none")
ggsave(paste0("plots/Fig3a_IFNg_Response_", Sys.Date(), ".pdf"), q, height = 7, width = 6)

#===================================================================================================
# Ifng Receptors
# SupFig 4A

d <- svz_alldata
features <- colnames(d)
ifngr_genes <- sort(features[grep("Ifngr.", features)])
ifngr_sig <- rowSums(d[,ifngr_genes])
d$ifngr_sig <- ifngr_sig

low_col <- "grey"
high_col <- "darkred"
size_range <- c(.8, 1.2)
alpha_range <- c(0.15, .95)

# Saturate color at 99.5% percentile of log-normalized expression.
threshhold <- quantile(d$ifngr_sig, .995) # caps at 4.5
ifngr_sig_cap <- d$ifngr_sig
ifngr_sig_cap[ifngr_sig_cap > threshold] <- threshold
d$ifngr_sig_cap <- ifngr_sig_cap
d$Age <- factor(d$Age, levels=c("y","o"), ordered=T)


# Plot IFN gamma receptor
q <- ggplot(data = d, aes(x = tSNE_1, y = tSNE_2,
			size = ifngr_sig_cap, alpha = ifngr_sig_cap, color = ifngr_sig_cap))
q <- q + geom_point() + scale_size(range = size_range) + scale_alpha(range = alpha_range)
q <- q + scale_colour_gradient(low = low_col, high = high_col)
q <- q + theme(legend.position="none")
ggsave(paste0("plots/Fig4a_IFNg_Receptors_", Sys.Date(), ".pdf"), q, height = 7, width = 6)

#===================================================================================================
# Bst2
# SupFig 4I

# Saturate color at 99.5% percentile of log-normalized expression.
threshold <- quantile(d$Bst2, .995) # caps at 3.07
bst2_sig_cap <- d$Bst2
bst2_sig_cap[bst2_sig_cap > threshold] <- threshold
d$bst2_sig_cap <- bst2_sig_cap

low_col <- "grey"
high_col <- "darkred"
size_range <- c(.8, 1.2)
alpha_range <- c(0.15, .95)

# Plot IFN gamma receptor
q <- ggplot(data = d, aes(x = tSNE_1, y = tSNE_2,
			size = bst2_sig_cap, alpha = bst2_sig_cap, color = bst2_sig_cap))
q <- q + geom_point() + scale_size(range = size_range) + scale_alpha(range = alpha_range)
q <- q + scale_colour_gradient(low = low_col, high = high_col)
q <- q + theme(legend.position="none")
ggsave(paste0("plots/Fig4h_Bst2_", Sys.Date(), ".pdf"), q, height = 7, width = 6)

