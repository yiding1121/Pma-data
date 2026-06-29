if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install("syntenet")

library("syntenet")

fasta_dir="~/Desktop/Pma_data/orthofinder_output/Pm_data/"

grangeslist <- gff2GRangesList(fasta_dir)
aastringsetlist <- fasta2AAStringSetlist(fasta_dir)
check_input(seq = aastringsetlist, annotation = grangeslist)
proteomes=aastringsetlist
annotation=grangeslist

pdata <- process_input(proteomes, annotation)

pdata$annotation
pdata


data(blast_list)
if(diamond_is_installed()) {
  blast_list <- run_diamond(seq = pdata$seq)
}

blast_list

names(blast_list)

net <- infer_syntenet(blast_list, pdata$annotation, anchors = 2)
net4 <- infer_syntenet(blast_list, pdata$annotation, anchors = 4)


head(net)
id_table <- create_species_id_table(names(proteomes))
clusters <- cluster_network(net4)
clusters
profiles <- phylogenomic_profile(clusters)

species_order <- setNames(
  # vector elements
  c("Pmm","Pc","Lt","Ps","Pl","Pa","Gm","Gc","Oid","Epu","Pn","Epi","Px","EnE","Bgt","Bgh","Bga"),
  # vector names
  c("Pmm","Pc","Lt","Ps","Pl","Pa","Gm","Gc","Oid","Epu","Pn","Epi","Px","EnE","Bgt","Bgh","Bga")
)
species_order

species_annotation <- data.frame(
  Species = species_order,
  Family = c(
    rep("Dicot",14),rep("Monocot",3)
  )
)

species_annotation

species_annotation$Species
species_annotation$Family <- paste0(species_annotation$Species, "_", species_annotation$Family)

head(species_annotation)

species_annotation


# phylogenomic profiles, but using Ruzicka distances
pdf("~/Desktop/Pma_data/orthofinder_output/Pm_data/heatmap_annot2_new.pdf", width = 6, height = 8)
plot_profiles(
  profiles, 
  species_annotation, 
  cluster_species = species_order, 
  dist_function = labdsv::dsvdis,
  dist_params = list(index= "ruzicka")
)
dev.off()

# Find group-specific clusters
gs_clusters <- find_GS_clusters(profiles, species_annotation)

head(gs_clusters)
head(gs_clusters2)

# How many family-specific clusters are there?
nrow(gs_clusters)
nrow(gs_clusters2)
nrow(gs_clusters3)

# Filter profiles matrix to only include group-specific clusters
idx <- rownames(profiles) %in% gs_clusters$Cluster
p_gs <- profiles[idx, ]

idx <- rownames(profiles) %in% gs_clusters2$Cluster
p_gs <- profiles[idx, ]
idx
p_gs


plot_profiles(
  p_gs, species_annotation, 
  cluster_species = species_order, 
  cluster_columns = TRUE
)
#species' uniqune cluster groups

gs_clusters$Cluster[gs_clusters$Group == "Bga_Monocot"]

gs_clusters3$Cluster[gs_clusters3$Group == "Monocot"]



clusters$Gene[which(clusters$Cluster %in% gs_clusters$Cluster[gs_clusters$Group == "Monocot"])]
Monocot_specific <- clusters$Gene[which(clusters$Cluster %in% gs_clusters$Cluster[gs_clusters$Group == "Monocot"])]

Monocot_specific <- clusters$Gene[which(clusters$Cluster %in% gs_clusters3$Cluster[gs_clusters3$Group == "Monocot"])]


Pma_unique = clusters$Gene[which(clusters$Cluster %in% gs_clusters$Cluster[gs_clusters$Group == "Bga_Monocot"])]
Pma_unique2 = clusters$Gene[which(clusters$Cluster %in% gs_clusters2$Cluster[gs_clusters2$Group == "Bga_Monocot"])]

Pma_unique2_clusters = clusters$Cluster[which(clusters$Cluster %in% gs_clusters2$Cluster[gs_clusters2$Group == "Bga_Monocot"])]

Pma_unique2
Pma_unique2_clusters

gs_clusters2$Cluster[gs_clusters2$Group == "Bga_Monocot"]

my_data_frame <- data.frame(Pma_unique)
my_data_frame2 <- data.frame(Pma_unique2)
my_data_frame3 <- data.frame(Monocot_specific)

write.table(my_data_frame, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/Pma_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)
write.table(Monocot_specific, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/Monocot_specific_id.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)
write.table(my_data_frame3, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/All_shared_id.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

my_clusters <- data.frame(clusters)
write.table(my_clusters, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/myclusters_all.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

#my_profiles <-data.frame(profiles)
#write.table(my_profiles, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/myprofiles.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

#ids <- read.table("~/Desktop/Pma_data/orthofinder_output/Pm_data/pre_results/nyclusters_sorted.txt")
#ids2 <- read.table("~/Desktop/Pma_data/orthofinder_output/Pm_data/pre_results/Pma_Pmt_shared_corted_id.txt")
#tmp <- ids[which(ids$V2 %in% ids2$V1), ]
#write.table(tmp, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/shared.gene.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

my_gs_clusters2_clusters = data.frame(gs_clusters2$Cluster[gs_clusters2$Group == "Bga_Monocot"])
write.table(my_gs_clusters2_clusters, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/myclusters_Bga_monocot.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

ids <- read.table("~/Desktop/Pma_data/orthofinder_output/Pm_data/myclusters_Bga_monocot.txt")
ids

ids2 <- read.table("~/Desktop/Pma_data/orthofinder_output/Pm_data/myclusters_all.txt")
ids2
tmp <- ids[which(ids$V1 %in% ids2$V2), ]
tmp
write.table(tmp, file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/shared.gene.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)

all_clusters=data.frame(clusters)
all_clusters

write.table(all_clusters,file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/all_clusters.id.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)
write.table(net4,file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/all_net_anchor4.id.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)
write.table(net,file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/all_net_anchor2.id.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)


id <- clusters$Cluster[8000:8100]
id
id <- gs_clusters2$Cluster

plot_network(net, clusters, cluster_id = Pma_secreted_dataframe$V2[1:100], color_by = gene_df, interactive=FALSE, dim_interactive = c(1000, 1000))

plot_network(net, clusters, cluster_id = id, color_by = gene_df)

Pma_secreted_id <- read.table("~/Desktop/Pma_data/Fasta_and_Annotation_files/Pma_final_secreted_syntenet_ids.txt")
Pma_secreted_dataframe=data.frame(Pma_secreted_id)

#genes <- unique(c(net$Anchor1, net$Anchor2))

gene_df <- data.frame(
  Gene = net,
  Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)

net

gene_df <- merge(gene_df, species_annotation)[, c("Gene", "Family")]
write.table(gene_df,file = "~/Desktop/Pma_data/orthofinder_output/Pm_data/all_gene_family.id.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = F)


head(gene_df)

plot_network(
  net, clusters, cluster_id = Pma_secreted_dataframe$V2, 
  interactive = TRUE, dim_interactive = c(1000, 600), color_by = gene_df
)


bt_mat <- binarize_and_transpose(profiles)
bt_mat[1:17, 1:17]

if(iqtree_is_installed()) {
    phylo <- infer_microsynteny_phylogeny(bt_mat, 
                                          threads = 1)
}

Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), "~Desktop/Scripts/iqtree-1.6.12-MacOSX/bin/", sep = ":"
  )
)

library(ggtree)
library(treeio)

phylo

tree=read.tree("/var/folders/kf/4xn__yc93d35stj58m8nc46m0000gn/T//RtmpUgeFtu/microsynteny_phylogeny_29_Nov_2023_23h43.phy.treefile")

ggtree(tree_net4) +
  geom_tiplab(size = 3)


