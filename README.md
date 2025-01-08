JASPAR_MOTIFscan: scans MOTIFs from JASPAR on inquire seuquences regardless of enrichment
MOTIF_score: Use the MOTIFscan csv file from JASPAR_MOTIFscan, and Sum the scores and counts based on MOTIFs
MOTIF_score_enrich: Use the MOTIFscan RDS file, can do enrichment analysis between groups, e.g. control and KO, negative and positive. Out put includes a "enrichment ratio", which is sum(rel_Score_group1)/sum (rel_Score_group2). 
