library(bio3d)
library(tidyverse)
library(RColorBrewer)

standard_SASA <- read_tsv("per_AA_standard_ASA_CWilke.txt", col_names = T)
switch_loops <- tibble(resi = c(19:25, 41:47, 69:77), "loop" = c(rep("ploop", 7), rep("switch1", 7), rep("switch2",9)))
pdb_dir <- "pdbs/clean"
chainA_dir <- "pdbs/clean_chainA/"
files_list <- list.files(path = pdb_dir)
all_chains_file_paths <- file.path(pdb_dir, files_list)
chainA_file_paths <- file.path(chainA_dir, files_list)
nucleotide_contacts_table <- data.frame()
interface_SASA_table <- tibble(file_path = character(), resno = integer(), resid = character(),
              SASA_difference = double(), percent_SASA_change = double(), interface = character())

for (i in seq_along(all_chains_file_paths)) {
  file_path <- all_chains_file_paths[i]
  chainA_path <- chainA_file_paths[i]
  pdb <- read.pdb(file_path, rm.alt = FALSE)
  pdbA <- read.pdb(chainA_path, rm.alt = FALSE)
  SASA_monomer <- dssp(pdbA)$acc
  SASA_complex <- dssp(pdb)$acc[1:length(SASA_monomer)]
  coord.temp <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
    filter(resid != "HOH" & resid != "SO4" & resid != "MG" & 
             resid != "GNP" & resid != "GDP"  & resid != "GTP" & resid != "AF3" & resid != "GOL" &
             resid != "EDO" & resid != "CL")
  chainA_deltaSASA <- coord.temp %>% 
    filter(chain == "A") %>% 
    select(resno, resid) %>% unique() %>% 
    add_column(., SASA_monomer, SASA_complex, file_path) %>%
    inner_join(., standard_SASA, by = c("resid" = "AA")) %>% 
    select(-A, -Residue) %>% 
    mutate("rASAm" = SASA_monomer/maxASA, "rASAc" = SASA_complex/maxASA) %>% 
    mutate("deltarASA" = rASAm - rASAc) %>% 
    filter(deltarASA > 0) %>% 
    mutate("interface" = ifelse((rASAm >= 0.25 & rASAc < 0.25 & deltarASA >= 0.25), "core", "rim")) %>% 
    mutate("interface" = ifelse(rASAm < 0.25, "support", interface)) %>% 
    select(file_path, resno, resid, deltarASA, interface)
  interface_SASA_table <- interface_SASA_table %>% bind_rows(., chainA_deltaSASA)
  nucleotide_coord <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
    filter(resid != "HOH" & resid != "SO4" & resid != "AF3" & resid != "GOL" &
             resid != "EDO" & resid != "CL" & (resid == "MG" | resid == "GNP" | resid == "GDP" | resid == "GTP")) %>% 
    mutate("unique" = str_c(chain, resno, elety, sep = "_"))
  if (nrow(nucleotide_coord)) {
  nucleotide_coord_matrix <- as.matrix(nucleotide_coord[, c("x", "y", "z")])
    rownames(nucleotide_coord_matrix) <- nucleotide_coord$unique
    Gsp1_coord <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
      filter(chain == "A" & resid != "HOH" & resid != "SO4" & resid != "MG" & 
               resid != "GNP" & resid != "GDP"  & resid != "GTP" & resid != "AF3" & resid != "GOL" &
               resid != "EDO" & resid != "CL") %>% 
      mutate("unique" = str_c(chain, resno, elety, sep = "_"))
    Gsp1_coord_matrix <- as.matrix(Gsp1_coord[, c("x", "y", "z")])
    rownames(Gsp1_coord_matrix) <- Gsp1_coord$unique
    distances <- dist.xyz(Gsp1_coord_matrix, nucleotide_coord_matrix)
    colnames(distances) <- rownames(nucleotide_coord_matrix)
    distances_df <- data.frame(distances)
    distances_df <- cbind(distances_df, "unique_Gsp1" = as.character(rownames(Gsp1_coord_matrix)))
    nucleotide_contacts <- as_tibble(distances_df) %>% 
      gather(unique_nucleotide, dist, -unique_Gsp1) %>%
      filter(dist < 5) %>%     ##### threshold for nucleotide contact is 5 A
      separate(col = unique_Gsp1, into = c("Gsp1", "resnum1", "atom1"), sep = "_", convert = T) %>%
      separate(col = unique_nucleotide, into = c("nucleotide", "resnum2", "atom2"), sep = "_", convert = T) %>%
      group_by(resnum1, resnum2) %>%
      summarise("number_of_contacts" = n()) %>%
      arrange(resnum1, resnum2)
    if (nrow(nucleotide_contacts) > 0) {
      nucleotide_contacts_table <- rbind(nucleotide_contacts_table, 
                              as_tibble(data.frame(file_path, nucleotide_contacts)))
    }
  }
  # chains <- unique(coord.temp$chain)
  # coords <- list()
  # for (ch in seq_along(chains)) {
  #   chain_id <- chains[ch]
  #   coords[[chain_id]] <- coord.temp %>% 
  #     filter(chain == chain_id) %>% 
  #     mutate("unique" = str_c(chain, resno, elety, sep = "_"))
  # }
  # 
  #pairs_of_chains <- combn(chains, 2)
  # for (p in seq_along(chains)) {
  #   chain1 <- "A"
  #   chain2 <- chains[p]
  #   if (chain2 != "A") {
  #     matrix1 <- as.matrix(coords[[chain1]][, c("x", "y", "z")])
  #     rownames(matrix1) <- coords[[chain1]]$unique
  #     matrix2 <- as.matrix(coords[[chain2]][, c("x", "y", "z")])
  #     rownames(matrix2) <- coords[[chain2]]$unique
  #     distances <- dist.xyz(matrix1, matrix2)
  #     colnames(distances) <- rownames(matrix2)
  #     distances_df <- data.frame(distances)
  #     distances_df <- cbind(distances_df, "unique_chain1" = as.character(rownames(matrix1)))
  #     contacts_df <- as_tibble(distances_df) %>% 
  #       gather(unique_chain2, dist, -unique_chain1) %>%
  #       filter(dist < 4) %>%
  #       separate(col = unique_chain1, into = c("chain1", "resnum1", "atom1"), convert = T) %>%
  #       separate(col = unique_chain2, into = c("chain2", "resnum2", "atom2"), convert = T) %>%
  #       group_by(resnum1, resnum2) %>%
  #       summarise("number_of_contacts" = n()) %>%
  #       arrange(resnum1, resnum2)
  #     if (nrow(contacts_df) > 0) {
  #       contacts_table <- rbind(contacts_table, as_tibble(data.frame(file_path, chain1, chain2, contacts_df)))
  #     }
  # }
  #}
}

#contacts_table <- as_tibble(contacts_table)
#write_delim(contacts_table, "4A_contacts.txt", delim = "\t")
index <- read_tsv("index.txt", col_names = T)
Ran_seq_index <- read_tsv("Gsp1_to_Ran.index", col_names = T)

### residues contacting the nucleotide (based on the 6 A threshold)
nuc_merged <- inner_join(nucleotide_contacts_table, index, by = c("file_path"))
nuc_merged_yeast <- nuc_merged %>%
  filter(species == "yeast") %>%
  select(partner, "yeastresnum" = resnum1, number_of_contacts, pdb_id)
nuc_merged_mamm <- nuc_merged %>%
  filter(species == "mammalian") %>%
  inner_join(., Ran_seq_index, by = c("resnum1" = "pdb_seq_num")) %>%
  select(partner, "yeastresnum" = ref_seq_num, number_of_contacts, pdb_id)
nucleotide_currated <- bind_rows(nuc_merged_yeast, nuc_merged_mamm)
write_tsv(nucleotide_currated, "Gsp1_residues_in_contact_with_nucleotides.txt")
#############################################################################


######## Gsp1 interfaces based on delta_rASA
merged_SASA <- interface_SASA_table %>% 
  inner_join(., index, by = "file_path") %>% 
  select(-chain)
merged_yeast_SASA <- merged_SASA %>% 
  filter(species == "yeast") %>% 
  select(partner, "yeastresnum" = resno, deltarASA, interface, pdb_id) %>% 
  mutate("deltarASA" = round(deltarASA, 2))
merged_mamm_SASA <- merged_SASA %>%
  filter(species == "mammalian") %>%  
  inner_join(., Ran_seq_index, by = c("resno" = "pdb_seq_num")) %>% 
  select(partner, "yeastresnum" = ref_seq_num, deltarASA, interface, pdb_id) %>% 
  mutate("deltarASA" = round(deltarASA, 2))
currated_SASA <- bind_rows(merged_yeast_SASA, merged_mamm_SASA) %>% 
  group_by(yeastresnum, partner) %>% 
  filter(deltarASA == max(deltarASA, na.rm = T)) %>% 
  filter(deltarASA > 0.05) %>% 
  select(-pdb_id) %>% 
  arrange(interface, yeastresnum, partner) %>% 
  unique() %>% 
  ungroup()
write_tsv(currated_SASA, "SASA_interfaces.txt")


#### merge yeast constructs of gsp1 mutants with the interface information
constructs <- read_tsv("constructs_table.txt")
merge <- currated_SASA %>%  inner_join(., constructs, by = "yeastresnum") %>% 
  mutate("interface_summary" = str_c(interface, deltarASA, sep = " / ")) %>% 
  select(-interface, -deltarASA)
constructs_interface_summary <- merge %>% 
  mutate("temp" = str_c(yeastresnum, `construct name`, sep = " / ")) %>% 
  select(temp, partner, interface_summary) %>% 
  spread(partner, interface_summary) %>% 
  separate(temp, into = c("yeastresnum", "construct name"), sep = " / ") %>% 
  mutate(yeastresnum = as.integer(yeastresnum)) %>% 
  arrange(yeastresnum)
write_tsv(constructs_interface_summary, "Gsp1_constructs_interface_summary.txt", na = "")


### cluster Gsp1 residues by interface delta_rASA
cluster_constructs_by_interface <- currated_SASA %>% 
  inner_join(constructs, by = "yeastresnum") %>% 
  select(`Gsp1 point mutation`, deltarASA, partner) %>% 
  unique() %>% 
  complete(`Gsp1 point mutation`, partner, fill = list(deltarASA = 0)) %>% 
  arrange(`Gsp1 point mutation`) #%>% 
  ##spread(partner, deltarASA)

get_dd <- function(data) {
  cormat <- cor(data[,-1], use = "pairwise.complete.obs")
  #dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  #cormat[is.na(cormat)] <- 0
  return(dd)
}

interface_heatmap <- function(data) {
  mutants <- data %>% 
    pull(`Gsp1 point mutation`) %>% unique()
  partners <- data %>% 
    pull(partner) %>% unique()
  mut_spread <- data %>% 
    spread(`Gsp1 point mutation`, deltarASA)
  dd <- get_dd(mut_spread)
  ordered_mutants <- mutants[hclust(dd)$order]
  partner_spread <- data %>% 
    spread(partner, deltarASA)
  dd <- get_dd(partner_spread)
  ordered_partners <- partners[hclust(dd)$order]
  data %>% 
    mutate("mutant" = factor(`Gsp1 point mutation`, ordered_mutants),
           "partner" = factor(partner, ordered_partners)) %>% 
    ggplot(aes(x = mutant, y = partner, fill = deltarASA)) +
    geom_tile() +
    scale_fill_gradient(low = "ivory", high = "seagreen4") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
interface_heatmap(cluster_constructs_by_interface)
ggsave("deltarASA_heatmap.pdf", width = 12)


#### attempted yeast transformations
res <- sort(c(34,58,78,79,80,84,101,102,105,108,112,115,129,132,137,139,141,143,147,148,154,157,169,180)) #,39,70,73,75))

### make a geom_tile plot that shows a union of all the interface residues and nucleotide contact residues, 
### and then all of those that I tried to make mutants of to show our coverage
### combine currated_SASA, nucleotide_currated and res
res_close_to_nuc <- nucleotide_currated %>% arrange(yeastresnum) %>% pull(yeastresnum) %>% unique()
combined_residue_data <- currated_SASA %>% 
  complete(yeastresnum, partner, fill = list(deltarASA = 0, interface = "none")) %>% 
  mutate("nucleotide" = ifelse(yeastresnum %in% res_close_to_nuc, "close_to_nucleotide", "not_close_to_nucleotide")) %>% 
  mutate("mutated" = ifelse(yeastresnum %in% res, "mutated", "not mutated")) %>% 
  mutate("switch" = ifelse(yeastresnum %in% switch_loops$resi, "switch", "outside of switch"))
  
res_that_are_core_at_least_once <- combined_residue_data %>% 
  filter(interface == "core") %>% 
  pull(yeastresnum) %>% unique()
combined_residue_data <- combined_residue_data %>% 
  filter(yeastresnum %in% res_that_are_core_at_least_once)
yeast_res <- combined_residue_data %>% 
  pull(yeastresnum) %>% unique()
partners <- combined_residue_data %>% 
  pull(partner) %>% unique()
res_spread <- combined_residue_data %>% 
  spread(yeastresnum, deltarASA) %>% select(-interface, -nucleotide, -mutated, -switch)
dd <- get_dd(res_spread)
ordered_res <- yeast_res[hclust(dd)$order]
partner_spread <- combined_residue_data %>% 
  select(-interface, -nucleotide, -mutated, -switch) %>% 
  unique() %>% 
  spread(partner, deltarASA)
dd <- get_dd(partner_spread)
ordered_partners <- partners[hclust(dd)$order]


cols <- c("not mutated" = "#fdb863", "mutated" = "#5e3c99")
shape <- c("not_close_to_nucleotide" = 16, "close_to_nucleotide" = 17)
size <- c("core" = 3, "rim" = 2, "support" = 2, "none" = 0.5)
alpha <- c("outside of switch" = 1, "switch" = 0.3)

combined_residue_data %>% 
  mutate("yeastresnum" = factor(yeastresnum, ordered_res),
         "partner" = factor(partner, ordered_partners)) %>% 
  ggplot(aes(x = yeastresnum, y = partner)) + 
  geom_point(aes(size = interface, shape = nucleotide, color = mutated, alpha = switch)) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shape) +
  scale_alpha_manual(values = alpha) +
  scale_size_manual(values = size) +
  xlab("Gsp1 residue number") + ylab("Gsp1 binding partner") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

ggsave("Overview_of_residue_selection_for_Gsp1_PM.pdf", width = 14, height = 8)

cols <- c("core" = "#e66101", "rim" = "#fdb863", "support" = "#5e3c99", "none" = "#b2abd2")
shape <- c("not_close_to_nucleotide" = 16, "close_to_nucleotide" = 17)
alpha <- c("mutated" = 1, "not mutated" = 0.3)
size <- c(c("mutated" = 2, "not mutated" = 1))




#merged <- inner_join(index, currated, by = c("chain" = "chain2", "structure" = "file_path"))
# interface_contact_n <- currated %>%
#   group_by(pdb_id, partner, yeastresnum) %>% 
#   summarise("per_resi_n_contacts" = round(sum(number_of_contacts), 2)) %>%
#   ungroup() %>% 
#   group_by(partner, yeastresnum) %>%
#   summarise("mean_n_contacts" = round(mean(per_resi_n_contacts), 2)) %>%
#   arrange(yeastresnum, desc(mean_n_contacts)) %>%
#   ungroup()
# 
# #contacts_df <- inner_join(contacts_table, index, by = c("file_path", "chain2" = "chain"))
# interface_size <- interface_contact_n %>%
#   group_by(partner) %>%
#   summarize("avg_interface_size" = round(sum(mean_n_contacts), 2)) %>%
#   arrange(partner) %>%
#   ungroup() 
# 
# interface_contact_n <- inner_join(interface_contact_n, interface_size, by = "partner") %>%
#   mutate("relative_contact" = round(mean_n_contacts / avg_interface_size, 4)) %>% 
#   mutate("percent_interface" = round(relative_contact * 100, 2))
# 
# write_tsv(interface_contact_n, path = "4A_contacts_and_interfaces.txt")
# 
# interface_contact_n %>% 
#   rename( `mean number of heavy\natoms within\n6 A from the partner` = mean_n_contacts,
#          `% of total contacts\nin the interface` = percent_interface) %>% 
#   ggplot(aes(partner, yeastresnum, size = `% of total contacts\nin the interface`, fill = `mean number of heavy\natoms within\n6 A from the partner`)) +
#   geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Gsp1/Ran interface residues") +
#   ylab("Gsp1 residue number")
# ggsave("4A_Gsp1_interfaces_point_plot.pdf", width = 12)
# 
nucleotide_contacts_table <- nucleotide_contacts_table %>%
  filter(number_of_contacts > 1) %>% 
  inner_join(., index, by = "file_path")
merged_mamm_nuc <- nucleotide_contacts_table %>%
  filter(species == "mammalian") %>%  
  inner_join(., Ran_seq_index, by = c("resnum1" = "pdb_seq_num")) %>% 
  select(partner, "yeastresnum" = ref_seq_num, number_of_contacts, pdb_id)
merged_yeast_nuc <- nucleotide_contacts_table %>%
  filter(species == "yeast") %>%  
  inner_join(., Ran_seq_index, by = c("resnum1" = "pdb_seq_num")) %>% 
  select(partner, "yeastresnum" = ref_seq_num, number_of_contacts, pdb_id)
nuc_binding_res <- bind_rows(merged_yeast_nuc, merged_mamm_nuc) %>% 
  arrange(yeastresnum) %>% 
  pull(yeastresnum) %>% unique()
switch_loops_and_nuc_binding <- switch_loops %>% 
  bind_rows(., tibble("resi" = nuc_binding_res, "loop" = "nuc_binding"))
interface_and_nuc_loop <- currated_SASA %>% 
  inner_join(., switch_loops_and_nuc_binding, by = c("yeastresnum" = "resi"))
### interface point plot - using deltaSASA
currated_SASA <- currated_SASA %>% 
  filter(partner != "NUP1" & partner != "NUP60")
core <- currated_SASA %>% 
  filter(interface == "core")
rim_support <- currated_SASA %>% 
  filter(interface != "core")
currated_SASA %>% 
  ggplot(aes(partner, yeastresnum, size = deltarASA), shape = 21, stroke = 0.1) +
  #geom_point(data = rim_support, alpha = 0.2, color = "#E69F00") +
  theme_bw(base_size = 14) +
  geom_point(data = core, color = "#D55E00") +
  ylab("Gsp1 residue number") + xlab("partner protein") +
  theme(axis.text.x = element_text(angle = -30, hjust = 0.5, vjust = 0.5)) +
  NULL
ggsave("interface_core_residues.pdf", height = 5, width = 10)
currated_SASA %>% 
  ggplot(aes(partner, yeastresnum, size = deltarASA), shape = 21, stroke = 0.1) +
  #geom_point(data = rim_support, alpha = 0.2, color = "#E69F00") +
  geom_point(data = core, color = "#D55E00") +
  geom_point(data = interface_and_nuc_loop, color = "gray", shape = 19) +
  NULL

