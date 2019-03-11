library(tidyverse)
yeast_seq <- read_tsv("Ran_conservation/yeast_seq.txt", col_names = F)
n <- nchar(yeast_seq$X2)
yeast_seq <- yeast_seq %>% 
  mutate(X2 = gsub('(?<=.)(?=.)', ',', X2, perl = TRUE)) %>% 
  separate(X2, into = paste0('res', 1:n), sep=",")
yeast_seq <- yeast_seq %>% separate(X2, str_c("res", 1:n), sep = "", remove = T)
uniprot_aln <- read_tsv("Ran_conservation/uniprot_alignment_seq.txt", col_names = F)
uniprot_aln <- uniprot_aln %>% 
  mutate(X2 = gsub('(?<=.)(?=.)', ',', X2, perl = TRUE)) %>% 
  separate(X2, into = paste0('res', 1:n), sep=",")
align <- bind_rows(yeast_seq, uniprot_aln) %>% 
  rename("seq" = X1) %>% 
  gather(., "res", "aa", -seq)
resi_to_keep <- align %>% 
  filter(seq == "yeast" & aa != "0") %>% 
  pull(res)
aln_to_keep <- align %>% 
  filter(res %in% resi_to_keep) %>% 
  spread(res, aa)
zoom_resi_to_keep <- resi_to_keep[100:110]
zoom_aln_to_keep <- align %>% 
  filter(res %in% zoom_resi_to_keep) %>% 
  spread(res, aa)
write_tsv(aln_to_keep, path = "Ran_conservation/chopped_aln.txt")
write_tsv(zoom_aln_to_keep, path = "Ran_conservation/chopped_aln_100-110.txt")

empiric <- read_tsv("Ran_conservation/BolonData.txt", col_names = F)
oldnames = str_c("X", 2:ncol(empiric))
newnames = 2:ncol(empiric)
aa_order <- empiric$X1
empiric <- empiric %>% 
  rename("aa" = X1) %>% 
  rename_at(vars(oldnames), ~ newnames) %>% 
  gather("residue", "score", -aa) %>% 
  mutate("residue" = as.numeric(residue)) %>% 
  filter(residue > 13 & residue < 170) %>% 
  mutate("aa" = factor(aa, aa_order)) %>% 
  arrange(residue, aa)
max_stop <- empiric %>% 
  filter(aa == "STOP") %>% 
  summarise("max" = max(score)) %>% pull(max)
empiric <- empiric %>% 
  mutate("fitness score" = ifelse(score <= max_stop, -1.5, score))
empiric %>% ggplot(aes(residue, aa, fill = `fitness score`)) + 
  geom_tile() + scale_fill_gradient2(low = "#0072B2", high = "#F0E442") +
  theme_bw() + ylab("amino acid") +
  theme_bw(base_size = 30)
  
  
ggsave("Ran_conservation/EMPIRICHeatmap.pdf", width = 25)

empiric_zoom <- empiric %>% 
  filter(residue > 100 & residue < 110) %>% 
  mutate("residue" = as.character(residue))
wt_resi <- tibble("resi" = 100:110, "aa" = str_split("YKNVPNWHRDL", pattern = "")[[1]]) %>% 
  mutate("temp" = str_c(resi, aa, sep = ""))
empiric_zoom_wt <- empiric_zoom %>% 
  mutate("temp" = str_c(residue, aa, sep = "")) %>% 
  filter(temp %in% wt_resi$temp)
empiric_zoom %>% ggplot(aes(residue, aa, fill = `fitness score`)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0072B2", high = "#F0E442") +
  #geom_tile(data = empiric_zoom_wt, aes(residue, aa, color = "black"), size = 1, width = 0.95, height = 0.95) +
  theme_bw() + ylab("amino acid") +
  theme_bw(base_size = 20)

ggsave("Ran_conservation/EMPIRIC_resi100-110.pdf", width = 10)
