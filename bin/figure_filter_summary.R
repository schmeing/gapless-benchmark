suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

args = commandArgs(trailingOnly=TRUE)
if(length(args)<4){
  stop("The four filter folders with the performance csv files and the output pdf must be supplied", call.=FALSE)
} else if(length(args)>1){
  stop("Only the four filter folders with the performance csv files and the output pdf must be supplied", call.=FALSE)
}
per_folder1 <- args[1]
per_folder2 <- args[2]
per_folder3 <- args[3]
per_folder4 <- args[4]
out_file <- args[5]
#per_folder1 <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-pacbio_clr_21x-supernova-chromium"
#per_folder2 <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-hifi_33x-supernova-T2T_10X_NovaSeq"
#per_folder3 <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-nanopore_30x-supernova-T2T_10X_NovaSeq"
#per_folder4 <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-nanopore_30x-flye-hifi_33x"
#out_file <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/minigap_benchmark/results/figures/supplementary/figure_filter_performance.pdf"

per_csv <- read_csv(paste0(per_folder4,'/performance.csv'), col_types = cols(), na = c("", "NA")) %>%
  rename(fdr4=fdr, completeness4=completeness)
per_csv <- read_csv(paste0(per_folder3,'/performance.csv'), col_types = cols(), na = c("", "NA")) %>%
  rename(fdr3=fdr, completeness3=completeness) %>%
  full_join(per_csv, by=c("fabs","fprem","frel","min_len"))
per_csv <- read_csv(paste0(per_folder2,'/performance.csv'), col_types = cols(), na = c("", "NA")) %>%
  rename(fdr2=fdr, completeness2=completeness) %>%
  full_join(per_csv, by=c("fabs","fprem","frel","min_len"))
per_csv <- read_csv(paste0(per_folder1,'/performance.csv'), col_types = cols(), na = c("", "NA")) %>%
  rename(fdr1=fdr, completeness1=completeness) %>%
  full_join(per_csv, by=c("fabs","fprem","frel","min_len"))

performance <- per_csv %>%
  mutate(fabs = as.factor(fabs), fprem=as.factor(fprem), frel=as.factor(frel), min_len=as.factor(min_len)) %>%
  mutate(Filter = if_else((fabs=="2") & (fprem=="0.001"), "Relative counts", "Unspecified") ) %>%
  mutate(Filter = if_else((fabs=="2") & (frel=="1e+50"), "Truncated distribution", Filter) ) %>%
  mutate(Filter = if_else((fprem=="0") & (frel=="1e+50"), "Absolute counts", Filter) ) %>%
  mutate(Filter = if_else(is.na(min_len), Filter, "Minimum Length")) %>%
  mutate(Filter = factor(Filter, levels=c("Unspecified","Minimum Length","Absolute counts","Truncated distribution","Relative counts"))) %>%
  mutate(Selected = ((fabs=="2") & (fprem=="0.001") & (frel=="10")) | 
                    ((fabs=="2") & (fprem=="0.001") & (frel=="1e+50")) |
                    ((fabs=="2") & (fprem=="0") & (frel=="1e+50")) |
                    (min_len == "500") ) %>%
  mutate(Selected = replace_na(Selected, FALSE)) %>%
  arrange(Selected,Filter)

text_size <- 20
tick_text_size <- 16

design <- list(
  geom_point(size=6),
  scale_colour_manual(values=c("grey","#D92120","#781C81","#664CFF","#7FB972")),
  xlab("FDR"),
  ylab("Completeness"),
  theme_bw(),
  theme(axis.text.x = element_text( size = tick_text_size),
        axis.text.y = element_text( size = tick_text_size),
        axis.title.x = element_text( size = text_size),
        axis.title.y = element_text( size = text_size),
        strip.text.x = element_text( size = text_size),
        strip.text.y = element_text( size = text_size),
        legend.title=element_text(size=text_size), 
        legend.text=element_text(size=tick_text_size),
        plot.title = element_text(hjust = 0.5, size=text_size),
        panel.grid.minor = element_blank())
)

p1 <- performance %>%
  ggplot(aes(x=fdr1, y=completeness1, color=Filter, shape=Selected)) +
    ggtitle("Dolphin PacBio CLR 21x") +
    design
p2 <- performance %>%
  ggplot(aes(x=fdr2, y=completeness2, color=Filter, shape=Selected)) +
  ggtitle("Human PacBio HiFi 33x") +
  design
p3 <- performance %>%
  ggplot(aes(x=fdr3, y=completeness3, color=Filter, shape=Selected)) +
  ggtitle("Human Nanopore sn 30x") +
  design
p4 <- performance %>%
  ggplot(aes(x=fdr4, y=completeness4, color=Filter, shape=Selected)) +
  ggtitle("Human Nanopore Flye 30x") +
  design

print((p1 + p2) / (p3 + p4) + plot_layout(guides = "collect"))

ggsave(out_file, width=297, height=210, units="mm")
