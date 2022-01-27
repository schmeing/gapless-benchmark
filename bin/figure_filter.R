suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  stop("The filter folder with the performance csv file and the plot title must be supplied", call.=FALSE)
} else if(length(args)>2){
  stop("Only the filter folder with the performance csv file and the plot title must be supplied", call.=FALSE)
}
per_folder <- args[1]
title <- args[2]
#per_folder <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-hifi_33x-supernova-T2T_10X_NovaSeq"
#title <- "Human PacBio HiFi 33x"
#per_folder <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-nanopore_30x-flye-hifi_33x"
#title <- "Human Nanopore Flye 30x"
#per_folder <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-nanopore_30x-supernova-T2T_10X_NovaSeq"
#title <- "Human Nanopore supernova 30x"
#per_folder <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/own_scaffolding_minimap2/filter_test-pacbio_clr_21x-supernova-chromium"
#title <- "Dolphin PacBio CLR 21x"

per_csv <- read_csv(paste0(per_folder,'/performance.csv'), col_types = cols(), na = c("", "NA"))

performance <- per_csv %>%
  filter(fabs > 0 | is.na(fabs)) %>%
  mutate(fabs = as.factor(fabs), fprem=as.factor(fprem), frel=as.factor(frel), min_len=as.factor(min_len))

text_size <- 20
tick_text_size <- 16

design <- list(
  geom_point(size=6),
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
        panel.grid.minor = element_blank()),
  guides(color=guide_legend(title=""))
)

p1 <- performance %>%
  ggplot(aes(x=fdr, y=completeness, color=fabs)) +
    design

p1 + ggtitle(title) + guides(color=guide_legend(title="Absolute counts"))

ggsave(paste0(per_folder,'/filter_fabs.pdf'), width=297, height=210, units="mm")

p2 <- performance %>%
  ggplot(aes(x=fdr, y=completeness, color=fprem)) +
  design

p2 + ggtitle(title) + guides(color=guide_legend(title="Truncated distribution"))

ggsave(paste0(per_folder,'/filter_fprem.pdf'), width=297, height=210, units="mm")

p3 <- performance %>%
  ggplot(aes(x=fdr, y=completeness, color=frel)) +
  design

p3 + ggtitle(title) + guides(color=guide_legend(title="Relative counts"))

ggsave(paste0(per_folder,'/filter_frel.pdf'), width=297, height=210, units="mm")

p4 <- performance %>%
  ggplot(aes(x=fdr, y=completeness, color=min_len)) +
  design

p4 + ggtitle(title) + guides(color=guide_legend(title="Minimum Length"))

ggsave(paste0(per_folder,'/filter_min_len.pdf'), width=297, height=210, units="mm")

(p4 + ggtitle("Minimum Length") + p1 + ggtitle("Absolute counts")) / (p2 + ggtitle("Truncated distribution") + p3 + ggtitle("Relative counts")) +
  plot_annotation(title = title, theme = theme(plot.title = element_text(hjust = 0.5, size=text_size)))

ggsave(paste0(per_folder,'/filter_performance.pdf'), width=297, height=210, units="mm")
