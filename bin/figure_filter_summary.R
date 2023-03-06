suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

args = commandArgs(trailingOnly=TRUE)
if(length(args)<9){
	  stop("The four filter performance csv files with the title (name title name title ...) and the output pdf must be supplied", call.=FALSE)
} else if(length(args)>9){
	  stop("Only the four filter performance csv files with the title (name title name title ...) and the output pdf must be supplied", call.=FALSE)
}

title1 <- args[2]
title2 <- args[4]
title3 <- args[6]
title4 <- args[8]
out_file <- args[9]

per_csv <- lapply(c(1,3,5,7), function(n_arg){
  read_csv(args[n_arg], col_types = cols(), na = c("", "NA")) %>%
  mutate(Dataset=args[n_arg+1])
}) %>%
  bind_rows()

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
  mutate(Dataset = factor(Dataset, levels=c(title1, title2, title3, title4))) %>%
  arrange(Dataset,Selected,Filter) %>%
  mutate(completeness = replace_na(completeness, 0))

performance_line <- performance %>%
  filter(Filter != "Unspecified") %>%
  mutate(copy="org")

performance_line <- bind_rows(list(performance_line, mutate(filter(performance_line, Selected), copy="copy"))) %>%
  mutate(Filter=if_else(copy=="copy" & Filter=="Truncated distribution", "Relative counts", as.character(Filter))) %>%
  mutate(Filter=if_else(copy=="copy" & Filter=="Absolute counts", "Truncated distribution", as.character(Filter))) %>%
  mutate(Filter=if_else(copy=="copy" & Filter=="Minimum Length", "Absolute counts", as.character(Filter))) %>%
  mutate(Filter = factor(Filter, levels=c("Unspecified","Minimum Length","Absolute counts","Truncated distribution","Relative counts"))) %>%
  arrange(Filter, min_len, fabs, fprem, frel)

text_size <- 20
tick_text_size <- 16

performance %>%
  ggplot(aes(x=fdr, y=completeness, color=Filter, shape=Selected)) +
    geom_vline(data=ungroup(filter(group_by(performance, Dataset), fdr==min(fdr))), mapping=aes(xintercept=fdr), size=1, color="black", linetype = "dashed") +
    geom_vline(data=filter(performance, Selected, Filter=="Relative counts"), mapping=aes(xintercept=fdr), size=1, color="black") +
    geom_point(size=5) +
    geom_path(data=performance_line, mapping=aes(group=Filter), size=1) +
    scale_colour_manual(values=c("grey","#D92120","#781C81","#664CFF","#7FB972")) +
    geom_point(data=filter(performance, Selected), shape=2, color="black", size=5) +
    facet_wrap(. ~ Dataset) +
    xlab("FDR") +
    ylab("Fraction Possibilities Kept") +
    theme_bw() +
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

ggsave(out_file, width=297, height=210, units="mm")

