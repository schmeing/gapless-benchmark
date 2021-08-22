suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
  stop("The input and output file must be supplied", call.=FALSE)
} else if(length(args)>2){
  stop("Only the input and output file must be supplied", call.=FALSE)
}
comp_res_csv <- args[1]
outfile <- args[2]

comp_res <- read_csv(comp_res_csv, col_types = cols(), na = c("", "NA", "-"))

data <- comp_res %>%
  mutate(cpu_time=cpu_time_s/3600, elapsed_time=elapsed_time_s/3600, memory=memory_kb/1024/1024) %>%
  select(-cpu_time_s,-elapsed_time_s,-memory_kb) %>%
  mutate(call=factor(call, levels=c("gapless_split","minimap2_repeats","minimap2_reads","gapless_scaffold","minimap2_extension","gapless_extend","gapless_finish","minimap2_consensus","racon"))) %>%
  mutate(call=recode(call, "gapless_split"="split", "minimap2_repeats"="mm2-rep", "minimap2_reads"="mm2-seq", "gapless_scaffold"="scaffold", "minimap2_extension"="mm2-ext", "gapless_extend"="extend", "gapless_finish"="finish", "minimap2_consensus"="mm2-rac")) %>%
  gather(measure,value,-pass,-call) %>%
  mutate(measure=recode(measure, "cpu_time"="CPU time (h)", "elapsed_time"="Elapsed time (h)", "memory"="Max. memory (Gb)"))

text_size <- 20
tick_text_size <- 16

data %>%
  ggplot(aes(x=call, y=value, color=pass, fill=pass)) +
    geom_col(position = "dodge") +
    facet_grid(measure ~ ., scales="free_y", switch="y") + 
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = tick_text_size),
          strip.text.y = element_text( size = tick_text_size),
          legend.title=element_blank(),
          legend.text=element_text(size=tick_text_size),
          legend.margin=margin(),
          legend.position = c(0.05, 0.94),
          panel.grid.minor = element_blank())


ggsave(outfile, width=297, height=210, units="mm")
