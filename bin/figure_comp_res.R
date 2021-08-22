suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)
if(length(args)<4){
  stop("The two input and two output files must be supplied", call.=FALSE)
} else if(length(args)>4){
  stop("Only the two input and two output files must be supplied", call.=FALSE)
}
comp_res_csv <- args[1]
assembly_csv <- args[2]
outfile1 <- args[3]
outfile2 <- args[4]

comp_res <- read_csv(comp_res_csv, col_types = cols(), na = c("", "NA", "-"))
assemblies <- read_csv(assembly_csv, col_types = cols(), na = c("", "NA", "-"))

data <- assemblies %>%
  select(-Quast) %>%
  left_join(rename(comp_res, Folder=folder), by="Folder") %>%
  mutate(folder=if_else(str_detect(Folder,"02-.*_gapcloser"), str_replace(Folder,"02-.*_gapcloser","01"), '')) %>%
  left_join(rename(comp_res, cpu_time_s1=cpu_time_s, elapsed_time_s1=elapsed_time_s, memory_kb1=memory_kb), by="folder") %>%
  mutate(folder=if_else(str_detect(Folder,"0[12]-"), str_replace(str_replace(str_replace(Folder,"0[12]-.+-supernova","00-supernova"),"0[12]-.+-shassembly","00-shassembly"),"0[12]-.+-flye","00-flye"), '')) %>%
  left_join(rename(comp_res, cpu_time_s0=cpu_time_s, elapsed_time_s0=elapsed_time_s, memory_kb0=memory_kb), by="folder") %>%
  replace_na(list(cpu_time_s1=0, elapsed_time_s1=0, memory_kb1=0, cpu_time_s0=0, elapsed_time_s0=0, memory_kb0=0)) %>%
  mutate(cpu_time_s=cpu_time_s+cpu_time_s1+cpu_time_s0, elapsed_time_s=elapsed_time_s+elapsed_time_s1+elapsed_time_s0, memory_kb=pmax(memory_kb,memory_kb1,memory_kb0)) %>%
  select(-folder,-cpu_time_s1,-elapsed_time_s1,-memory_kb1,-cpu_time_s0,-elapsed_time_s0,-memory_kb0) %>%
  mutate(memory_gb = memory_kb/1024/1024, cpu_time = if_else(Species == "E. Coli", cpu_time_s/60, cpu_time_s/3600/24), elapsed_time = if_else(Cores == 48, if_else(Species == "E. Coli", elapsed_time_s/60, elapsed_time_s/3600/24), as.numeric(NA))) %>%
  mutate(Coverage = replace_na(Coverage, "NA"), nrows=if_else((Assembly == "supernova" | (Assembly == "Flye" & Coverage == "high" & Data == "PacBio HiFi")) & Species == "Human", 2, 1)) %>%
  uncount(nrows, .id="id") %>%
  mutate(Data=if_else(Data != "Illumina" & Data != "10X Chromium", Data, if_else(Species != "Human", "PacBio CLR", if_else(id == 1, "PacBio HiFi", "Nanopore")))) %>%
  mutate(Data=if_else(id == 1, Data, "Nanopore"), Coverage=factor(if_else(id == 1, Coverage, "NA"), levels=c("very high","high","medium","low","very low","NA")), Fold=if_else(id == 1, Fold, as.numeric(NA))) %>%
  select(-id) %>%
  unite(Category, Species, Data, sep = " ", remove = TRUE) %>%
  mutate(Category = factor(Category, levels = unique(Category))) %>%
  mutate(Assembly = factor(Assembly, levels = unique(Assembly)))

text_size <- 20
tick_text_size <- 16

#minute <- 60
#hour <- 60*minute
#day <- 24*hour
#week <- 7*day
#time_ticks <- c(c(seq(2,18,2),seq(30,120,30))*minute,seq(1,6,1)*day,seq(1,41,2)*week)
#time_labels <- c(str_c(c(seq(2,18,2),seq(30,120,30)),"m"),str_c(seq(1,6,1),"d"),str_c(seq(1,9,2),"w"),if_else((1:16)%%2 == 0,str_c(seq(11,41,2),"w"), ''))

finish_plot <- list(
  scale_color_manual(values=c("#488BC2","#7FB972","#B5BD4C","#D92120","#E6642C","#781C81","#664CFF","#BBBBBB","#D9AD3C")),
  scale_shape_manual(values=c(17, 15, 19, 18, 25, 7)),
  #scale_x_continuous(breaks=time_ticks, labels=time_labels),
  facet_wrap(. ~ Category, scales="free", ncol=2),
  theme_bw(),
  theme(axis.text.x = element_text( size = tick_text_size),
        axis.text.y = element_text( size = tick_text_size),
        axis.title.x = element_text( size = text_size),
        axis.title.y = element_text( size = text_size),
        strip.text.x = element_text( size = tick_text_size),
        strip.text.y = element_text( size = tick_text_size),
        legend.title=element_text( size = text_size),
        legend.text=element_text(size=tick_text_size),
        legend.position = "top",
        legend.key.size = unit(8, 'mm'),
        legend.box = "vertical",
        legend.margin=margin(),
        panel.grid.minor = element_blank())
)

data %>%
  ggplot(aes(x=cpu_time, y=memory_gb, color=Assembly, shape=Coverage)) +
    geom_point(size=6) +
    #scale_x_continuous(limits = c(0, NA), expand=c(0.07,0.05)) +
    xlab("CPU time (min/day)") +
    ylab("Max. memory (Gb)") +
    finish_plot


ggsave(outfile1, width=297, height=210, units="mm")

data %>%
  ggplot(aes(x=elapsed_time, y=memory_gb, color=Assembly, shape=Coverage)) +
    geom_point(size=6, na.rm=TRUE) +
    xlab("Elapsed time (min/day)") +
    ylab("Max. memory (Gb)") +
    finish_plot

ggsave(outfile2, width=297, height=210, units="mm")

