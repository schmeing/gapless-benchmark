suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)
if(length(args)<4){
  stop("The base directory, the input information and two output files must be supplied", call.=FALSE)
} else if(length(args)>4){
  stop("Only the base directory, the input information and two output files must be supplied", call.=FALSE)
}
base_dir <- args[1]
assembly_csv <- args[2]
outfile1 <- args[3]
outfile2 <- args[4]

assemblies <- read_csv(paste0(base_dir,"/",assembly_csv), col_types = cols(), na = c("", "NA", "-"))

quast_tsv <- lapply(1:length(assemblies[["Folder"]]), function(i){
  folder <- assemblies[["Folder"]][i]
  quast_file <- paste0(base_dir,"/",folder,"/",assemblies[["Quast"]][i],"/transposed_report.tsv")
  read_tsv(quast_file, col_types = cols()) %>% select(Assembly,`# contigs`,`Largest contig`,`Total length`,`NG50`,`# misassemblies`, `Genome fraction (%)`, `Duplication ratio`, NGA50, `Reference length`) %>% mutate(Folder=folder)
})
quast_tsv <- bind_rows(quast_tsv)

data <- quast_tsv %>%
  mutate(Type = if_else(substr(Assembly, nchar(Assembly)-6, nchar(Assembly))  == "_broken", "Contig", "Scaffold"), misassemblies=`# misassemblies`/`Reference length`*1e6) %>%
  select(Folder, Type, NG50, misassemblies, "Genome fraction (%)", "Duplication ratio")

data <- assemblies %>%
  select(-Quast) %>%
  uncount(2, .id="Type") %>%
  mutate(Type=if_else(Type == 1, "Contig", "Scaffold")) %>%
  full_join(data, by=c("Folder","Type")) %>%
  group_by(Folder) %>%
  mutate(NG50 = if_else(is.na(NG50), min(NG50, na.rm=TRUE), NG50), misassemblies = if_else(is.na(misassemblies), min(misassemblies, na.rm=TRUE), misassemblies), `Genome fraction (%)` = if_else(is.na(`Genome fraction (%)`), min(`Genome fraction (%)`, na.rm=TRUE), `Genome fraction (%)`), `Duplication ratio` = if_else(is.na(`Duplication ratio`), min(`Duplication ratio`, na.rm=TRUE), `Duplication ratio`)) %>%
  ungroup %>%
  mutate(NG50 = NG50/1e6, Coverage = replace_na(Coverage, "NA"), nrows=if_else((Assembly == "supernova" | (Assembly == "Flye" & Coverage == "high" & Data == "PacBio HiFi")) & Species == "Human", 2, 1)) %>%
  uncount(nrows, .id="id") %>%
  mutate(Data=if_else(Data != "Illumina" & Data != "10X Chromium", Data, if_else(Species != "Human", "PacBio CLR", if_else(id == 1, "PacBio HiFi", "Nanopore")))) %>%
  mutate(Data=if_else(id == 1, Data, "Nanopore"), Coverage=factor(if_else(id == 1, Coverage, "NA"), levels=c("very high","high","medium","low","very low","NA")), Fold=if_else(id == 1, Fold, as.numeric(NA))) %>%
  select(-id) %>%
  unite(Category, Species, Data, sep = " ", remove = TRUE) %>%
  mutate(Category = factor(Category, levels = unique(Category))) %>%
  mutate(Assembly = factor(Assembly, levels = unique(Assembly)))

text_size <- 20
tick_text_size <- 16

finish_plot <- list(
  scale_color_manual(values=c("#488BC2","#7FB972","#B5BD4C","#D92120","#E6642C","#781C81","#664CFF","#BBBBBB","#D9AD3C")),
  scale_shape_manual(values=c(17, 15, 19, 18, 25, 7)),
  facet_wrap(Type ~ Category, scales="free", ncol=4),
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
  ggplot(aes(x=`misassemblies`, y=NG50, color=Assembly, shape=Coverage)) +
    geom_point(size=6) +
    scale_x_continuous(limits = c(0, NA), expand=c(0.07,0.05)) +
    xlab("misassemblies / Mbp") +
    ylab("NG50 (Mbp)") +
    finish_plot


ggsave(outfile1, width=297, height=210, units="mm")

data %>%
  ggplot(aes(x=`Duplication ratio`, y=`Genome fraction (%)`, color=Assembly, shape=Coverage)) +
    geom_point(size=6) +
    xlab("Duplication ratio") +
    ylab("Genome fraction (%)") +
    finish_plot

ggsave(outfile2, width=297, height=210, units="mm")

