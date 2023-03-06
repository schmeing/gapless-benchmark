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
  read_tsv(quast_file, col_types = cols()) %>% select(Assembly,`# contigs`,`Largest contig`,`Total length`,`NGA50`,`# misassemblies`, `Genome fraction (%)`, `Duplication ratio`, NGA50, `Reference length`) %>% mutate(Folder=folder)
})
quast_tsv <- bind_rows(quast_tsv)

data <- quast_tsv %>%
  mutate(Type = if_else(substr(Assembly, nchar(Assembly)-6, nchar(Assembly))  == "_broken", "Contig", "Scaffold"), misassemblies=`# misassemblies`/`Reference length`*1e6) %>%
  select(Folder, Type, NGA50, misassemblies, "Genome fraction (%)", "Duplication ratio")

data <- assemblies %>%
  select(-Quast) %>%
  uncount(2, .id="Type") %>%
  mutate(Type=if_else(Type == 1, "Contig", "Scaffold")) %>%
  full_join(data, by=c("Folder","Type")) %>%
  group_by(Folder) %>%
  mutate(NGA50 = if_else(is.na(NGA50), min(NGA50, na.rm=TRUE), NGA50), misassemblies = if_else(is.na(misassemblies), min(misassemblies, na.rm=TRUE), misassemblies), `Genome fraction (%)` = if_else(is.na(`Genome fraction (%)`), min(`Genome fraction (%)`, na.rm=TRUE), `Genome fraction (%)`), `Duplication ratio` = if_else(is.na(`Duplication ratio`), min(`Duplication ratio`, na.rm=TRUE), `Duplication ratio`)) %>%
  ungroup %>%
  mutate(NGA50 = NGA50/1e6, Coverage = replace_na(Coverage, "other"), nrows=if_else((Assembly == "supernova" | (Assembly == "Flye" & Coverage == "high" & Data == "PacBio HiFi")) & Species == "Human", 2, 1)) %>%
  uncount(nrows, .id="id") %>%
  mutate(Data=if_else(Data != "Illumina" & Data != "10X Chromium", Data, if_else(Species != "Human", "PacBio CLR", if_else(id == 1, "PacBio HiFi", "Nanopore")))) %>%
  mutate(Data=if_else(id == 1, Data, "Nanopore"), Coverage=factor(if_else(id == 1, Coverage, "other"), levels=c("very high","high","medium","low","very low","other")), Fold=if_else(id == 1, Fold, as.numeric(NA))) %>%
  select(-id) %>%
  mutate(Assembly = factor(Assembly, levels = unique(Assembly))) %>%
  unite(Category, Species, Data, sep = " ", remove = TRUE) %>%
  mutate(CType = if_else(Type == "Contig", "(Contig)", "(Scaffold)")) %>%
  unite(Category2, Category, CType, sep = " ", remove = FALSE) %>%
  select(-CType) %>%
  mutate(Category = factor(Category, levels = unique(Category))) %>%
  arrange(Category) %>%
  mutate(Category2 = factor(Category2, levels = unique(Category2)))

text_size <- 20
tick_text_size <- 16
#"#488BC2"',"#4065B1"
finish_plot <- list(
  scale_color_manual(values=c("#7FB972","#B5BD4C","#D92120","#E6642C","#488BC2","#664CFF","#D9AD3C","#B15928","#781C81","#BBBBBB","black")),
  scale_shape_manual(values=c(24, 22, 21, 23, 25, 7)),
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
        legend.box = "vertical",
        legend.margin=margin(l=-8, t = -1.5, b=-2, unit='mm'),
        panel.grid.minor = element_blank())
)

data %>%
  #filter(Category == "E. Coli PacBio CLR" | Category == "Dolphin PacBio CLR") %>%
  #filter(Category == "Human PacBio HiFi" | Category == "Human Nanopore") %>%
  ggplot(aes(x=`misassemblies`, y=NGA50, color=Assembly)) +
    geom_point(aes(shape=Coverage), size=4) +
    geom_path(data=filter(data,Coverage!="other"), size=2, show.legend=FALSE) +
    #geom_path(data=filter(data,Coverage!="other", Category == "E. Coli PacBio CLR" | Category == "Dolphin PacBio CLR"), size=2, show.legend=FALSE) +
    #geom_path(data=filter(data,Coverage!="other", Category == "Human PacBio HiFi" | Category == "Human Nanopore"), size=2, show.legend=FALSE) +
    facet_wrap(Category2 ~ ., scales="free", ncol=2) +
    scale_x_continuous(limits = c(0, NA), expand=c(0.07,0.05)) +
    xlab("misassemblies (1/Mbp)") +
    ylab("NGA50 (Mbp)") +
    finish_plot +
    #scale_color_manual(values=c("#7FB972","#D92120","#488BC2","#D9AD3C","#781C81","#BBBBBB","black")) +
    guides(color=guide_legend(nrow=4, order=1))

ggsave(outfile1, width=210, height=240, units="mm")

data %>%
  filter(Type == "Scaffold") %>%
  ggplot(aes(x=`Duplication ratio`, y=`Genome fraction (%)`, color=Assembly)) +
    geom_point(aes(shape=Coverage), size=4) +
    geom_path(data=filter(data,Coverage!="other",Type == "Scaffold"), size=2, show.legend=FALSE) +
    facet_wrap(Category ~ ., scales="free", ncol=2) +
    xlab("Duplication ratio") +
    ylab("Genome fraction (%)") +
    finish_plot +
    guides(color=guide_legend(nrow=4, order=1))

ggsave(outfile2, width=297, height=210, units="mm")

