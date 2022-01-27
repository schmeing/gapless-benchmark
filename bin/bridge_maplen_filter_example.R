suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(patchwork))

args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  stop("The example csv file must be supplied", call.=FALSE)
} else if(length(args)>1){
  stop("Only the example csv file must be supplied", call.=FALSE)
}
csv_file <- args[1]
#csv_file <- "/run/user/1000/gvfs/smb-share:server=imlsportmacquarie.uzh.ch,share=shared/schmeing/minigap_benchmark/input/csv/bridge_maplen_filter_example.csv"

full_csv <- read_csv(csv_file, col_types = cols())

data <- full_csv %>%
  select(from, from_side, to, to_side, mean_dist, from_extlen, to_gaplen, to_extlen) %>%
  rename(gaplen = to_gaplen) %>%
  filter(from == 1514 & from_side == "l")

groups <- data %>%
  group_by(from, from_side, to, to_side, mean_dist) %>%
  summarise(count=n()) %>%
  filter(count > 5) %>%
  ungroup() %>%
  mutate(group = 1, group = cumsum(group))

data <- data %>%
  inner_join(groups, by=c("from","from_side","to","to_side","mean_dist")) %>%
  select(-from, -from_side, -to, -to_side, -mean_dist) %>%
  arrange(desc(from_extlen), desc(gaplen), desc(to_extlen)) %>%
  mutate(from_pos=1, from_pos=cumsum(from_pos)) %>%
  arrange(desc(to_extlen), desc(gaplen), desc(from_extlen)) %>%
  mutate(to_pos=1, to_pos=cumsum(to_pos))

pdata <- data %>%
  mutate(group = if_else(group==1, "Test", "Control")) %>%
  rename(Group = group) %>%
  gather(type,Length,from_extlen,gaplen,to_extlen)

text_size <- 20
tick_text_size <- 16

p1 <- pdata %>%
  filter(type != "from_extlen") %>%
  ggplot(aes(x=-to_pos, y=-Length, color=Group, fill=Group)) +
    geom_col() +
    geom_col(aes(x=-to_pos, y=Length, color=Group, fill=Group), filter(pdata, type == "from_extlen")) +
    geom_col(aes(x=-to_pos, y=-Length), filter(pdata, type == "gaplen"), color="black", fill="black") +
    geom_hline(yintercept=0.0) +
    coord_flip() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_color_manual(values=c("#7FB972","#488BC2")) +
    scale_fill_manual(values=c("#7FB972","#488BC2")) +
    ggtitle("Nlonger = 4    p = 0.056") +
    ylab("Length") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size ),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text( size = text_size ),
          axis.title.y = element_blank(),
          plot.title = element_text( size=text_size, hjust=0.5 ),
          strip.text.x = element_text( size = text_size ),
          strip.text.y = element_text( size = text_size ),
          legend.title=element_text( size=text_size ),
          legend.text=element_text( size=tick_text_size ),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())

p2 <- pdata %>%
  filter(type != "to_extlen") %>%
  ggplot(aes(x=-from_pos, y=Length, color=Group, fill=Group)) +
    geom_col() +
    geom_col(aes(x=-from_pos, y=-Length, color=Group, fill=Group), filter(pdata, type == "to_extlen")) +
    geom_col(aes(x=-from_pos, y=Length), filter(pdata, type == "gaplen"), color="black", fill="black") +
    geom_hline(yintercept=0.0) +
    coord_flip() +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_color_manual(values=c("#7FB972","#488BC2")) +
    scale_fill_manual(values=c("#7FB972","#488BC2")) +
    ggtitle("Nlonger = 6    p = 0.009") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size ),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text( size = text_size ),
          axis.title.y = element_blank(),
          plot.title = element_text( size=text_size, hjust=0.5 ),
          strip.text.x = element_text( size = text_size ),
          strip.text.y = element_text( size = text_size ),
          legend.title=element_text( size=text_size ),
          legend.text=element_text( size=tick_text_size ),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())

p1 + p2 + plot_layout(guides = "collect") + plot_annotation(title = "Ncontrol = 12    Ntest = 11", theme = theme(plot.title = element_text( size = text_size, hjust=0.42 )))
