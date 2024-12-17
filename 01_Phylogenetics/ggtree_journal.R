# Set the dir

main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Load/install required packages

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  ggtree,
  ggimage,
  ggtext,
  phangorn,
  ggplot2,
  treeio,
  ggnewscale,
  viridis,
  phytools,
  patchwork
)

####################
# METADATA PARSING #
####################

# Load `metadata`

metadata <- read.table('metadata/metadata.tsv', sep = '\t', header = T)

# `Year` column

meta.year <- as.data.frame(metadata[, 'Year'])
colnames(meta.year) <- 'Year'
rownames(meta.year) <- metadata$Name
meta.year$Year[meta.year$Year == "ND"] <- NA
meta.year["DQ648858.1", "Year"] <- "2005"
meta.year["EF203064.1", "Year"] <- "2006"
meta.year["EU375869.1", "Year"] <- "2007"
meta.year["EU375864.1", "Year"] <- "2007"

# `Host` column

meta.host <- as.data.frame(metadata[, 'Host'])
colnames(meta.host) <- 'Host'
rownames(meta.host) <- metadata$Name
meta.host$Host[meta.host$Host == "ND"] <- NA
meta.host$Host[grepl("^Myotis brandtii", meta.host$Host)] <- "Myotis brandtii"
meta.host$Host[meta.host$Host == "unknown"] <- "Bat"
meta.host$Host[meta.host$Host == "bat"] <- "Bat"
meta.host$Host[meta.host$Host == "P.nathusii"] <- "Pipistrellus nathusii"
meta.host["EU375869.1", "Host"] <- "Pipistrellus nathusii"
meta.host["EU375864.1", "Host"] <- "Pipistrellus nathusii"
meta.host["DQ648858.1", "Host"] <- "Scotophilus kuhlii"
meta.host["OQ725983.1", "Host"] <- "Myotis brandtii"
meta.host["MG923572.2", "Host"] <- "Myotis brandtii"
meta.host["AF304460.1", "Host"] <- "Homo sapiens"
meta.host["AY567487.2", "Host"] <- "Homo sapiens"
meta.host["PQ450473.1", "Host"] <- "Myotis sibiricus"
meta.host["PQ450475.1", "Host"] <- "Myotis sibiricus"
meta.host["PQ450480.1", "Host"] <- "Myotis sibiricus"

# `Country` column

meta.loc <- as.data.frame(metadata[, 'Country'])
colnames(meta.loc) <- 'Country'
rownames(meta.loc) <- metadata$Name
meta.loc$Country[meta.loc$Country == "ND"] <- NA
meta.loc$Country[grepl("^Russia: Nov", meta.loc$Country)] <- "Russia: Novosibirsk"
meta.loc["OR052074.1", "Country"] <- "Russia: Rostov-on-Don"
meta.loc["OR052075.1", "Country"] <- "Russia: Rostov-on-Don"
meta.loc["OR052076.1", "Country"] <- "Russia: Rostov-on-Don"
meta.loc["EF203064.1", "Country"] <- "China: Guangdong"
meta.loc["EU375869.1", "Country"] <- "Germany: Bad Segeberg"
meta.loc["EU375864.1", "Country"] <- "Germany: Bad Segeberg"
meta.loc["DQ648858.1", "Country"] <- "China: Hainan"
colnames(meta.loc) <- c("Location")


######################
# TREE VISUALIZATION #
######################

# Read the tree file

tree <- read.tree("tree/tree_ufb.treefile")

# Midpoint root the tree

midpoint.root(tree)

# Draft tree

tree_fig <- ggtree(tree) %<+% metadata +
  xlim(0, 4.7) +
  geom_tiplab(
    aes(
      label = Full.Name,
      #fontface = ifelse(grepl("^PQ", Full.Name), "bold", "plain"),
      # Bold if starts with "PQ"
      color = ifelse(
        grepl(paste0("^(", paste(
          c(
            "PQ439331.1",
            "PQ439332.1",
            "PQ439333.1",
            "PQ439334.1",
            "PQ450477.1",
            "PQ450472.1",
            "PQ450475.1",
            "PQ450471.1",
            "PQ450476.1",
            "PQ450479.1",
            "PQ450474.1",
            "PQ450473.1",
            "PQ450480.1"
          ),
          collapse = "|"
        ), ")"), Full.Name),
        "#d335f2",
        ifelse(
          grepl(paste0("^(", paste(
            c("DQ648858.1", "MG923574.2"), collapse = "|"
          ), ")"), Full.Name),
          "#6f1c80",
          ifelse(
            grepl(paste0("^(", paste(
              c("PQ450482.1", "PQ450478.1", "PQ450481.1", "PQ450483.1"),
              collapse = "|"
            ), ")"), Full.Name),
            "#35b3f2",
            ifelse(
              grepl("^KJ473806.1", Full.Name),
              "#1c5e80",
              ifelse(
                grepl("^PQ439335.1", Full.Name),
                "#f2a035",
                ifelse(grepl(paste0(
                  "^(", paste(c("MK472068.1", "MK720944.1"), collapse = "|"), ")"
                ), Full.Name), "#80541c", "black")
              )
            )
          )
        )
      )
    ),
    align = TRUE,
    geom = "label",
    fill = "white",
    label.size = 0
  ) +
  scale_color_identity()

# Onehot encode bootstrap values (<70 = 0; >70 = 1)

tree_boot <- tree_fig$data
tree_boot <- tree_boot[!tree_boot$isTip, ]
tree_boot$label <- as.numeric(tree_boot$label)
tree_boot$bootstrap <- '0'
tree_boot$bootstrap[tree_boot$label >= 70] <- '1'
tree_boot$bootstrap[is.na(tree_boot$label)] <- '1'

# Add bootstrap values to the tree (black branches = bootstrap >70; grey branches = bootstrap <70)

tree_fig <- tree_fig + new_scale_color() +
  geom_tree(data = tree_boot, aes(color = bootstrap == '1')) +
  scale_color_manual(name = 'Bootstrap',
                     values = setNames(c("black", "grey"), c(T, F)),
                     guide = "none")

# Add `Location` heatmap

tree_fig <- gheatmap(
  tree_fig,
  meta.loc,
  width = 0.05,
  offset = 1.5,
  color = "black",
  font.size = 4,
  colnames_offset_y = -0.05
) +
  scale_fill_manual(
    breaks = c(
      "Russia",
      "Russia: Rostov-on-Don",
      "Russia: Novosibirsk",
      "Finland: Espoo",
      "Finland: Mustasaari",
      "Denmark: Vadum",
      "Denmark: Moensted",
      "Germany: Bad Segeberg",
      "Italy",
      "China: Guangdong",
      "China: Hainan",
      "Hong Kong",
      "China",
      "Australia",
      "Kenya",
      "Peru",
      "USA: Wisconsin",
      "USA:Indiana",
      "USA: Boulder County, CO",
      NA
    ),
    labels = c(
      "Russia: Moscow",
      "Russia: Rostov-on-Don",
      "Russia: Novosibirsk",
      "Finland: Espoo",
      "Finland: Mustasaari",
      "Denmark: Vadum",
      "Denmark: Moensted",
      "Germany: Bad Segeberg",
      "Italy",
      "China: Guangdong",
      "China: Hainan",
      "China: Hong-Kong",
      "China",
      "Australia",
      "Kenya",
      "Peru",
      "USA: Wisconsin",
      "USA: Indiana",
      "USA: Boulder County",
      "NA"
    ),
    values = c(
      "#8B2323",
      "#FF4040",
      "#FF7256",
      "#FFF8DC",
      "#CDC8B1",
      "#C1FFC1",
      "#9BCD9B",
      "#EEC900",
      "#698B22",
      "#FFB5C5",
      "#EEA9B8",
      "#CD919E",
      "#8B636C",
      "#473C8B",
      "#9FB6CD",
      "#FF7F24",
      "#7EC0EE",
      "#6CA6CD",
      "#4A708B"
    ),
    na.translate = TRUE,
    na.value = "white",
    name = "Location"
  ) +
  new_scale_fill()

# Add `Host` heatmap

tree_fig <- gheatmap(
  tree_fig,
  meta.host,
  width = 0.05,
  offset = 1.7,
  color = "black",
  font.size = 4,
  colnames_offset_y = -0.05
) +
  scale_fill_manual(
    breaks = c(
      "Homo sapiens",
      "Nyctalus noctula",
      "Nyctalus velutinus",
      "Pipistrellus kuhlii",
      "Pipistrellus nathusii",
      "Pipistrellus pygmaeus",
      "Myotis brandtii",
      "Myotis dasycneme",
      "Myotis daubentonii",
      "Myotis lucifugus",
      "Myotis petax",
      "Myotis ricketti",
      "Myotis sibiricus",
      "Vespertilio murinus",
      "Scotophilus kuhlii",
      "Tylonycteris robustula",
      "Rhinolophus ferrumequinum",
      "Triaenops afer",
      "Hipposideros larvatus",
      "Desmodus rotundus (common vampire bat)",
      "Bat",
      "microbat",
      "Rattus norvegicus",
      #Грызун
      "Mustela vison",
      #Норка
      "Sorex araneus",
      #Землеройка
      "Suncus murinus",
      #Землеройка
      "pig",
      NA
    ),
    labels = c(
      "*Homo sapiens*",
      "*Nyctalus noctula*",
      "*Nyctalus velutinus*",
      "*Pipistrellus kuhlii*",
      "*Pipistrellus nathusii*",
      "*Pipistrellus pygmaeus*",
      "*Myotis brandtii*",
      "*Myotis dasycneme*",
      "*Myotis daubentonii*",
      "*Myotis lucifugus*",
      "*Myotis petax*",
      "*Myotis ricketti*",
      "*Myotis sibiricus*",
      "*Vespertilio murinus*",
      "*Scotophilus kuhlii*",
      "*Tylonycteris robustula*",
      "*Rhinolophus ferrumequinum*",
      "*Triaenops afer*",
      "*Hipposideros larvatus*",
      "*Desmodus rotundus*",
      "Bat*",
      "Microbat*",
      "*Rattus norvegicus*",
      #Грызун
      "*Mustela vison*",
      #Норка
      "*Sorex araneus*",
      #Землеройка
      "*Suncus murinus*",
      #Землеройка
      "Pig*",
      "NA"
    ),
    values = c(
      "#D2B48C",
      # "Homo_sapiens"
      "#00FF7F",
      # "Nyctalus_noctula"
      "#008B45",
      # "Nyctalus_velutinus"
      "#FFFF00",
      # "Pipistrellus_kuhlii"
      "#EEEE00",
      # "Pipistrellus_nathusii"
      "#8B8B00",
      # "Pipistrellus_pygmaeus"
      "#8B475D",
      # "Myotis_brandtii"
      "#8B636C",
      # "Myotis_dasycneme"
      "#CD6889",
      # "Myotis_daubentonii"
      "#CD919E",
      # "Myotis_lucifugus"
      "#CD96CD",
      # "Myotis_petax"
      "#FFB5C5",
      # "Myotis_ricketti"
      "#FFBBFF",
      # "Myotis_sibiricus"
      "#8B2252",
      # "Vespertilio_murinus"
      "#D02090",
      # "Scotophilus_kuhlii"
      "#5D478B",
      # "Tylonycteris_robustula"
      "#36648B",
      # "Rhinolophus_ferrumequinum"
      "#4F94CD",
      # "Triaenops_afer"
      "#5CACEE",
      # "Hipposideros_larvatus"
      "#63B8FF",
      # "Desmodus_rotundus"
      "#838B83",
      # "Bat"
      "#C1CDC1",
      # "Microbat"
      "#00868B",
      # "Rattus_norvegicus"
      "#668B8B",
      # "Mustela_vison"
      "#FFC125",
      # "*Sorex araneus*"
      "#8B4500",
      # "*Suncus murinus*"
      "#EEE9E9"
    ),
    # "Pig"
    na.translate = TRUE,
    na.value = "white",
    name = "Host"
  ) +
  guides(fill = guide_legend(title.theme = element_markdown(), label.theme = element_markdown())) +
  theme(legend.title = element_text(), legend.text = element_text()) +
  new_scale_fill()

# Add `Year` heatmap

tree_fig <- gheatmap(
  tree_fig,
  meta.year,
  width = 0.05,
  offset = 1.9,
  color = "black",
  font.size = 4,
  colnames_offset_y = -0.05
) +
  scale_fill_viridis(
    option = "D",
    name = "Year",
    discrete = TRUE,
    na.translate = TRUE
  )

# Save the fig

ggsave(
  'images/tree.png',
  tree_fig,
  width = 20,
  height = 16,
  dpi = 600
)
