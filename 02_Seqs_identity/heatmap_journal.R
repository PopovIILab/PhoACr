# Set the dir

main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Load/install required packages

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(tidyverse, ggtext, patchwork, ggnewscale)

################################
# NUCL. SEQS. IDENTITY HEATMAP #
################################

# Datasets of nucleotide sequences identity were provided by our collaborators

# Read the data

data <- read.table('data.csv',
                   sep = ';',
                   comment = '',
                   head = T)

data_rows <- data %>%
  pivot_longer(cols = -GENBANK_ID,
               names_to = "OUR_ID",
               values_to = "identity")

# Read the data on our seqs

data_our <- read.table('data_our.csv',
                       sep = ';',
                       comment = '',
                       head = T)

data_our_rows <- data_our %>%
  pivot_longer(cols = -GENBANK_ID,
               names_to = "OUR_ID",
               values_to = "identity")

data_rows$identity <- data_rows$identity * 100

data_general <- merge(data_our_rows, data_rows, all = TRUE)

# Draft heatmap

p1 <- data_general %>%
  ggplot(aes(x = OUR_ID, fill = identity, y = GENBANK_ID)) +
  geom_tile() +
  geom_text(aes(label = identity), size = 3) +
  theme_classic() +
  scale_fill_gradient2(
    name = "Nucleotide<br>sequence<br>identity<br>",
    low = "white",
    mid = "lightblue",
    high = "darkblue",
    na.value = "white",
    expand = c(0, 0),
    limits = c(50, 100),
    midpoint = 85,
    position = "left",
    breaks = seq(50, 100, by = 5),
    labels = paste0(seq(50, 100, by = 5), "%")
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(
    labels = function(labels)
      gsub("^PK", "", labels),
    position = "left"
  ) +
  theme(
    legend.title = element_markdown(size = 14),
    legend.text = element_text(size = 14),
    legend.key.height = unit(80, "pt"),
    legend.position = "left",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x.top = element_text(angle = 45, hjust = 0)
  )

##########################
# ANNOTATING THE HEATMAP #
##########################

# Load `metadata` from the previous step

metadata <- read.table('../01_Phylogenetics/metadata/metadata.tsv',
                       sep = '\t',
                       header = T)
metadata <- rbind(metadata, c(Name = "PK", rep(NA, ncol(metadata) - 1)))

# Delete the row with "OQ725987.1". Do not ask why. Collaborators told us to.

metadata <- metadata[metadata$Name != "OQ725987.1", ]

# Create columns "X", "Y", "Z", "Location_title" & "Host_title" - will be needed for better visualization

metadata$X <- 'a'
metadata$Y <- 'b'
metadata$Z <- 'c'
metadata$Location_title <- NA
metadata$Location_title[metadata$Name == "AF353511.1"] <- "ND"
metadata$Location_title[metadata$Name == "AF304460.1"] <- "ND"
metadata$Location_title[metadata$Name == "AY567487.2"] <- "ND"
metadata$Host_title <- NA
metadata$Host_title[metadata$Name == "AF353511.1"] <- "ND"
metadata$Host_title[metadata$Name == "EF203064.1"] <- "ND"
metadata$Host_title[metadata$Name == "EU420139.1"] <- "ND"
metadata$Host_title[metadata$Name == "EU420138.1"] <- "ND"

# `Year` column

metadata$Year[metadata$Name == "DQ648858.1"] <- "2005"
metadata$Year[metadata$Name == "EF203064.1"] <- "2006"
metadata$Year[metadata$Name == "EU375869.1"] <- "2007"
metadata$Year[metadata$Name == "EU375864.1"] <- "2007"

# `Host` column

metadata$Host[metadata$Host == "ND"] <- NA
metadata$Host[grepl("^Myotis brandtii", metadata$Host)] <- "Myotis brandtii"
metadata$Host[metadata$Host == "unknown"] <- "Bat"
metadata$Host[metadata$Host == "bat"] <- "Bat"
metadata$Host[metadata$Host == "P.nathusii"] <- "Pipistrellus nathusii"
metadata$Host[metadata$Name == "EU375869.1"] <- "Pipistrellus nathusii"
metadata$Host[metadata$Name == "EU375864.1"] <- "Pipistrellus nathusii"
metadata$Host[metadata$Name == "DQ648858.1"] <- "Scotophilus kuhlii"
metadata$Host[metadata$Name == "OQ725983.1"] <- "Myotis brandtii"
metadata$Host[metadata$Name == "MG923572.2"] <- "Myotis brandtii"
metadata$Host[metadata$Name == "AF304460.1"] <- "Homo sapiens"
metadata$Host[metadata$Name == "AY567487.2"] <- "Homo sapiens"
metadata$Host[metadata$Name == "PQ450473.1"] <- "Myotis sibiricus"
metadata$Host[metadata$Name == "PQ450475.1"] <- "Myotis sibiricus"
metadata$Host[metadata$Name == "PQ450480.1"] <- "Myotis sibiricus"

# `Country` column

metadata$Country[metadata$Country == "ND"] <- NA
metadata$Country[grepl("^Russia: Nov", metadata$Country)] <- "Russia: Novosibirsk"
metadata$Country[metadata$Name == "OR052074.1"] <- "Russia: Rostov-on-Don"
metadata$Country[metadata$Name == "OR052075.1"] <- "Russia: Rostov-on-Don"
metadata$Country[metadata$Name == "OR052076.1"] <- "Russia: Rostov-on-Don"
metadata$Country[metadata$Name == "EF203064.1"] <- "China: Guangdong"
metadata$Country[metadata$Name == "EU375869.1"] <- "Germany: Bad Segeberg"
metadata$Country[metadata$Name == "EU375864.1"] <- "Germany: Bad Segeberg"
metadata$Country[metadata$Name == "DQ648858.1"] <- "China: Hainan"
names(metadata)[names(metadata) == "Country"] <- "Location"

# Add `Location` and `Host` heatmaps

p2 <- metadata %>%
  ggplot() +
  geom_tile(aes(x = Y, fill = Host, y = Name)) +
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
      "pig"
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
      "Pig*"
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
    na.value = "white"
  ) +
  geom_text(aes(x = Y, y = Name, label = Host_title), size = 3) +
  new_scale_fill() +
  geom_tile(aes(x = X, fill = Location, y = Name)) +
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
      "USA: Boulder County, CO"
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
      "USA: Boulder County"
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
    na.value = "white"
  ) +
  geom_text(aes(x = X, y = Name, label = Location_title), size = 3) +
  new_scale_fill() +
  geom_tile(aes(x = Z, fill = Year, y = Name), show.legend = F) + # Add `Year` information
  scale_fill_manual(
    values = c(
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white",
      "white"
    ),
    na.value = "white"
  ) +
  geom_text(aes(x = Z, y = Name, label = as.factor(Year)), size = 3) +
  scale_x_discrete(labels = c("Location", "Host", "Year"),
                   position = "top") +
  new_scale_fill() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x.top = element_text(
      angle = 45,
      hjust = 0,
      size = 12
    ),
    legend.text = element_markdown(size = 12),
    legend.title = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # Transparent plot background
    legend.background = element_rect(fill = "transparent", color = NA)  # Transparent legend background
  )

# Combine plots

combined <- p1 + p2 + plot_layout(design = "AAAAAAAABB")

# Save the fig

ggsave(
  "images/HM.png",
  plot = combined,
  width = 16,
  height = 12,
  dpi = 600
)