# Set working directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Install or call libraries
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(dplyr, stringr, tidyr, tidyverse, ggnewscale, ggtext, patchwork)


# Data management ---------------------------------------------------------

# Read the entire file as lines
file_path <- "results/selectivity.txt"
lines <- readLines(file_path)

# Convert lines into a data frame for easier handling
data <- data.frame(line = lines, stringsAsFactors = FALSE)

# Extract primer names
data <- data %>%
  mutate(Primer = ifelse(
    str_detect(line, "^Primer name"),
    str_extract(line, "(?<=Primer name ).+"),
    NA
  )) %>%
  fill(Primer, .direction = "down")

# Extract sequence names
data <- data %>%
  mutate(Sequence = ifelse(
    str_detect(line, "^\\s*Sequence:"),
    str_extract(line, "(?<=Sequence: )[^ ]+"),
    NA
  )) %>%
  fill(Sequence, .direction = "down")

# Extract forward and reverse mismatches
data <- data %>%
  mutate(
    Forward_mismatches = ifelse(
      str_detect(line, "hits forward strand"),
      as.numeric(str_extract(line, "\\d+(?= mismatches)")),
      NA
    ),
    Reverse_mismatches = ifelse(
      str_detect(line, "hits reverse strand"),
      as.numeric(str_extract(line, "\\d+(?= mismatches)")),
      NA
    )
  )

# Filter rows to keep only lines where mismatches are recorded
mismatch_data <- data %>%
  filter(!is.na(Forward_mismatches) |
           !is.na(Reverse_mismatches)) %>%
  select(Primer, Sequence, Forward_mismatches, Reverse_mismatches)

# Remove duplicate entries if they exist, by taking the first non-NA value in each column
mismatch_data <- mismatch_data %>%
  group_by(Primer, Sequence) %>%
  summarise(
    Forward_mismatches = first(na.omit(Forward_mismatches)),
    Reverse_mismatches = first(na.omit(Reverse_mismatches)),
    .groups = 'drop'
  )

# Pivot data so each primer has a Forward and Reverse column next to each other
mismatch_df <- mismatch_data %>%
  pivot_wider(
    names_from = Primer,
    values_from = c(Forward_mismatches, Reverse_mismatches),
    names_glue = "{Primer}_{.value}"
  )

# Reorder columns to group Forward and Reverse columns for each primer
ordered_columns <- c("Sequence", sort(colnames(mismatch_df)[-1], method = "radix"))
mismatch_df <- mismatch_df %>%
  select(all_of(ordered_columns))

# Define a function to convert numbers to text
number_to_text <- function(x) {
  if (is.na(x)) {
    return("No_hit")
  }
  
  # Define text representations for specific numbers
  words <- c("Zero", "One", "Two", "Three", "Four", "Five")
  
  # Convert numeric values to text if within range; otherwise, keep as numeric
  if (x >= 0 && x <= 5) {
    return(words[x + 1])
  } else {
    return(as.character(x))
  }
}

# Apply the conversion to each numeric column in mismatch_df
mismatch_df <- mismatch_df %>%
  mutate(across(-Sequence, ~ sapply(.x, number_to_text)))

# Display the final dataset
print(mismatch_df)

#Save the file
write.csv(mismatch_df, file = "mismatch_df.csv", row.names = FALSE)



# Heatmap_metadata -----------------------------------------------------------------

# Metadata was created manually by Igor Popov and consisted of `Sequence_id;Name;Host;Genera` columns

metadata <- read.table("virus_metadata.csv", sep = ";", header = T)

p1 <- metadata %>%
  ggplot() +
  geom_tile(aes(
    x = Axis_A,
    fill = Sequence_id,
    y = factor(
      Name,
      level = c(
        "Bat coronavirus",
        "Shrew coronavirus",
        "Rodent coronavirus",
        "Bat coronavirus isolate PREDICT/PDF-2180",
        "Duck coronavirus",
        "Infectious bronchitis virus isolate Ind-TN92-03",
        "Canada goose coronavirus",
        "Turkey coronavirus",
        "Beluga whale coronavirus SW1",
        "Avian infectious bronchitis virus",
        "Porcine coronavirus HKU15",
        "Common moorhen coronavirus HKU21",
        "Wigeon coronavirus HKU20",
        "Night heron coronavirus HKU19",
        "Magpie-robin coronavirus HKU18",
        "Sparrow coronavirus HKU17",
        "White-eye coronavirus HKU16",
        "Munia coronavirus HKU13-3514",
        "Thrush coronavirus HKU12-600",
        "Bulbul coronavirus HKU11-934",
        "Severe acute respiratory syndrome coronavirus 2",
        "Betacoronavirus England 1",
        "Rousettus bat coronavirus",
        "Middle East respiratory syndrome-related coronavirus",
        "Rabbit coronavirus HKU14",
        "Bat coronavirus BM48-31/BGR/2008",
        "Rousettus bat coronavirus HKU9",
        "Pipistrellus bat coronavirus HKU5",
        "Tylonycteris bat coronavirus HKU4",
        "SARS coronavirus Tor2",
        "Human coronavirus HKU1",
        "Murine hepatitis virus strain JHM",
        "Murine hepatitis virus strain A59",
        "Betacoronavirus Erinaceus/VMC/DEU/2012",
        "Betacoronavirus HKU24",
        "Bat Hp-betacoronavirus/Zhejiang2013",
        "Rat coronavirus Parker",
        "Human coronavirus OC43",
        "Bovine coronavirus",
        "Murine hepatitis virus strain MHV-A59",
        "Bat alphacoronavirus isolate AMA_L_F",
        "NL63-related bat coronavirus strain BtKYNL63-9b",
        "Wencheng Sm shrew coronavirus isolate Xingguo-74",
        "Alphacoronavirus Bat-CoV/P.kuhlii/Italy/3398-19/2015",
        "Wencheng Sm shrew coronavirus isolate Xingguo-101",
        "Coronavirus AcCoV-JC34",
        "Ferret coronavirus",
        "BtNv-AlphaCoV/SC2013",
        "BtRf-AlphaCoV/YN2012",
        "Bat coronavirus CDPHE15/USA/2006",
        "Human coronavirus NL63",
        "Lucheng Rn rat coronavirus",
        "Camel alphacoronavirus",
        "Rousettus bat coronavirus HKU10",
        "Human coronavirus 229E",
        "Alphacoronavirus sp. isolate WA3607",
        "Alphacoronavirus sp. isolate WA2028",
        "Bat alphacoronavirus isolate BtCoV/020_16/M.dau/FIN/2016",
        "Tylonycteris bat coronavirus HKU33",
        "Hipposideros pomona bat coronavirus CHB25",
        "Alphacoronavirus sp. isolate WA1087",
        "Transmissible gastroenteritis virus",
        "NL63-related bat coronavirus strain BtKYNL63-9a",
        "BtRf-AlphaCoV/HuB2013",
        "BtMr-AlphaCoV/SAX2011",
        "Swine enteric coronavirus",
        "Mink coronavirus strain WD1127",
        "Miniopterus bat coronavirus HKU8",
        "Bat coronavirus 1A",
        "Rhinolophus bat coronavirus HKU2",
        "Scotophilus bat coronavirus 512",
        "Porcine epidemic diarrhea virus",
        "Feline infectious peritonitis virus"
      )
    )
  ), show.legend = FALSE) +
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
    )
  ) +
  geom_text(aes(x = Axis_A, y = Name, label = Sequence_id), size = 3) +
  new_scale_fill() +
  geom_tile(aes(x = Axis_C, fill = Host, y = Name)) +
  scale_fill_manual(
    breaks = c(
      "Humans",
      "Bats",
      "Birds",
      "Rodents",
      "Pigs",
      "Shrews",
      "Cattle",
      "Camels",
      "Hedgehogs",
      "Ferrets",
      "Minks",
      "Cats",
      "Belugas",
      "Rabbits"
    ),
    values = c(
      "grey",
      "darkred",
      "forestgreen",
      "#8B814C",
      "#EEA2AD",
      "#FFFACD",
      "#CD9B1D",
      "#FFD39B",
      "#1874CD",
      "#00CDCD",
      "#EEDC82",
      "#FFFAFA",
      "#EED2EE",
      "#FFA07A"
    )
  ) +
  new_scale_fill() +
  geom_tile(aes(x = Axis_B, fill = Genera, y = Name)) +
  scale_fill_manual(
    breaks = c(
      "Alphacoronavirus",
      "Betacoronavirus",
      "Deltacoronavirus",
      "Gammacoronavirus"
    ),
    values = c("#FFF8DC", "#66CD00", "#1874CD", "#FFAEB9"),
    na.value = "white"
  ) +
  geom_text(aes(x = Axis_B, y = Name, label = Genera_title), size = 3) +
  scale_x_discrete(labels = c("a" = "GenBank ID", "b" = "CoV<br>genus", "c" = "Host")) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_markdown(
      angle = 45,
      hjust = 1,
      size = 16
    ),
    legend.text = element_markdown(size = 16),
    legend.title = element_text(size = 18),
    legend.position = "left",
    panel.background = element_blank()
  )

ggsave(
  "HM.png",
  plot = p1,
  width = 10,
  height = 12,
  dpi = 600
)

# Heatmap_PCR -------------------------------------------------------------

PCR <- read.table("mismatch_df.csv", sep = ",", header = T) %>%
  pivot_longer(cols = -c(Sequence),
               names_to = "Primers",
               values_to = "Mismatch")

p2 <- PCR %>%
  ggplot() +
  geom_tile(aes(
    x = Primers,
    fill = Mismatch,
    y = factor(
      Sequence,
      level = c(
        "NC_048212.1",
        "NC_046955.1",
        "NC_046954.1",
        "NC_034440.1",
        "NC_048214.1",
        "NC_048213.1",
        "NC_046965.1",
        "NC_010800.1",
        "NC_010646.1",
        "NC_001451.1",
        "NC_039208.1",
        "NC_016996.1",
        "NC_016995.1",
        "NC_016994.1",
        "NC_016993.1",
        "NC_016992.1",
        "NC_016991.1",
        "NC_011550.1",
        "NC_011549.1",
        "NC_011547.1",
        "NC_045512.2",
        "NC_038294.1",
        "NC_030886.1",
        "NC_019843.3",
        "NC_017083.1",
        "NC_014470.1",
        "NC_009021.1",
        "NC_009020.1",
        "NC_009019.1",
        "NC_004718.3",
        "NC_006577.2",
        "AC_000192.1",
        "NC_048217.1",
        "NC_039207.1",
        "NC_026011.1",
        "NC_025217.1",
        "NC_012936.1",
        "NC_006213.1",
        "NC_003045.1",
        "NC_001846.1",
        "NC_055953.1",
        "NC_048216.1",
        "NC_048211.1",
        "NC_046964.1",
        "NC_035191.1",
        "NC_034972.1",
        "NC_030292.1",
        "NC_028833.1",
        "NC_028824.1",
        "NC_022103.1",
        "NC_005831.2",
        "NC_032730.1",
        "NC_028752.1",
        "NC_018871.1",
        "NC_002645.1",
        "NC_076685.1",
        "NC_076684.1",
        "NC_076629.1",
        "NC_054015.1",
        "NC_054004.1",
        "NC_054003.1",
        "NC_038861.1",
        "NC_032107.1",
        "NC_028814.1",
        "NC_028811.1",
        "NC_028806.1",
        "NC_023760.1",
        "NC_010438.1",
        "NC_010437.1",
        "NC_009988.1",
        "NC_009657.1",
        "NC_003436.1",
        "NC_002306.3"
      )
    )
  )) +
  scale_fill_manual(
    name = "Mismatches",
    breaks = c("Four", "No_hit", "Two", "Three", "One", "Zero"),
    labels = c("No hit", "No hit", "2", "No hit", "1", "0"),
    values = c(
      "#CD5555",
      "#CD5555",
      "#CAFF70",
      "#CD5555",
      "#A2CD5A",
      "#458B00"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "aCovbat1_Forward_mismatches" = "PrimA F",
      "aCovbat1_Reverse_mismatches" = "PrimA R",
      "aCovbat2_Forward_mismatches" = "PrimB F",
      "aCovbat2_Reverse_mismatches" = "PrimB R",
      "bVijgen_Forward_mismatches" = "PrimC F",
      "bVijgen_Reverse_mismatches" = "PrimC R",
      "cWatanabe_Forward_mismatches" = "PrimD F",
      "cWatanabe_Reverse_mismatches" = "PrimD R",
      "dHolbrook1_Forward_mismatches" = "PrimE F",
      "dHolbrook1_Reverse_mismatches" = "PrimE R",
      "dHolbrook1_Forward_mismatches" = "PrimF F",
      "dHolbrook1_Reverse_mismatches" = "PrimF R",
      "dHolbrook2_Forward_mismatches" = "PrimG F",
      "dHolbrook2_Reverse_mismatches" = "PrimG R",
      "dHolbrook3_Forward_mismatches" = "PrimH F",
      "dHolbrook3_Reverse_mismatches" = "PrimH R"
    )
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_markdown(
      angle = 55,
      hjust = 1,
      size = 14
    ),
    axis.text.y = element_blank(),
    legend.text = element_markdown(size = 16),
    legend.title = element_text(size = 18),
    legend.position = "right",
    panel.background = element_blank()
  )

# Combine plots

combined <- p1 + p2 + plot_layout(design = "AABB")

# Save the fig

ggsave(
  "images/HM.png",
  plot = combined,
  width = 14,
  height = 12,
  dpi = 600
)
