library(tidyverse)
library(here)
library(cooccur)
knitr::opts_chunk$set(dpi = 300, dev = "png")
options(list(dplyr.summarise.inform = FALSE, future.rng.onMisuse = "ignore"))
daphnia <- readxl::read_excel(
  here("Data/Empirical/Parasite_data.xlsx"),
  col_types = c(
    "numeric",
    "date",
    "numeric",
    "skip",
    "skip",
    "date",
    "skip",
    "skip",
    "text",
    "text",
    "numeric",
    "skip",
    "skip",
    "text",
    "numeric",
    "numeric",
    "numeric",
    "numeric",
    "numeric",
    "numeric",
    "numeric",
    "numeric",
    "text",
    "numeric",
    "text",
    "numeric",
    "text",
    "text",
    "numeric",
    "text",
    "numeric",
    "text",
    "numeric",
    "text",
    "numeric",
    "text",
    "text",
    "text",
    "text",
    "skip",
    "text",
    "text",
    "skip",
    "text",
    "text",
    "text",
    "text"
  )
)
snail <- read_csv(
  here("Data/Empirical/Trematode Infections.csv"),
  show_col_types = FALSE
)
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
trematodes <- snail |>
  select(starts_with("trematode")) |>
  flatten_chr() |>
  unique() |>
  discard(
    .x = _,
    ~ is.na(.x)
  ) |>
  discard(
    .x = _,
    ~ str_detect(.x, "\\?")
  ) |>
  discard(
    .x = _,
    ~ str_detect(.x, "\\*")
  ) |>
  discard(
    .x = _,
    ~ str_detect(.x, "\\.")
  ) |>
  discard(
    .x = _,
    ~ str_detect(.x, "^im.+")
  ) |>
  discard(.x = _, ~ .x == "u")
snail_data_prep <- function(df) {
  all_snails <- df |>
    rownames_to_column("ID") |>
    select(ID) |>
    distinct()

  processed_data <- df |>
    rownames_to_column("ID") |>
    pivot_longer(
      cols = starts_with("trematode"),
      names_to = "trematode_survey",
      values_to = "trematode"
    ) |>
    filter(trematode %in% trematodes) |>
    select(ID, trematode) |>
    distinct() |>
    mutate(presence = TRUE)

  all_snails |>
    cross_join(tibble(trematode = unique(processed_data$trematode))) |>
    left_join(processed_data, by = c("ID", "trematode")) |>
    mutate(presence = coalesce(presence, FALSE)) |>
    pivot_wider(
      names_from = trematode,
      values_from = presence,
      id_cols = ID
    ) |>
    column_to_rownames(var = "ID") |>
    as.matrix() |>
    t()
}
cooccurrence_matrices <- snail |>
  group_by(site) |>
  group_split() |>
  map(snail_data_prep)
cooccurrence_analyses <- cooccurrence_matrices |>
  map(\(m) quiet(cooccur(m, thresh = F, spp_names = TRUE, prob = "comb")))
competitive_pairs <- cooccurrence_analyses |>
  map("results") |>
  map(\(df) {
    df |>
      filter(exp_cooccur >= 1, p_lt < 0.05) |>
      select(sp1_name, sp2_name, obs_cooccur, exp_cooccur, p_lt)
  })
facilitative_pairs <- cooccurrence_analyses |>
  map("results") |>
  map(\(df) {
    df |>
      filter(exp_cooccur >= 1, p_gt < 0.05) |>
      select(sp1_name, sp2_name, obs_cooccur, exp_cooccur, p_gt)
  })
trematode_stats <- snail |>
  pivot_longer(
    cols = starts_with("trematode"),
    names_to = "trematode_survey",
    values_to = "trematode"
  ) |>
  group_by(`Snail no.`, site) |>
  reframe(infection_intensity = 1 - str_count(trematode, "^u$")) |>
  filter(!is.na(infection_intensity)) |>
  group_by(site) |>
  summarize(sum_infection_intensity = sum(infection_intensity))
df_cf <- snail |>
  pivot_longer(
    cols = starts_with("trematode"),
    names_to = "trematode_survey",
    values_to = "trematode"
  ) |>
  filter(trematode %in% trematodes) |>
  group_by(site) |>
  summarize(true_species_richness = n_distinct(trematode)) |>
  ungroup() |>
  mutate(
    facilitative_pairs = map_int(facilitative_pairs, nrow),
    competitive_pairs = map_int(competitive_pairs, nrow)
  ) |>
  left_join(trematode_stats, by = join_by(site))
p1 <- df_cf |>
  ggplot(aes(x = site)) +
  geom_col(
    aes(y = facilitative_pairs),
    fill = "#007BC3FF",
    alpha = 0.7,
    width = 0.6
  ) +
  geom_col(
    aes(y = -1 * competitive_pairs),
    fill = "#EF7C12FF",
    alpha = 0.7,
    width = 0.6
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  labs(
    title = "A. Species Associations in C. californica",
    x = "Site",
    y = "Count (Facilitative above, Competitive below)"
  ) +
  theme_minimal() +
  geom_text(
    aes(
      label = paste0("Species:\n", true_species_richness),
      y = facilitative_pairs + 1.5
    ),
    position = position_dodge(0.9)
  )

daphnia_data_prep <- function(df) {
  df |>
    pivot_longer(
      cols = contains("status"),
      names_to = "Parasite",
      values_to = "presence"
    ) |>
    select(id = `data order`, Parasite, presence) |>
    distinct() |>
    mutate(
      presence = presence == "u",
      Parasite = str_split_i(Parasite, " ", 1)
    ) |>
    pivot_wider(
      names_from = Parasite, # The species names become columns
      values_from = presence, # TRUE/FALSE based on presence
      values_fill = FALSE, # Fill absent species with FALSE
      id_cols = id,
      id_expand = TRUE
    ) |>
    column_to_rownames(var = "id") |>
    as.matrix() |>
    t()
}
cooccurrence_matrices <- daphnia |>
  group_by(lake, `date collected`) |>
  group_split() |>
  map(daphnia_data_prep)
cooccurrence_analyses <- cooccurrence_matrices |>
  map(\(m) quiet(cooccur(m, thresh = T, spp_names = TRUE, prob = "comb")))
p2 <- daphnia |>
  group_by(lake, `date collected`) |>
  group_keys() |>
  mutate(
    positive = map_int(cooccurrence_analyses, "positive"),
    negative = map_int(cooccurrence_analyses, "negative")
  ) |>
  pivot_longer(
    cols = c(positive, negative),
    names_to = "Association",
    values_to = "count"
  ) |>
  ggplot() +
  geom_histogram(
    aes(x = count, fill = Association),
    bins = 2,
    position = "dodge"
  ) +
  xlab("Number of Associations per Site/Month") +
  ggtitle(
    "B. Species Associations in Daphnia (blue = positive, orange = negative)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#EF7C12FF", "#007BC3FF")) +
  scale_x_continuous(breaks = c(0, 1))
ggsave("Figures/Fig4a_snail.png", p1, width = 8, height = 4)
ggsave("Figures/Fig4b_daphnia.png", p2, width = 8, height = 4)
