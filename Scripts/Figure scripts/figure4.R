library(tidyverse)
library(here)
library(cooccur)
library(patchwork)
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
df_cf <- data.frame(
  facilitative = map_int(facilitative_pairs, nrow),
  competitive = -map_int(competitive_pairs, nrow)
) |>
  bind_cols(
    snail |>
      group_by(site) |>
      summarize(sample_size = n())
  ) |>
  mutate(
    site = factor(
      site,
      levels = c(
        "c2",
        "c3",
        "c6",
        "c7",
        "c8",
        "f1",
        "f2",
        "f3",
        "f4",
        "f5",
        "f1 (original)",
        "f5 (original)"
      ),
      labels = c(
        "c2",
        "c3",
        "c6",
        "c7",
        "c8",
        "f1",
        "f2",
        "f3",
        "f4",
        "f5",
        "f1\n(original)",
        "f5\n(original)"
      )
    )
  )
p1 <- df_cf |>
  ggplot(aes(x = site)) +
  geom_col(
    aes(y = facilitative),
    fill = "#007BC3FF",
    alpha = 0.7,
    width = 0.6
  ) +
  geom_col(
    aes(y = competitive),
    fill = "#EF7C12FF",
    alpha = 0.7,
    width = 0.6
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  ggtitle(
    bquote(bold("A.") ~ "Species Associations in" ~ italic("C. californica"))
  ) +
  labs(
    x = "Site",
    y = "Count (Facilitative above, Competitive below)"
  ) +
  theme_minimal() +
  geom_text(
    aes(
      label = paste0("n=", sample_size),
      y = facilitative + 1.5
    ),
    position = position_dodge(0.9),
    vjust = 0,
    size = 3
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
      presence = presence != "u",
      Parasite = str_split_i(Parasite, " ", 1)
    ) |>
    filter(
      Parasite != "UGP",
      Parasite != "GPB"
    ) |>
    pivot_wider(
      names_from = Parasite, # The species names become columns
      values_from = presence, # TRUE/FALSE based on presence
      values_fill = FALSE, # Fill absent species with FALSE
      id_cols = id,
      id_expand = TRUE
    ) |>
    mutate(
      `B. paedophthorum` = BP | BC,
      `Spider_oomycete` = POD | SM
    ) |>
    select(-BP, -BC, -POD, -SM) |>
    column_to_rownames(var = "id") |>
    as.matrix() |>
    t()
}
cooccurrence_matrices <- daphnia |>
  group_by(lake, `date collected`) |>
  group_split() |>
  map(daphnia_data_prep)
cooccurrence_analyses <- cooccurrence_matrices |>
  map(possibly(
    \(m) quiet(cooccur(m, thresh = T, spp_names = TRUE, prob = "comb")),
    otherwise = NULL
  ))
p2 <- daphnia |>
  group_by(lake, `date collected`) |>
  group_keys() |>
  mutate(
    idx = seq_len(n()),
    positive = map_int(
      idx,
      ~ if_else(
        !is.null(cooccurrence_analyses[[.x]]),
        nrow(cooccurrence_analyses[[.x]]$results[
          cooccurrence_analyses[[.x]]$results$p_gt < 0.05,
        ]) %||%
          0L,
        0L
      )
    ),
    negative = map_int(
      idx,
      ~ if_else(
        !is.null(cooccurrence_analyses[[.x]]),
        nrow(cooccurrence_analyses[[.x]]$results[
          cooccurrence_analyses[[.x]]$results$p_lt < 0.05,
        ]) %||%
          0L,
        0L
      )
    )
  ) |>
  pivot_longer(
    cols = c(positive, negative),
    names_to = "Association",
    values_to = "Count"
  ) |>
  ggplot() +
  geom_histogram(
    aes(x = Count, fill = Association),
    bins = 3,
    position = "dodge"
  ) +
  xlab("Number of Associations per Site/Time") +
  ggtitle(
    bquote(bold("B.") ~ "Species Associations in" ~ italic("Daphnia"))
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#EF7C12FF", "#007BC3FF")) +
  scale_x_continuous(breaks = c(0, 1, 2))
p1 / p2
ggsave(filename = "Figures/Fig4_empirical_summary.png", width = 10, height = 10)
