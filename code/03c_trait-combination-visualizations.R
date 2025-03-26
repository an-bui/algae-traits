# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "03a_correlation-and-tradeoffs.R"))

combo_2traits <- read_rds(file = here("rds-objects",
                                      "combo_2-traits_2025-03-25.rds")) 

combo_3traits <- read_rds(file = here("rds-objects",
                                      "combo_3-traits_2025-03-25.rds")) 

combo_3traits_full <- read_rds(file = here("rds-objects",
                                           "combo_3-traits_2024-12-03.rds"))

combo_4traits <- read_rds(file = here("rds-objects",
                                      "combo_4-traits_2025-03-25.rds")) 

combo_4traits_full <- read_rds(file = here("rds-objects",
                                           "combo_4-traits_2024-12-03.rds"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 1. functions and aesthetics ------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to extract cumulative proportion, conserved pairwise differences,
# and traits
keep_trait_function <- function(df) {
  df %>% 
    select(contains("trait"), cumu_prop, pairwise_padj_conserved_euc) %>% 
    unnest(cols = c(cumu_prop, pairwise_padj_conserved_euc)) %>% 
    rownames_to_column("combo") # %>% 
  # mutate(across(contains("trait"), 
  #               \(x) fct_relevel(x, trait_factor))) 
}

# function to organize the data frame for the heat map (bottom of upset plot)
keep_traits_heatmap_function <- function(df) {
  df %>% 
    # arrange in descending order of cumulative proportion explained
    arrange(-cumu_prop) %>% 
    # make the order of combinations the order in which they appear
    mutate(combo = fct_inorder((combo))) %>% 
    # get rid of the cumulative proportion column
    select(!cumu_prop) %>% 
    # make the data frame longer for plotting
    pivot_longer(cols = contains("trait")) %>% 
    # rename the colume "value" to be trait
    rename(trait = value) %>% 
    # mark the trait as present
    mutate(present = "TRUE") %>% 
    # get rid of the "name" column, which contains trait1, trait2, etc.
    select(!name) %>% 
    # complete all traits within a combination
    complete(combo, trait) %>% 
    # if the trait is not present in the combination, mark "FALSE"
    replace_na(list(present = "FALSE")) # %>% 
  # reorder the traits in factor order
  # mutate(trait = fct_relevel(trait, rev(trait_factor)))
}

# aesthetics for the bottom of the upset plot
upset_plot_bottom <- list(
  geom_tile(color = "black"), 
  scale_x_discrete(expand = c(0, 0)),
  scale_y_discrete(expand = c(0, 0),
                   limits = rev),
  labs(y = "Trait"),
  theme_minimal(),
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20))
)

# aesthetics for the top of the upset plot
upset_plot_top <- list(
  coord_cartesian(ylim = c(0.78, 1)),
  scale_x_discrete(expand = c(0, 0)),
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(0.80, 0.85, 0.90, 0.95, 1.00)),
  labs(y = "Cumulative proportion"),
  theme_minimal(),
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20))
)

# aesthetics for the tally plot
tally_plot_aesthetics <- list(
  labs(x = "Number of occurrences"),
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))),
  theme(panel.background = element_blank(),
        # plot.background = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20),
        text = element_text(size = 20),
        plot.title.position = "plot"
        # plot.margin = margin(0, 0, 0, 0, "pt")
        )
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------- 2. model combinations  ------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ i. 4 traits -----------------------------------------------------------

# assigning each combination "same", "all", or "everything else"
keep_4traits_categories <- combo_4traits %>% 
  mutate(combo = rownames(.)) %>% 
  # count TRUE and FALSE occurrences
  mutate(count = map(
    pairwise_padj_significant_euc,
    ~ .x %>% unlist() %>% table() 
  )) %>% 
  select(combo, pairwise_padj_conserved_euc) %>% 
    mutate(category = case_when(
      pairwise_padj_conserved_euc == "yes" ~ "same",
      TRUE ~ "everything else"
    )) %>% 
  unnest(cols = everything())

# only keep combinations where the pairwise comparisons are conserved
keep_4traits <- combo_4traits %>% 
  keep_trait_function() %>% 
  left_join(., keep_4traits_categories, by = "combo") %>% 
  mutate(alpha = case_when(
    category %in% c("same", "all") ~ "opaque",
    category %in% c("everything else") ~ "transparent"
  ))

keep_4traits_tally <- keep_4traits %>%
  keep_traits_heatmap_function() %>%
  filter(category %in% c("same", "all")) %>% 
  group_by(trait) %>%
  tally() %>%
  ungroup() %>%
  arrange(-n) 

keep_4traits_heatmap <- keep_4traits %>% 
  keep_traits_heatmap_function() %>% 
  mutate(trait = fct_relevel(trait, pull(keep_4traits_tally, trait))) %>% 
  unite("category_present", category, present, sep = "/", remove = FALSE) %>% 
  mutate(alpha = case_when(
    category_present %in% c("same/TRUE", "all/TRUE", "NA/FALSE") ~ "opaque",
    category_present %in% c("everything else/TRUE") ~ "transparent"
  ))

keep_4traits %>% 
  filter(if_any(contains("trait"), 
                ~ str_detect(., "Surface area:perimeter"))) %>% 
  group_by(category) %>% 
  tally()

keep_4traits %>% 
  filter(if_any(contains("trait"), 
                ~ str_detect(., "Surface area|Height"))) %>% 
  group_by(category) %>% 
  tally()

keep_4traits %>% 
  filter(if_any(contains("trait"), 
                ~ str_detect(., "Thickness|Dry:wet weight"))) %>% 
  group_by(category) %>% 
  tally()

bottom_4traits <- keep_4traits_heatmap %>% 
  ggplot(aes(x = combo,
             y = trait,
             fill = category_present,
             alpha = alpha)) +
  scale_fill_manual(values = c("NA/FALSE" = "white",
                               "same/TRUE" = "#356B7A",
                               "everything else/TRUE" = "#9fa0a1",
                               "all/TRUE" = "#059EE6")) +
  scale_alpha_manual(values = c("transparent" = 0.4,
                                "opaque" = 1)) +
  upset_plot_bottom

top_4traits <- ggplot(data = keep_4traits,
                      aes(x = reorder(combo, -cumu_prop),
                          y = cumu_prop,
                          alpha = alpha,
                          fill = category)) +
  geom_col() +
  scale_alpha_manual(values = c("transparent" = 0.4,
                                "opaque" = 1)) +
  scale_fill_manual(values = c("same" = "#356B7A",
                               "everything else" = "#9fa0a1",
                               "all" = "#059EE6")) +
  upset_plot_top

tally_plot_4traits <- keep_4traits_tally %>%
  ggplot(aes(x = n,
             y = reorder(trait, n))) +
  geom_point(size = 1,
             color = "#356B7A") +
  geom_segment(aes(x = 0, xend = n),
               color = "#356B7A") +
  geom_text(aes(x = n + 3.5,
                label = n),
            size = 6) +
  tally_plot_aesthetics +
  theme(axis.title.x = element_blank()) +
  labs(title = "(a) 4 trait combination")

upset_plot_4traits <- (top_4traits / bottom_4traits) # | (plot_spacer() / right_4traits)
upset_plot_4traits

# ggsave(here("figures",
#             "trait-selection",
#             paste0("upset-plot_4traits_all-combinations_euclidean_", today(), ".jpg")),
#        upset_plot_4traits,
#        height = 8,
#        width = 10,
#        units = "cm",
#        dpi = 400)

# ⟞ ii. 3 traits ----------------------------------------------------------

# assigning each combination "same", "all", or "everything else"
keep_3traits_categories <- combo_3traits %>% 
  mutate(combo = rownames(.)) %>% 
  # count TRUE and FALSE occurrences
  mutate(count = map(
    pairwise_padj_significant_euc,
    ~ .x %>% unlist() %>% table() 
  )) %>% 
  select(combo, pairwise_padj_conserved_euc) %>% 
  mutate(category = case_when(
    pairwise_padj_conserved_euc == "yes" ~ "same",
    TRUE ~ "everything else"
  )) %>% 
  unnest(cols = everything())

# only keep combinations where the pairwise comparisons are conserved
keep_3traits <- combo_3traits %>% 
  keep_trait_function() %>% 
  left_join(., keep_3traits_categories, by = "combo") %>% 
  mutate(alpha = case_when(
    category %in% c("same", "all") ~ "opaque",
    category %in% c("everything else") ~ "transparent"
  ))

keep_3traits_tally <- keep_3traits %>%
  keep_traits_heatmap_function() %>%
  filter(category %in% c("same", "all")) %>% 
  group_by(trait) %>%
  tally() %>%
  ungroup() %>%
  arrange(n) 

keep_3traits_heatmap <- keep_3traits %>% 
  keep_traits_heatmap_function() %>% 
  mutate(trait = fct_relevel(trait, pull(keep_3traits_tally, trait))) %>% 
  unite("category_present", category, present, sep = "/", remove = FALSE) %>% 
  mutate(alpha = case_when(
    category_present %in% c("same/TRUE", "all/TRUE", "NA/FALSE") ~ "opaque",
    category_present %in% c("everything else/TRUE") ~ "transparent"
  ))

keep_3traits %>% 
  filter(if_any(contains("trait"), 
                ~ str_detect(., "Surface area|Height"))) %>% 
  group_by(category) %>% 
  tally()

keep_3traits %>% 
  filter(if_any(contains("trait"), 
                ~ str_detect(., "Thickness|Dry:wet weight"))) %>% 
  group_by(category) %>% 
  tally()

bottom_3traits <- keep_3traits_heatmap %>% 
  ggplot(aes(x = combo,
             y = trait,
             fill = category_present,
             alpha = alpha)) +
  scale_fill_manual(values = c("NA/FALSE" = "white",
                               "same/TRUE" = "#077A54",
                               "everything else/TRUE" = "#9fa0a1",
                               "all/TRUE" = "#05E67D")) +
  scale_alpha_manual(values = c("transparent" = 0.4,
                                "opaque" = 1)) +
  upset_plot_bottom

top_3traits <- ggplot(data = keep_3traits,
                      aes(x = reorder(combo, -cumu_prop),
                          y = cumu_prop,
                          alpha = alpha,
                          fill = category)) +
  geom_col() +
  scale_alpha_manual(values = c("transparent" = 0.4,
                                "opaque" = 1)) +
  scale_fill_manual(values = c("same" = "#077A54",
                               "everything else" = "#9fa0a1",
                               "all" = "#05E67D")) +
  upset_plot_top

tally_plot_3traits <- keep_3traits_tally %>%
  ggplot(aes(x = n,
             y = reorder(trait, n))) +
  geom_point(size = 1,
             color = "#077A54") +
  geom_segment(aes(x = 0, xend = n),
               color = "#077A54") +
  geom_text(aes(x = n + 1.5,
                label = n),
            size = 6) +
  tally_plot_aesthetics +
  labs(title = "(b) 3 trait combination")

upset_plot_3traits <- (top_3traits / bottom_3traits) # | (plot_spacer() / right_3traits)

# upset_plot_3traits <- bottom_3traits + plot_spacer() + right_3traits +
#   plot_layout(widths = c(10, -1, 10))

# upset_plot_3traits <- (top_3traits | plot_spacer() | plot_spacer()) /
#   (bottom_3traits | plot_spacer() | right_3traits) +
#   plot_layout(widths = c(4  -10, 4.5))

# upset_plot_3traits <- top_3traits / bottom_3traits
upset_plot_3traits

# ggsave(here("figures",
#             "trait-selection",
#             paste0("upset-plot_3traits_all-combinations_euclidean_", today(), ".jpg")),
#        upset_plot_3traits,
#        height = 8,
#        width = 10,
#        units = "cm",
#        dpi = 400)

tally_plots <- tally_plot_4traits / tally_plot_3traits

# ggsave(here("figures",
#             "trait-selection",
#             paste0("tally-plots_", today(), ".jpg")),
#        tally_plots,
#        height = 8,
#        width = 6,
#        units = "cm",
#        dpi = 400)


# ⟞ iii. 2 traits ---------------------------------------------------------

# only keep combinations where the pairwise comparisons are conserved
keep_2traits <- combo_2traits %>%
  keep_trait_function()
# 
keep_2traits_tally <- keep_2traits %>%
  keep_traits_heatmap_function() %>%
  filter(pairwise_padj_conserved_euc == "yes" & present == TRUE) %>%
  group_by(trait) %>%
  tally() %>%
  ungroup() %>%
  arrange(n) %>%
  mutate(trait = fct_inorder(trait))
# 
keep_2traits_heatmap <- keep_2traits %>%
  keep_traits_heatmap_function() %>%
  mutate(trait = fct_relevel(trait, as.character(pull(keep_2traits_tally, trait))))
# 
bottom_2traits <- keep_2traits_heatmap %>%
  ggplot(aes(x = combo,
             y = trait,
             fill = present,
             alpha = pairwise_padj_conserved_euc)) +
  scale_fill_manual(values = c("TRUE" = "orange",
                               "FALSE" = "white")) +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_bottom
# 
# 
top_2traits <- keep_2traits %>%
  mutate(combo = fct_relevel(combo, as.character(unique(pull(keep_2traits_heatmap, combo))))) %>%
  # filter(pairwise_conserved == "yes") %>%
  ggplot(aes(x = combo,
             y = cumu_prop,
             alpha = pairwise_padj_conserved_euc)) +
  geom_col(fill = "orange") +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_top
# 
# right_2traits <- keep_2traits_tally %>%
#   ggplot(aes(x = n,
#              y = trait)) +
#   geom_point(size = 1,
#              color = "darkgreen") +
#   geom_segment(aes(x = 0, xend = n),
#                color = "darkgreen") +
#   geom_text(aes(x = n + 0.5,
#                 label = n),
#             size = 6) +
#   upset_plot_right +
#   coord_cartesian(xlim = c(9.5, 19.75))
# # 
# # right_2traits
# 
upset_plot_2traits <- (top_2traits / bottom_2traits) # | (plot_spacer() / right_2traits)
# 
# upset_plot_2traits
# 
# ggsave(here("figures",
#             "trait-selection",
#             paste0("upset-plot_2traits_all-combinations_euclidean", today(), ".jpg")),
#        upset_plot_2traits,
#        height = 5,
#        width = 14,
#        units = "cm",
#        dpi = 400)