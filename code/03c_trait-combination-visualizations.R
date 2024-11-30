
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "03a_correlation-and-tradeoffs.R"))

combo_3traits <- read_rds(file = here("rds-objects",
                                      "combo_3-traits_2024-11-29.rds")) %>% 
  rename("pairwise_padj_conserved" = "pairwise_padj_")

combo_4traits <- read_rds(file = here("rds-objects",
                                      "combo_4-traits_2024-11-29.rds")) %>% 
  rename("pairwise_padj_conserved" = "pairwise_padj_")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 1. functions and aesthetics ------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to extract cumulative proportion, conserved pairwise differences,
# and traits
keep_trait_function <- function(df) {
  df %>% 
    select(contains("trait"), cumu_prop, pairwise_padj_conserved) %>% 
    unnest(cols = c(cumu_prop, pairwise_padj_conserved)) %>% 
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
  scale_y_discrete(expand = c(0, 0)),
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

# aesthetics for the right of the upset plot
upset_plot_right <- list(
  labs(x = "Number of occurrences"),
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))),
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------- 2. model combinations  ------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ i. 4 traits -----------------------------------------------------------

# only keep combinations where the pairwise comparisons are conserved
keep_4traits <- combo_4traits %>% 
  keep_trait_function()

keep_4traits_tally <- keep_4traits %>% 
  keep_traits_heatmap_function() %>% 
  filter(pairwise_padj_conserved == "yes" & present == TRUE) %>% 
  group_by(trait) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(n) %>% 
  mutate(trait = fct_inorder(trait))

keep_4traits_heatmap <- keep_4traits %>% 
  keep_traits_heatmap_function() %>% 
  mutate(trait = fct_relevel(trait, as.character(pull(keep_4traits_tally, trait))))

bottom_4traits <- keep_4traits_heatmap %>% 
  ggplot(aes(x = combo,
             y = trait,
             fill = present,
             alpha = pairwise_padj_conserved)) +
  scale_fill_manual(values = c("TRUE" = "cornflowerblue",
                               "FALSE" = "white")) +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_bottom

top_4traits <- ggplot(data = keep_4traits,
                      aes(x = reorder(combo, -cumu_prop),
                          y = cumu_prop,
                          alpha = pairwise_padj_conserved)) +
  geom_col(fill = "cornflowerblue") +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_top

right_4traits <- keep_4traits_tally %>% 
  ggplot(aes(x = n,
             y = trait)) +
  geom_point(size = 1,
             color = "cornflowerblue") +
  geom_segment(aes(x = 0, xend = n),
               color = "cornflowerblue") +
  geom_text(aes(x = n + 0.5,
                label = n),
            size = 8) +
  upset_plot_right +
  coord_cartesian(xlim = c(31.5, 46.5)) +
  theme(axis.title.x = element_text(size = 14))

upset_plot_4traits <- (top_4traits / bottom_4traits) | (plot_spacer() / right_4traits)
upset_plot_4traits

ggsave(here("figures",
            "trait-selection",
            paste0("upset-plot_4traits_all-combinations_", today(), ".jpg")),
       upset_plot_4traits,
       height = 5,
       width = 19,
       units = "cm",
       dpi = 400)

# top combination: SA:V, SA:DW, H, H:WW
# PCA summary
summary(combo_real[[7]][[44]])
# PERMANOVA summary
combo_4traits[[9]][[44]]
# pairwise summary
combo_4traits[[10]][[44]]

# ⟞ ii. 3 traits ----------------------------------------------------------

# only keep combinations where the pairwise comparisons are conserved
keep_3traits <- combo_3traits %>% 
  keep_trait_function() 

keep_3traits_tally <- keep_3traits_heatmap %>% 
  filter(present == TRUE) %>% 
  group_by(trait, pairwise_padj_conserved) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(trait) %>% 
  mutate(total = sum(n)) %>% 
  ungroup() %>% 
  filter(pairwise_padj_conserved == "yes") %>% 
  arrange(n) %>% 
  mutate(trait = fct_inorder(trait)) %>% 
  mutate(label = paste0(n, "/", total))

keep_3traits_heatmap <- keep_3traits %>% 
  keep_traits_heatmap_function() %>% 
  mutate(trait = fct_relevel(trait, as.character(pull(keep_3traits_tally, trait))))



bottom_3traits <- keep_3traits_heatmap %>% 
  ggplot(aes(x = combo,
             y = trait,
             fill = present,
             alpha = pairwise_padj_conserved)) +
  scale_fill_manual(values = c("TRUE" = "darkgreen",
                               "FALSE" = "white")) +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_bottom


top_3traits <- keep_3traits %>% 
  # filter(pairwise_conserved == "yes") %>% 
  ggplot(aes(x = reorder(combo, -cumu_prop),
             y = cumu_prop,
             alpha = pairwise_padj_conserved)) +
  geom_col(fill = "darkgreen") +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_top

right_3traits <- keep_3traits_tally %>% 
  ggplot(aes(x = n,
             y = trait)) +
  geom_point(size = 1,
             color = "darkgreen") +
  geom_segment(aes(x = 0, xend = n),
               color = "darkgreen") +
  geom_text(aes(x = n + 0.5,
                label = n),
            size = 6) +
  upset_plot_right +
  coord_cartesian(xlim = c(9.5, 19.75)) 

right_3traits
# theme_void() +
# # theme_void() +
# # upset_plot_right +
# theme(legend.position = "right",
#       axis.title.x = element_text(),
#       panel.background = element_blank())

upset_plot_3traits <- (top_3traits / bottom_3traits) | (plot_spacer() / right_3traits)

# upset_plot_3traits <- bottom_3traits + plot_spacer() + right_3traits +
#   plot_layout(widths = c(10, -1, 10))

# upset_plot_3traits <- (top_3traits | plot_spacer() | plot_spacer()) /
#   (bottom_3traits | plot_spacer() | right_3traits) +
#   plot_layout(widths = c(4  -10, 4.5))

# upset_plot_3traits <- top_3traits / bottom_3traits
upset_plot_3traits

ggsave(here("figures",
            "trait-selection",
            paste0("upset-plot_3traits_all-combinations_", today(), ".jpg")),
       upset_plot_3traits,
       height = 5,
       width = 14,
       units = "cm",
       dpi = 400)

# top combination: SA:V, SA:DW, H:WW
# PCA summary
summary(combo_3traits[[6]][[33]])
# PERMANOVA summary
combo_3traits[[8]][[33]]
# pairwise summary
combo_3traits[[9]][[33]]

