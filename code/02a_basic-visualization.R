
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))

# basic theme
boxplot_theme <- list(
    theme_bw(),
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_text(size = 26),
          plot.title = element_text(size = 30),
          axis.title = element_blank(),
          plot.title.position = "plot")
)

distribution_theme <- list(
  theme_bw(),
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        plot.title.position = "plot")
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. functions -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This is a function to get the pairwise comparisons from the Tukey HSD into
# a matrix for easier comparisons within pairs. Thanks to IRTFM for this
# solution on how to fill in a matrix with the values from a vector:
# https://stackoverflow.com/questions/24598710/presenting-tukey-hsd-pairwise-p-values-in-a-table

# The function expects the subsetted output of a `TukeyHSD()` call. It's
# easier to see in the downstream section where I do the ANOVAs with the post-
# hoc comparisons, but the function's output creates this structure where the
# actual information is in a separate column, like `TukeyHSD(data)$name`.
# In any case, this function expects the `$` object, not the full `TukeyHSD()`
# output.

pairwise_sig_fxn <- function(pairwise) {
  # create an empty matrix
  sig_mat <- matrix(NA, 9, 9)
  
  # turn the pairwise object into a data frame
  sig_df <- data.frame(pairwise) %>% 
    # create a nicer p-value display
    mutate(p_display = case_when(
      # if p-value < 0.001, then display < 0.001
      between(p.adj, 0, 0.001) ~ "<0.001",
      # if p-value is between 0.001 and 0.01, round to 3 digits
      between(p.adj, 0.001, 0.01) ~ as.character(round(p.adj, digits = 3)),
      # if p-value is between 0.01 and 1, round to 2 digits
      between(p.adj, 0.01, 1) ~ as.character(round(p.adj, digits = 2))
    ))
  
  # fill in the matrix with the p-value display
  sig_mat[lower.tri(sig_mat)] <- sig_df$p_display
  
  # make the columns and row names the species codes
  colnames(sig_mat) <- algae_spcode_factors
  rownames(sig_mat) <- algae_spcode_factors
  
  # make the matrix a data frame (for easy table making later)
  sig_mat_to_df <- data.frame(sig_mat) %>% 
    # replace the NAs with a -
    mutate(across(where(is.character), ~replace_na(.,"-")))
  
  # return the final data frame
  return(sig_mat_to_df)
}

sig_table_fxn <- function(pairwise_sig_df) {
  table <- gt(pairwise_sig_df, rownames_to_stub = T) %>%
    cols_width(everything() ~ px(45)) %>% 
    tab_style(
      style = cell_text(align = "center"),
      locations = list(cells_body(columns = everything(),
                                  rows = everything()),
                       cells_column_labels(columns = everything()),
                       cells_stub(rows = everything()))) %>% 
    tab_options(
      table.border.top.color = "white",
      table.font.size = px(26)
    )
  
  return(table)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------- 2. ANOVAs, post-hoc comparisons, and plots --------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to do one-way ANOVAs comparing trait values 
# between species. I also do a Tukey HSD post-hoc comparison after each ANOVA.
# The trait values used in the correlation analysis, standardized major axis
# regression, and PCA are log transformed; thus, this code includes ANOVAs for
# log transformed (log) and raw (raw) trait values. 
# 
# The code also includes visualizations of logged and raw trait values in 
# boxplots and means and 95% CI. 

# ⟞ a. ANOVA --------------------------------------------------------------

sp_anovas <- ind_traits_filtered_full %>% 
  filter(!(specimen_ID %in% c("20210721-BULL-023", "20210719-IVEE-009"))) %>% 
  filter(sp_code %in% c("BO", "CC", "CO", "BF", "DP", 
                        "LAFA", "PTCA", "R", "CYOS")) %>% 
  mutate(sp_code = fct_relevel(sp_code, algae_spcode_factors),
         scientific_name = fct_relevel(scientific_name, algae_factors)) %>% 
  #  H, T, SA, H:WW, DW:WW, H:V, SA:V, SA:DW, and SA:P
  drop_na(specimen_ID, scientific_name, sp_code, 
          maximum_height, 
          thickness_mm_mean, 
          frond_area_scaled,
          height_ww,
          total_dmc, 
          height_vol,
          sav_scaled, 
          sta_scaled,
          sap_mean) %>% 
  select(specimen_ID, scientific_name, sp_code, date_collected, year, site,
         maximum_height, 
         thickness_mm_mean, 
         frond_area_scaled,
         height_ww,
         total_dmc, 
         height_vol,
         sav_scaled, 
         sta_scaled,
         sap_mean
  ) %>% 
  mutate(sp_label = case_match(
    scientific_name,
    "Corallina officinalis" ~ "Corallina officinalis",
    "Bossiella orbigniana" ~ "Bossiella orbigniana",
    "Cryptopleura ruprechtiana" ~ "Cryptopleura ruprechtiana",
    "Chondracanthus corymbiferus; Chondracanthus exasperatus" ~ "Chondracanthus spp.",
    "Rhodymenia californica" ~ "Rhodymenia californica",
    "Stephanocystis osmundacea" ~ "Stephanocystis osmundacea",
    "Dictyota binghamiae; Dictyota flabellata; Dictyota coriacea" ~ "Dictyota spp.",
    "Pterygophora californica" ~ "Pterygophora californica",
    "Laminaria farlowii" ~ "Laminaria farlowii"
  )) %>% 
  pivot_longer(cols = maximum_height:sap_mean,
               names_to = "trait",
               values_to = "value") %>% 
  nest(.by = trait,
       data = everything()) %>% 
  mutate(units = map(
    trait,
    ~ case_match(
      .x,
      "maximum_height" ~ "Maximum height (cm)", 
      "thickness_mm_mean" ~ "Thickness (mm)", 
      "frond_area_scaled" ~ "Surface area (mm\U00B2)",
      "height_ww" ~ "Height:Wet weight (cm/wet mg)",
      "total_dmc" ~ "Dry:wet weight", 
      "height_vol" ~ "Height:volume (cm/mL)",
      "sav_scaled" ~ "Surface area:volume (mm\U00B2/mL)", 
      "sta_scaled" ~ "Surface area:dry weight (mm\U00B2/dry mg)",
      "sap_mean" ~ "Surface area:perimeter (mm\U00B2/mm)"
    )
  )) %>% 
  # log transform the trait value
  mutate(data_log = map(
    data,
    ~ .x %>% 
      mutate(log_value = log(value)) %>% 
      mutate(sp_code = fct_relevel(sp_code, algae_spcode_factors))
  )) %>% 
  mutate(data = map(
    data,
    ~ .x %>% 
      mutate(sp_code = fct_relevel(sp_code, algae_spcode_factors))
  )) %>% 
  
  # calculate variances (ANOVA assumption check)
  # I don't do this for the raw trait values because the variances are all off
  mutate(variance_log = map(
    data_log,
    ~ .x %>% 
      group_by(sp_code) %>% 
      summarize(var = var(log_value))
  )) %>% 
  mutate(variance_test_untransform = map(
    data,
    ~ car::leveneTest(value ~ sp_code,
          data = .x,
          center = "mean") %>% 
      tibble::as_tibble(rownames = "term")
  )) %>% 
  mutate(variance_test_log = map(
    data_log,
    ~ car::leveneTest(log_value ~ sp_code,
                      data = .x,
                      center = "mean") %>% 
      tibble::as_tibble(rownames = "term")
  )) %>% 
  
  # ANOVA using log transformed trait values
  mutate(log_anova = map(
    data_log,
    ~ aov(log_value ~ sp_code,
          data = .x) 
  )) %>% 
  # ANOVA using raw trait values
  mutate(raw_anova = map(
    data,
    ~ aov(value ~ sp_code,
          data = .x) 
  )) %>% 
  
  # tidy ANOVA outputs for table making later
  mutate(log_anova_table = map(
    log_anova,
    ~ tidy(.x)
  )) %>% 
  mutate(raw_anova_table = map(
    raw_anova,
    ~ tidy(.x)
  )) %>% 
  
  # pairwise comparisons
  mutate(log_pairwise = map(
    log_anova,
    ~ TukeyHSD(.x)$sp_code
  )) %>% 
  mutate(raw_pairwise = map(
    raw_anova,
    ~ TukeyHSD(.x)$sp_code
  )) %>% 
  
  # CLD and means
  mutate(log_cld = map(
    log_pairwise,
    ~ .x %>% 
      as_tibble(rownames = "comparisons") %>% 
      clean_names() %>% 
      cldList(p_adj ~ comparisons,
              data = .,
              threshold = 0.05,
              remove.space = FALSE) %>% 
      rename_at(vars(Letter, MonoLetter), ~ paste0(., '_log')) %>% 
      rename("sp_code" = "Group")
  )) %>% 
  mutate(raw_cld = map(
    raw_pairwise,
    ~ .x %>% 
      as_tibble(rownames = "comparisons") %>% 
      clean_names() %>% 
      cldList(p_adj ~ comparisons,
              data = .,
              threshold = 0.05,
              remove.space = FALSE) %>% 
      rename_at(vars(Letter, MonoLetter), ~ paste0(., '_raw')) %>% 
      rename("sp_code" = "Group")
  )) %>% 
  mutate(means_cld = pmap(
    list(x = data_log, y = log_cld, z = raw_cld),
    function(x, y, z) x %>% 
      group_by(sp_code) %>% 
      summarize(mean_log = mean(log_value, na.rm = TRUE),
                max_log = max(log_value, na.rm = TRUE),
                mean_raw = mean(value, na.rm = TRUE),
                max_raw = max(value, na.rm = TRUE)) %>% 
      ungroup() %>% 
      mutate(max_log = max(max_log),
             max_raw = max(max_raw)) %>% 
      left_join(., y, by = "sp_code") %>% 
      left_join(., z, by = "sp_code")
  )) %>% 
  
  # getting the p-values from pairwise comparisons
  mutate(log_sig_df = map(
    log_pairwise,
    ~ pairwise_sig_fxn(.x)
  )) %>% 
  mutate(raw_sig_df = map(
    raw_pairwise,
    ~ pairwise_sig_fxn(.x)
  )) %>% 
  
  # getting the p-values into a table (using gt)
  mutate(log_sig_table = map(
    log_sig_df,
    ~ sig_table_fxn(.x)
  )) %>% 
  mutate(raw_sig_table = map(
    raw_sig_df,
    ~ sig_table_fxn(.x)
  )) %>% 
  
  # log transformed boxplot
  mutate(log_boxplot = map2(
    data, units,
    ~ .x %>% 
      ggplot(aes(x = sp_code,
                 y = log(value))) +
      geom_boxplot(aes(fill = sp_code),
                   outliers = FALSE) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      scale_fill_manual(values = algae_spcode_colors) +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = paste0(.y, " (log transform)")) +
      boxplot_theme
  )) %>% 
  # log transformed means plot (mean + 95% CI)
  mutate(log_means_plot = pmap(
    list(x = data, y = means_cld), 
    function(x, y) ggplot(
      data = x,
      aes(x = sp_code,
          y = log(value),
          group = sp_code,
          color = sp_code)) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      stat_summary(geom = "pointrange",
                   fun.data = mean_cl_normal) +
      geom_text(data = y,
                aes(x = sp_code,
                    y = max_log*1.1,
                    label = Letter_log),
                color = "#000000",
                size = 8) +
      scale_color_manual(values = algae_spcode_colors) +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = "(a) log transformed") +
      boxplot_theme
  )) %>% 
  # raw boxplot
  mutate(raw_boxplot = map2(
    data, units,
    ~ .x %>% 
      ggplot(aes(x = sp_code,
                 y = value)) +
      geom_boxplot(aes(fill = sp_code),
                   outliers = FALSE) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      scale_fill_manual(values = algae_spcode_colors) +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = .y) +
      boxplot_theme
  )) %>% 
  # raw means plot (means + 95% CI)
  mutate(raw_means_plot = pmap(
    list(x = data, y = means_cld), 
    function(x, y) ggplot(
      data = x,
      aes(x = sp_code,
          y = value,
          group = sp_code,
          color = sp_code)
    ) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      stat_summary(geom = "pointrange",
                   fun.data = mean_cl_normal) +
      geom_text(data = y,
                aes(x = sp_code,
                    y = max_raw*1.1,
                    label = Letter_raw),
                color = "#000000",
                size = 8) +
      scale_color_manual(values = algae_spcode_colors) +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = "(b) untransformed") +
      boxplot_theme
  )) %>% 
  mutate(log_means_plot_table = map2(
    log_means_plot, log_sig_table,
    ~ (.x +
      theme(# axis.text = element_text(size = 10),
            axis.text.x = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), units = "cm"))) / wrap_table(.y, space = "fixed")
  )) %>% 
  mutate(raw_means_plot_table = map2(
    raw_means_plot, raw_sig_table,
    ~ (.x +
      theme(# axis.text = element_text(size = 10),
            axis.text.x = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), units = "cm"))) / wrap_table(.y, space = "fixed")
  )) %>% 
  mutate(means_plots_together = pmap(
    list(x = log_means_plot_table, y = raw_means_plot_table, z = units),
    function(x, y, z) (x | y) +
      plot_annotation(
        title = z
      ) & 
      theme(plot.title = element_text(size = 32))
  ))

# h:ww not sig with raw values, sig with log because of huge BO outlier

# maximum height
log_h_means_plot_table <- pluck(sp_anovas, 25, 1)
raw_h_means_plot_table <- pluck(sp_anovas, 26, 1)
h_means_plots_together <- pluck(sp_anovas, 27, 1)

# thickness
log_t_means_plot_table <- pluck(sp_anovas, 25, 2)
raw_t_means_plot_table <- pluck(sp_anovas, 26, 2)
t_means_plots_together <- pluck(sp_anovas, 27, 2)

# surface area
log_sa_means_plot_table <- pluck(sp_anovas, 25, 3)
raw_sa_means_plot_table <- pluck(sp_anovas, 26, 3)
sa_means_plots_together <- pluck(sp_anovas, 27, 3)

# height:wet weight
log_h_ww_means_plot_table <- pluck(sp_anovas, 25, 4)
raw_h_ww_means_plot_table <- pluck(sp_anovas, 26, 4)
h_ww_means_plots_together <- pluck(sp_anovas, 27, 4)

# dry:wet weight
log_dw_ww_means_plot_table <- pluck(sp_anovas, 25, 5)
raw_dw_ww_means_plot_table <- pluck(sp_anovas, 26, 5)
dw_ww_means_plots_together <- pluck(sp_anovas, 27, 5)

# height:volume
log_h_v_means_plot_table <- pluck(sp_anovas, 25, 6)
raw_h_v_means_plot_table <- pluck(sp_anovas, 26, 6)
h_v_means_plots_together <- pluck(sp_anovas, 27, 6)

# SA:V
log_sa_v_means_plot_table <- pluck(sp_anovas, 25, 7)
raw_sa_v_means_plot_table <- pluck(sp_anovas, 26, 7)
sa_v_means_plots_together <- pluck(sp_anovas, 27, 7)

# surface area:dry weight
log_sa_dw_means_plot_table <- pluck(sp_anovas, 25, 8)
raw_sa_dw_means_plot_table <- pluck(sp_anovas, 26, 8)
sa_dw_means_plots_together <- pluck(sp_anovas, 27, 8)

# SA:P
log_sa_p_means_plot_table <- pluck(sp_anovas, 25, 9)
raw_sa_p_means_plot_table <- pluck(sp_anovas, 26, 9)
sa_p_means_plots_together <- pluck(sp_anovas, 27, 9)


# ⟞ b. ANOVA tables -------------------------------------------------------

# ⟞ ⟞ i. wrangling functions ----------------------------------------------

# solution from TarJae: https://stackoverflow.com/questions/72451868/flextable-scientific-formats-for-a-table-that-has-both-very-large-and-very-smal

table_fxn <- function(df) {
  df %>% 
  mutate(term = case_when(
    term == "sp_code" ~ "Species",
    TRUE ~ term
  )) %>% 
    # create a new column that indicates whether an effect is significant
    mutate(signif = case_when(
      p.value < 0.05 ~ "yes",
      TRUE ~ "no"
    )) %>% 
    # create a p-value column that converts very small values to < 0.001
    # and rounds all other values to relevant digits
    mutate(p.value = case_when(
      between(p.value, 0, 0.001) ~ "<0.001",
      between(p.value, 0.001, 0.01) ~ as.character(round(p.value, digits = 3)),
      between(p.value, 0.01, 1) ~ as.character(round(p.value, digits = 2))
    )) %>% 
    # round other numeric values to two digits
    mutate(across(where(is.numeric), ~ round(., digits = 2))) %>% 
    # format sum squares and mean squares in scientific format if > 1000000
    mutate(across(sumsq:meansq, ~
                    case_when(. < 1000000 ~ format(., scientific = FALSE),
                              . >= 1000000 ~ format(., scientific = TRUE, digits = 2))))
}

variance_table_fxn <- function(df) {
  df %>% 
    mutate(term = case_when(
      term == "group" ~ "Species",
      TRUE ~ "Residuals"
    )) %>% 
    # create a new column that indicates whether an effect is significant
    mutate(signif = case_when(
      `Pr(>F)` < 0.05 ~ "yes",
      TRUE ~ "no"
    )) %>% 
    # create a p-value column that converts very small values to < 0.001
    # and rounds all other values to relevant digits
    mutate(p.value = case_when(
      between(`Pr(>F)`, 0, 0.001) ~ "<0.001",
      between(`Pr(>F)`, 0.001, 0.01) ~ as.character(round(`Pr(>F)`, digits = 3)),
      between(`Pr(>F)`, 0.01, 1) ~ as.character(round(`Pr(>F)`, digits = 2))
    )) %>%
    # # round other numeric values to two digits
    mutate(across(where(is.numeric), ~ round(., digits = 2)))
}

sp_anova_tables <- sp_anovas %>% 
  select(units, log_anova_table, raw_anova_table) %>% 
  mutate(log_anova_table = map2(
    log_anova_table, units,
    ~ .x %>% 
      mutate(trait = .y) %>% 
      relocate(trait, .before = term) %>% 
      table_fxn()
  )) %>% 
  mutate(raw_anova_table = map2(
    raw_anova_table, units,
    ~ .x %>% 
      mutate(trait = .y) %>% 
      relocate(trait, .before = term) %>% 
      table_fxn()
  ))


# ⟞ ⟞ ii. generating tables -----------------------------------------------

# ANOVA using raw trait values
raw_anova_table_full <- sp_anova_tables %>% 
  select(raw_anova_table) %>% 
  unnest(cols = raw_anova_table) %>% 
  flextable(col_keys = c("trait",
                         "term",
                         "df",
                         "sumsq",
                         "meansq",
                         "statistic", 
                         "p.value")) %>% 
  # merge group cells to create a grouping column
  merge_v(j = ~ trait) %>% 
  valign(j = ~ trait,
         i = NULL,
         valign = "top") %>% 
  # change the column names
  set_header_labels("trait" = "Trait",
                    "term" = "Term",
                    "df" = "Degrees of freedom",
                    "sumsq" = "Sum of squares",
                    "meansq" = "Mean squares",
                    "statistic" = "F-statistic", 
                    "p.value" = "p-value") %>% 
  # bold p values if they are significant
  style(i = ~ signif == "yes",
        j = c("p.value"),
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  autofit() %>% 
  fit_to_width(7) %>% 
  font(fontname = "Times New Roman",
       part = "all")

# ANOVA using log transformed trait values

log_anova_table_full <- sp_anova_tables %>% 
  select(log_anova_table) %>% 
  unnest(cols = log_anova_table) %>% 
  flextable(col_keys = c("trait",
                         "term",
                         "df",
                         "sumsq",
                         "meansq",
                         "statistic", 
                         "p.value")) %>% 
  # merge group cells to create a grouping column
  merge_v(j = ~ trait) %>% 
  valign(j = ~ trait,
         i = NULL,
         valign = "top") %>% 
  # change the column names
  set_header_labels("trait" = "Trait",
                    "term" = "Term",
                    "df" = "Degrees of freedom",
                    "sumsq" = "Sum of squares",
                    "meansq" = "Mean squares",
                    "statistic" = "F-statistic", 
                    "p.value" = "p-value") %>% 
  # bold p values if they are significant
  style(i = ~ signif == "yes",
        j = c("p.value"),
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  autofit() %>% 
  fit_to_width(7) %>% 
  font(fontname = "Times New Roman",
       part = "all")

# variance comparisons

sp_variance_tables_full <- sp_anovas %>% 
  select(units, variance_test_untransform, variance_test_log) %>% 
  mutate(variance_test_log_table = map2(
    variance_test_log, units,
    ~ .x %>% 
      mutate(trait = .y) %>% 
      relocate(trait, .before = term) %>% 
      variance_table_fxn()
  )) %>% 
  mutate(variance_test_untransform_table = map2(
    variance_test_untransform, units,
    ~ .x %>% 
      mutate(trait = .y) %>% 
      relocate(trait, .before = term) %>% 
      variance_table_fxn()
  ))

raw_variance_table_full <- sp_variance_tables %>% 
  select(variance_test_untransform_table) %>% 
  unnest(cols = variance_test_untransform_table) %>% 
  flextable(col_keys = c("trait",
                         "term",
                         "Df",
                         "F value",
                         "p.value")) %>% 
  # merge group cells to create a grouping column
  merge_v(j = ~ trait) %>% 
  valign(j = ~ trait,
         i = NULL,
         valign = "top") %>% 
  # change the column names
  set_header_labels("trait" = "Trait",
                    "term" = "Term",
                    "Df" = "Degrees of freedom",
                    "F value" = "F-statistic", 
                    "p.value" = "p-value") %>% 
  # bold p values if they are significant
  style(i = ~ signif == "yes",
        j = c("p.value"),
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  # autofit() %>% 
  fit_to_width(5, unit = "in") %>% 
  font(fontname = "Times New Roman",
       part = "all")

log_variance_table_full <- sp_variance_tables %>% 
  select(variance_test_log_table) %>% 
  unnest(cols = variance_test_log_table) %>% 
  flextable(col_keys = c("trait",
                         "term",
                         "Df",
                         "F value",
                         "p.value")) %>% 
  # merge group cells to create a grouping column
  merge_v(j = ~ trait) %>% 
  valign(j = ~ trait,
         i = NULL,
         valign = "top") %>% 
  # change the column names
  set_header_labels("trait" = "Trait",
                    "term" = "Term",
                    "Df" = "Degrees of freedom",
                    "F value" = "F-statistic", 
                    "p.value" = "p-value") %>% 
  # bold p values if they are significant
  style(i = ~ signif == "yes",
        j = c("p.value"),
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  # autofit() %>% 
  fit_to_width(5, unit = "in") %>% 
  font(fontname = "Times New Roman",
       part = "all")


# ⟞ d. Kruskal-Wallis and Dunn --------------------------------------------

# thickness, SA:P, and DW:WW

sp_kw <- sp_anovas %>% 
  select(trait, units, data_log) %>% 
  filter(trait %in% c("thickness_mm_mean", "sap_mean", "total_dmc")) %>% 
  mutate(log_kw = map(
    data_log,
    ~ kruskal.test(value ~ sp_code,
                   data = .x)
  )) %>% 
  mutate(log_dunn = map(
    data_log,
    ~ dunn_test(value ~ sp_code,
               data = .x,
               p.adjust.method = "BH")
  )) %>% 
  mutate(cld = map(
    log_dunn,
    ~ .x %>% 
      unite("Comparison", group1, group2, sep = "-", remove = FALSE) %>%
      cldList(p.adj ~ Comparison, 
              data = .,
              threshold = 0.05,
              remove.space = FALSE) %>% 
      rename("sp_code" = "Group")
  )) %>% 
  mutate(medians_cld = pmap(
    list(x = data_log, y = cld),
    function(x, y) x %>% 
      group_by(sp_code) %>% 
      summarize(median_log = median(log_value, na.rm = TRUE),
                max_log = max(log_value, na.rm = TRUE),
                median_raw = median(value, na.rm = TRUE),
                max_raw = max(value, na.rm = TRUE)) %>%
      ungroup() %>% 
      mutate(max_log = max(max_log),
             max_raw = max(max_raw)) %>% 
      left_join(., y, by = "sp_code")
  )) %>% 
  # log transformed means plot (median as line)
  mutate(log_medians_plot = pmap(
    list(x = data_log, y = medians_cld), 
    function(x, y) ggplot(
      data = x, 
      aes(x = sp_code,
          y = log_value,
          group = sp_code,
          color = sp_code)) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      geom_errorbar(data = y,
                    aes(y = median_log,
                        ymax = median_log,
                        ymin = median_log),
                    linewidth = 0.5) +
      geom_text(data = y,
                aes(x = sp_code,
                    label = Letter,
                    y = max_log*1.1),
                color = "#000000",
                size = 8) +
      scale_color_manual(values = algae_spcode_colors) +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = "(a) log transformed") +
      boxplot_theme
  )) %>% 
  # raw transformed means plot (median as line)
  mutate(raw_medians_plot = pmap(
    list(x = data_log, y = medians_cld), 
    function(x, y) ggplot(
      data = x,
      aes(x = sp_code,
          y = value,
          group = sp_code,
          color = sp_code)
    ) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      geom_errorbar(data = y,
                    aes(y = median_raw,
                        ymax = median_raw,
                        ymin = median_raw),
                    linewidth = 0.5) +
      geom_text(data = y,
                aes(x = sp_code,
                    label = Letter,
                    y = max_raw*1.1),
                color = "#000000",
                size = 8) +
      scale_color_manual(values = algae_spcode_colors) +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = "(b) untransformed") +
      boxplot_theme
  )) %>% 
  mutate(dunn_table = map(
    log_dunn,
    ~ .x %>% 
      select(group1, group2, p.adj) %>% 
      # create a nicer p-value display
      mutate(p_display = case_when(
        # if p-value < 0.001, then display < 0.001
        between(p.adj, 0, 0.001) ~ "<0.001",
        # if p-value is between 0.001 and 0.01, round to 3 digits
        between(p.adj, 0.001, 0.01) ~ as.character(round(p.adj, digits = 3)),
        # if p-value is between 0.01 and 1, round to 2 digits
        between(p.adj, 0.01, 1) ~ as.character(round(p.adj, digits = 2))
      )) %>% 
      mutate(group1 = fct_relevel(group1, algae_spcode_factors),
             group2 = fct_relevel(group2, algae_spcode_factors)) %>%
      select(!p.adj) %>% 
      pivot_wider(names_from = group1, values_from = p_display) %>% 
      add_row(CO = NA, .before = 1) %>% 
      mutate(group2 = case_when(
        is.na(group2) ~ "CO",
        TRUE ~ group2
      )) %>% 
      # replace the NAs with a -
      mutate(across(where(is.character), ~ replace_na(., "-"))) %>%
      column_to_rownames("group2") %>%
      mutate(LAFA = "-") %>% 
      sig_table_fxn()
  )) %>% 
  mutate(log_medians_plot_table = pmap(
    list(w = log_medians_plot, x = raw_medians_plot, y = dunn_table, z = units),
    function(w, x, y, z) w + x + y +
      plot_layout(design = "AABB
                            #CC#") +
      plot_annotation(
        title = z
      ) & 
      theme(plot.title = element_text(size = 32))
  ))

# throws a warning from the fct_relevel - ok to ignore

t_median_plots_together <- pluck(sp_kw, 11, 1)
dw_ww_median_plots_together <- pluck(sp_kw, 11, 2)
sa_p_median_plots_together <- pluck(sp_kw, 11, 3)


# ⟞ e. KW tables ----------------------------------------------------------

kw_tables <- sp_kw %>% 
  select(units, log_kw) %>% 
  mutate(log_kw_table = map(
    log_kw,
    ~ tidy(.x)
  )) %>% 
  select(units, log_kw_table) %>% 
  unnest(cols = c(units, log_kw_table)) %>% 
  select(!method) %>% 
  mutate(p.value = case_when(
    p.value < 0.001 ~ "< 0.001"
  )) %>% 
  mutate(statistic = round(statistic, digits = 2)) %>% 
  relocate(parameter, .before = statistic) %>% 
  flextable() %>% 
  set_header_labels("units" = "Trait",
                    "parameter" = "Degrees of freedom",
                    "statistic" = "\U01D712\U00B2",
                    "p.value" = "p-value") %>% 
  # bold p values if they are significant
  style(j = c("p.value"),
        pr_t = officer::fp_text(bold = TRUE),
        part = "body") %>% 
  autofit() %>% 
  fit_to_width(5, unit = "in") %>% 
  font(fontname = "Times New Roman",
       part = "all")

kw_tables

# ⟞ f. saving outputs -----------------------------------------------------

# ⟞ ⟞ i. plots ------------------------------------------------------------

trait_file_names <- list(
  "h",
  "t",
  "sa",
  "h-ww",
  "dw-ww",
  "h-v",
  "sa-v",
  "sa-dw",
  "sa-p"
)

log_means_plot_tables <- list(
  log_h_means_plot_table,
  log_t_means_plot_table,
  log_sa_means_plot_table,
  log_h_ww_means_plot_table,
  log_dw_ww_means_plot_table,
  log_h_v_means_plot_table,
  log_sa_v_means_plot_table,
  log_sa_dw_means_plot_table,
  log_sa_p_means_plot_table
)

raw_means_plot_tables <- list(
  raw_h_means_plot_table,
  raw_t_means_plot_table,
  raw_sa_means_plot_table,
  raw_h_ww_means_plot_table,
  raw_dw_ww_means_plot_table,
  raw_h_v_means_plot_table,
  raw_sa_v_means_plot_table,
  raw_sa_dw_means_plot_table,
  raw_sa_p_means_plot_table
)

means_plots_together <- list(
  h_means_plots_together,
  t_means_plots_together,
  sa_means_plots_together,
  h_ww_means_plots_together,
  dw_ww_means_plots_together,
  h_v_means_plots_together,
  sa_v_means_plots_together,
  sa_dw_means_plots_together,
  sa_p_means_plots_together
)

# for(i in 1:length(log_means_plot_tables)) {
#   ggsave(here("figures",
#               "basic-visualizations",
#               "means-plots",
#               paste0("log_", trait_file_names[[i]], "_plot-table_", today(), ".jpg")),
#          plot = log_means_plot_tables[[i]],
#          width = 21,
#          height = 21,
#          units = "cm",
#          dpi = 300)
# }
# 
# for(i in 1:length(raw_means_plot_tables)) {
#   ggsave(here("figures",
#               "basic-visualizations",
#               "means-plots",
#               paste0("raw_", trait_file_names[[i]], "_plot-table_", today(), ".jpg")),
#          plot = raw_means_plot_tables[[i]],
#          width = 21,
#          height = 21,
#          units = "cm",
#          dpi = 300)
# }
# 
# for(i in 1:length(means_plots_together)) {
#   ggsave(here("figures",
#               "basic-visualizations",
#               "means-plots",
#               paste0("means-plots-together_", trait_file_names[[i]], "_", today(), ".jpg")),
#          plot = means_plots_together[[i]],
#          width = 25,
#          height = 16,
#          units = "cm",
#          dpi = 250)
# }

kw_trait_file_names <- list(
  "t",
  "dw-ww",
  "sa-p"
)

kw_plots_together <- list(
  t_median_plots_together,
  dw_ww_median_plots_together, 
  sa_p_median_plots_together 
)

# for(i in 1:length(kw_plots_together)) {
#   ggsave(here("figures",
#               "basic-visualizations",
#               "means-plots",
#               paste0("medians-plots-together_", kw_trait_file_names[[i]], "_", today(), ".jpg")),
#          plot = kw_plots_together[[i]],
#          width = 25,
#          height = 14,
#          units = "cm",
#          dpi = 250)
# }



# ⟞ ⟞ ii. tables ----------------------------------------------------------

# save_as_docx(path = here("tables",
#                   "ANOVA",
#                   paste0("raw-trait-ANOVA_", today(), ".docx")),
#              raw_anova_table_full)
# 
# save_as_docx(path = here("tables",
#                          "ANOVA",
#                          paste0("log-trait-ANOVA_", today(), ".docx")),
#              log_anova_table_full)

# save_as_docx(path = here("tables",
#                   "ANOVA",
#                   paste0("raw-trait-ANOVA_variances_", today(), ".docx")),
#              raw_variance_table_full)
# 
# save_as_docx(path = here("tables",
#                          "ANOVA",
#                          paste0("log-trait-ANOVA_variances_", today(), ".docx")),
#              log_variance_table_full)

# save_as_docx(path = here("tables",
#                          "ANOVA",
#                          paste0("kw-tables_", today(), ".docx")),
#              kw_tables)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------------- 2. distributions ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to generate histograms and QQ plots for each 
# trait in three forms: 1) un-transformed, 2) log transformed, and 3) square
# root transformed. Because many of the traits follow a skewed distribution, we
# wanted to know which traits would be candidates for transformation to follow
# a normal distribution.

# ⟞ a. generating plots ---------------------------------------------------

distributions <- sp_anovas %>%
  select(data, units, trait) %>% 
  mutate(length = map(
    data,
    ~ .x %>% 
      drop_na(value) %>% 
      nrow(.)
  )) %>% 
  mutate(bins = map(
    length,
    ~ round((.x^(1/3))*2)
  )) %>% 
  # colors for figures
  mutate(colors = map(
    trait,
    ~ case_match(
      .x,
      "maximum_height" ~ height_col, 
      "thickness_mm_mean" ~ thickness_col, 
      "frond_area_scaled" ~ surface_area_col,
      "height_ww" ~ tradeoff_col,
      "total_dmc" ~ tradeoff_col, 
      "height_vol" ~ tradeoff_col,
      "sav_scaled" ~ tradeoff_col, 
      "sta_scaled" ~ tradeoff_col,
      "sap_mean" ~ tradeoff_col
    )
  )) %>% 
  mutate(distribution_plot = pmap(
    list(w = data, x = bins, y = colors, z = units),
    function(w, x, y, z) ggplot(data = w,
                             aes(x = value)) +
      geom_histogram(bins = x,
                     fill = y,
                     color = "black") +
      labs(title = paste0(z, " (no transformation)"),
           sample = "Value",
           y = "Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      distribution_theme
  )) %>% 
  mutate(distribution_plot_log = pmap(
    list(w = data, x = bins, y = colors, z = units),
    function(w, x, y, z) ggplot(data = w,
                             aes(x = log(value))) +
      geom_histogram(bins = x,
                     fill = y,
                     color = "black") +
      labs(title = paste0(z, " (log transform)"),
           sample = "natural log(value)",
           y = "Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      distribution_theme
  )) %>% 
  mutate(distribution_plot_sqrt = pmap(
    list(w = data, x = bins, y = colors, z = units),
    function(w, x, y, z) ggplot(data = w,
                             aes(x = sqrt(value))) +
      geom_histogram(bins = x,
                     fill = y,
                     color = "black") +
      labs(title = paste0(z, " (square root transform)"),
           sample = "square root (value)",
           y = "Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      distribution_theme
  )) %>% 
  mutate(qq = pmap(
    list(x = data, y = units),
    function(x, y) ggplot(data = x,
                          aes(sample = value)) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21,
              size = 0.5, 
              alpha = 0.5) +
      labs(title = paste0(y, " (no transformation)")) +
      distribution_theme
  )) %>% 
  mutate(qq_log = pmap(
    list(x = data, y = units),
    function(x, y) ggplot(data = x,
                          aes(sample = log(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21,
              size = 0.5, 
              alpha = 0.5) +
      labs(title = paste0(y, " (log transform)")) +
      distribution_theme
  )) %>% 
  mutate(qq_sqrt = pmap(
    list(x = data, y = units),
    function(x, y) ggplot(data = x,
                          aes(sample = sqrt(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21,
              size = 0.5, 
              alpha = 0.5) +
      labs(title = paste0(y, " (square root transform)")) +
      distribution_theme
  ))

# max height
pluck(distributions, 7, 1) # right skewed
pluck(distributions, 10, 1)
pluck(distributions, 8, 1) # log transform looks normal (ish)
pluck(distributions, 11, 1)
pluck(distributions, 9, 1)
pluck(distributions, 12, 1)

# thickness
pluck(distributions, 7, 2) # right skewed
pluck(distributions, 10, 2)
pluck(distributions, 8, 2) # log transform looks normal
pluck(distributions, 11, 2)
pluck(distributions, 9, 2)
pluck(distributions, 12, 2)

# surface area
pluck(distributions, 7, 3) # right skewed
pluck(distributions, 10, 3)
pluck(distributions, 8, 3) # log transform looks normal
pluck(distributions, 11, 3)
pluck(distributions, 9, 3)
pluck(distributions, 12, 3)


# height:wet weight
pluck(distributions, 7, 4) # right skewed
pluck(distributions, 10, 4)
pluck(distributions, 8, 4) # log transform looks normal
pluck(distributions, 11, 4)
pluck(distributions, 9, 4)
pluck(distributions, 12, 4)


# dry:wet weight
pluck(distributions, 7, 5) # right skewed
pluck(distributions, 10, 5)
pluck(distributions, 8, 5) # log transform looks better
pluck(distributions, 11, 5)
pluck(distributions, 9, 5)
pluck(distributions, 12, 5)

# height:volume
pluck(distributions, 7, 6) # right skewed
pluck(distributions, 10, 6)
pluck(distributions, 8, 6) # log transform is normal
pluck(distributions, 11, 6) 
pluck(distributions, 9, 6)
pluck(distributions, 12, 6)


# surface area:volume
pluck(distributions, 7, 7) # right skewed
pluck(distributions, 10, 7)
pluck(distributions, 8, 7) 
pluck(distributions, 11, 7) # log transform is normal
pluck(distributions, 9, 7)
pluck(distributions, 12, 7)


# surface area:dry weight
pluck(distributions, 7, 8) 
pluck(distributions, 10, 8)
pluck(distributions, 8, 8) # log transform is normal
pluck(distributions, 11, 8)
pluck(distributions, 9, 8)
pluck(distributions, 12, 8)


# surface area:perimeter
pluck(distributions, 7, 9) 
pluck(distributions, 10, 9)
pluck(distributions, 8, 9) 
pluck(distributions, 11, 9) # normal
pluck(distributions, 9, 9)
pluck(distributions, 12, 9)

# transform all traits

# ⟞ b. multipanel plot ----------------------------------------------------

# 1 = height
# 2 = thickness
# 3 = SA
# 4 = height:wet weight
# 5 = dry:wet weight
# 6 = height:volume
# 7 = SA:V
# 8 = SA:DW
# 9 = SA:P
# 
# main text: height, thickness, height:volume, height:wet weight, dry:wet weight
# supp: SA, SA:V, SA:DW, SA:P

distributions_log_multipanel_main <- wrap_plots(
  pluck(distributions, 8, 1), 
  pluck(distributions, 8, 2), 
  plot_spacer(),
  pluck(distributions, 8, 6), 
  pluck(distributions, 8, 4), 
  pluck(distributions, 8, 5)
)

distributions_raw_multipanel_main <- wrap_plots(
  pluck(distributions, 7, 1), 
  pluck(distributions, 7, 2), 
  plot_spacer(),
  pluck(distributions, 7, 6), 
  pluck(distributions, 7, 4), 
  pluck(distributions, 7, 5)
)

qq_log_multipanel_main <- wrap_plots(
  pluck(distributions, 11, 1), 
  pluck(distributions, 11, 2), 
  plot_spacer(),
  pluck(distributions, 11, 6), 
  pluck(distributions, 11, 4), 
  pluck(distributions, 11, 5)
)

qq_raw_multipanel_main <- wrap_plots(
  pluck(distributions, 10, 1), 
  pluck(distributions, 10, 2), 
  plot_spacer(),
  pluck(distributions, 10, 6), 
  pluck(distributions, 10, 4), 
  pluck(distributions, 10, 5)
)

# ⟞ c. saving outputs -----------------------------------------------------

# ggsave(here("figures",
#             "basic-visualizations",
#             "distributions",
#             paste0("multipanel_hist_", today(), ".jpg")),
#        distributions_log_multipanel_main,
#        width = 14,
#        height = 8,
#        units = "cm",
#        dpi = 300)
# 
# ggsave(here("figures",
#             "basic-visualizations",
#             "distributions",
#             paste0("multipanel_qq_", today(), ".jpg")),
#        qq_log_multipanel_main,
#        width = 14,
#        height = 8,
#        units = "cm",
#        dpi = 300)
# 
# ggsave(here("figures",
#             "basic-visualizations",
#             "distributions",
#             paste0("multipanel_hist_raw_", today(), ".jpg")),
#        distributions_raw_multipanel_main,
#        width = 14,
#        height = 8,
#        units = "cm",
#        dpi = 300)
# 
# ggsave(here("figures",
#             "basic-visualizations",
#             "distributions",
#             paste0("multipanel_qq_raw_", today(), ".jpg")),
#        qq_raw_multipanel_main,
#        width = 14,
#        height = 8,
#        units = "cm",
#        dpi = 300)

