library(tidyverse)

data <- read_csv("ind-trait-values_2024-11-19.csv")

plot <- ggplot(data = data,
               aes(x = sp_code,
                   y = thickness_mm_mean)) +
  geom_point()

write_rds(x = plot,
          file = "thickness.rds")

# read_plot <- read_rds(here("code", "misc", "height.rds"))
# read_plot
