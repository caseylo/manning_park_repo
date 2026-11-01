require(data.table)
require(tidyverse)
plots <- fread("./data/raw/plot_info_for_test.csv")
library(data.table)
library(ggplot2)

plots <- fread("./data/raw/plot_info_for_test.csv")
plots[, SNR_num :=
        c(A = 1, B = 2, C = 3, D = 4, E = 5)[gsub("[+-]", "", SNR)]
]
plot_long <- melt(
  plots,
  id.vars = "year",
  measure.vars = c("SMR", "elev", "SNR_num"),
  variable.name = "Variable",
  value.name = "Value"
)

ggplot(
  plot_long,
  aes(
    x = factor(year),
    y = Value
  )
) +
  geom_boxplot(
    fill = "lightblue",
    colour = "black",
    outlier.alpha = 0.5
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    x = "Survey Year",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )
