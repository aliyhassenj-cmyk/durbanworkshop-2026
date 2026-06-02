FROM ghcr.io/rocker-org/devcontainer/tidyverse:4.5

# Curated Workshop 1 R packages not already in the tidyverse image.
RUN R -q -e 'install.packages(c("here", "janitor", "gt", "palmerpenguins", "ggbeeswarm", "gtsummary", "ggpubr", "patchwork", "broom.helpers"))'
