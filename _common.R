library(simhelpers)
library(rsimsum)
library(metafor)
library(clubSandwich)
library(robumeta)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(knitr)

theme_set(theme_minimal(base_size = 13))

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 7,
  fig.height = 5,
  out.width = "85%",
  message = FALSE,
  warning = FALSE
)
