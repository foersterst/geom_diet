# R Script: Gradual evolution of larval dietary breadth and its correlated effects
# on body size diversification in geometrid moths
# 13 Sep 2026


# Libraries & Data --------------------------------------------------------

library(tidyverse)
library(phytools)
library(flextable)
library(RColorBrewer)
library(patchwork)
library(ggsci)
library(ggdist)
library(coda)

# read in trait data
dat <- read_csv("data/geom-traits.csv")

# read in tree
phy <- read.nexus("data/geom_tree.nex")

# must be TRUE
identical(dat$species, phy$tip.label)

# vectors for ancestral states
spp <- setNames(dat$specialization, dat$species)
spp <- spp[phy$tip.label]
identical(names(spp), phy$tip.label)


# 1. DESCRIPTIVE STATS ----------------------------------------------------

# proportion species in each dietary state
round(table(dat$specialization) / nrow(dat), 2)

# mean body size per state
dat |>
  group_by(specialization) |>
  summarise(
    bs_x = mean(m_wingspan),
    bs_sd = sd(m_wingspan)
  )


# 2. ANCESTRAL STATES -----------------------------------------------------

## Mk Models ----

# Test different transition matrices (Q)

# Define a backbone Q matrix
qm <- matrix(
  c(
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
  ),
  3, 3,
  byrow = TRUE
)
rownames(qm) <- colnames(qm) <- sort(unique(spp))
qm

# sequential different rates (equal rates is a case of ARD):
q1 <- qm
q1[1, 2] <- 1
q1[2, 3] <- 2

# reverse sequential different rates:
q2 <- qm
q2[2, 1] <- 1
q2[3, 2] <- 2

# bifurcated at mono, irreversible:
q3 <- qm
q3[1, 2] <- 1
q3[1, 3] <- 2

# bifurcated at mono, reversible:
q4 <- q3
q4[2, 1] <- 3
q4[3, 1] <- 4

# bifurcated at olig, irreversible:
q5 <- qm
q5[2, 1] <- 1
q5[2, 3] <- 2

# bifurcated at olig, reversible:
q6 <- q5
q6[1, 2] <- 3
q6[3, 2] <- 4

# bifurcated at poly, irreversible:
q7 <- qm
q7[3, 1] <- 1
q7[3, 2] <- 2

# bifurcated at poly, reversible:
q8 <- q7
q8[1, 3] <- 3
q8[2, 3] <- 4

# release from mono: once you gain the ability to feed on more plants, it gets
# easier to change from oligo to poly and you do not go back to mono:
q9 <- qm
q9[1, 2] <- 1
q9[2, 3] <- 2
q9[3, 2] <- 3

# becoming mono: reverse of q9:
q10 <- qm
q10[3, 2] <- 1
q10[2, 1] <- 2
q10[1, 2] <- 3

# q11: equal rates
# q12: ARD

# list of Q matrices
list_Q <- list(q1, q2, q3, q4, q5, q6, q7, q8, q9, q10)
names(list_Q) <- paste("q", 1:10, sep = "")
rm(qm, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10)

# Rescale tree (% total time): the transition rates (q_ij) represent
# the expected number of state changes per unit of relative tree depth.
phy_rel <- phy
phy_rel$edge.length <- phy_rel$edge.length / max(node.depth.edgelength(phy))

# define the function
fit_custom_Mk <- function(Q) {
  fitMk(tree = phy_rel, x = spp, model = Q, pi = "fitzjohn")
}

# fit the custom models
list_mk <- purrr::map(list_Q, fit_custom_Mk)

# add the ARD and ER models
list_mk[[11]] <- fitMk(tree = phy_rel, x = spp, model = "ARD")
list_mk[[12]] <- fitMk(tree = phy_rel, x = spp, model = "ER")
names(list_mk)[11] <- "ARD"
names(list_mk)[12] <- "ER"

# AIC
sort(sapply(list_mk, AIC))

# AIC Weights
round(sort(aic.w(sapply(list_mk, AIC)), decreasing = TRUE), 3)


## Best-fit MK Model ----

# The best MK model was the one with matrix Q6, followed by ARD.
# I will fit the two models and average them when doing the ancestral states.

# design the matrix q6
q6 <- matrix(
  c(
    0, 3, 0,
    1, 0, 2,
    0, 4, 0
  ),
  3, 3,
  byrow = TRUE
)
rownames(q6) <- colnames(q6) <- sort(unique(spp))

# and fit the Mk model for matrix q6
q6m <- fitMk(tree = phy_rel, x = spp, model = q6, pi = "fitzjohn")

# as well as the Mk model with ARD Q matrix
qARD <- fitMk(tree = phy_rel, x = spp, model = "ARD", pi = "fitzjohn")

# the two best-fit Mk models tell more or less the same story
par(mfrow = c(1, 2))
plot(q6m, color = TRUE, width = TRUE, main = "Matrix Q6")
plot(qARD, color = TRUE, width = TRUE, main = "ARD")
dev.off()

# export
# pdf(file = "out/fitMk_ARD.pdf", width = 7, height = 6)
cols <- pal_jama()(3)
colsC <- viridisLite::viridis(n = 10)
qARD$states <- substr(qARD$states, 1, 1)
xy <- plot(
  qARD,
  show.zeros = TRUE,
  spacer = 0.15,
  mar = rep(0.5, 4),
  color = TRUE,
  lwd = 2.2,
  text = F,
  palette = colsC
)

xy$states <- c("M", "O", "P")

invisible(
  mapply(
    draw.circle,
    xy$x,
    xy$y,
    col = cols,
    border = NA,
    radius = 0.1,
    nv = 10000
  )
)

text(
  xy$x,
  xy$y,
  xy$states,
  cex = 1.2,
  col = "white",
  font = 2
)

# dev.off()

# pdf(file = "out/fitMk_Q6.pdf", width = 7, height = 6)
cols <- pal_jama()(3)
colsC <- viridisLite::viridis(n = 10)
q6m$states <- substr(q6m$states, 1, 1)
xy <- plot(
  q6m,
  show.zeros = TRUE,
  spacer = 0.15,
  mar = rep(0.5, 4),
  color = TRUE,
  lwd = 2.2,
  text = F,
  palette = colsC
)

xy$states <- c("M", "O", "P")

invisible(
  mapply(
    draw.circle,
    xy$x,
    xy$y,
    col = cols,
    border = NA,
    radius = 0.1,
    nv = 10000
  )
)

text(xy$x, xy$y, xy$states,
  cex = 1.2, col = "white",
  font = 2
)

# dev.off()


## ACE ----

# Ancestral state estimation averaging across the two best-fit Mk models
qAVG <- anova(q6m, qARD)
spp_anc <- ancr(qAVG)

# Graph
# Ancestral States
co <- setNames(pal_jama(alpha = 0.75)(3), c("mono", "olig", "poly"))
ts <- to.matrix(spp, sort(unique(spp)))
arv <- setNames(dat$log10_ws, dat$species)

# save
# svg("out/ase-tree-nodes.svg", width = 13/2.54, height = 13/2.54)
plot.phylo(
  x = phy_rel,
  type = "fan",
  show.tip.label = FALSE,
  edge.width = 0.6,
  no.margin = TRUE
)

par(fg = NA)

nodelabels(
  node = 1:phy_rel$Nnode + Ntip(phy_rel),
  pie = spp_anc$ace,
  piecol = co,
  cex = 0.42
)

tiplabels(pie = ts[phy_rel$tip.label, ], piecol = co, cex = 0.2)

# dev.off()

# Tibble to store ancestral states
n <- colnames(spp_anc$ace)
n <- n[apply(spp_anc$ace, 1, which.max)]
p <- apply(spp_anc$ace, 1, max)
dat_ace <- tibble(node = rownames(spp_anc$ace), ace = n, prob = p)
rm(n, p)
# The ancestral states of each node is the one with the highest probability
dat_ace


## Stochastic Maps ----

# stochastic character mapping (time spent in each state)
set.seed(16)
q6Map <- make.simmap(tree = phy, x = spp, model = q6, nsim = 1000, pi = "fitzjohn")
describe.simmap(q6Map)
# 1000 trees with a mapped discrete character with states:
#   mono, olig, poly
#
# trees have 4688.305 changes between states on average
#
# changes are of the following types:
#      mono,olig mono,poly olig,mono olig,poly poly,mono poly,olig
# x->y   2284.91         0   2286.11    60.645         0     56.64
#
# mean total time spent in each state is:
#              mono         olig         poly    total
# raw  1953.8123227 2374.9340602 3195.2898498 7524.036
# prop    0.2596761    0.3156463    0.4246776    1.000

# stochastic character mapping (time spent in each state)
set.seed(16)
qARDmap <- make.simmap(tree = phy, x = spp, model = "ARD", nsim = 1000, pi = "fitzjohn")
describe.simmap(qARDmap)
# 1000 trees with a mapped discrete character with states:
#   mono, olig, poly
#
# trees have 349.806 changes between states on average
#
# changes are of the following types:
#      mono,olig mono,poly olig,mono olig,poly poly,mono poly,olig
# x->y   112.921     8.948   114.899    54.777         0    58.261
#
# mean total time spent in each state is:
#              mono         olig         poly    total
# raw  2062.9447453 2267.4683377 3193.6231497 7524.036
# prop    0.2741806    0.3013633    0.4244561    1.000


# 3. BAYESTRAITS --------------------------------------------------------------------

## Five Partitions ----

# Post-processing of analysis with five branch partitions: mono, oligo, poly, increase,
# decrease.

# read in BayesTraits output
part_BT <- read_table("bt/output/res-5p-log.txt")

# convergence (ESS)
part_BT |>
  select(-c(Iteration, Lh, Tree_No)) |>
  as.mcmc() |>
  effectiveSize()

# data for plotting
part_BT %>%
  mutate(
    remain_mono = log10(Sigma2_1),
    remain_olig = log10(Sigma2_1 * Olig_Branch),
    remain_poly = log10(Sigma2_1 * Poly_Branch),
    incr = log10(Sigma2_1 * Incr_Branch),
    decr = log10(Sigma2_1 * Decr_Branch)
  ) %>%
  select(Iteration, remain_mono, remain_olig, remain_poly, decr, incr) %>%
  pivot_longer(!Iteration, names_to = "part", values_to = "value") -> dat_PART

dat_PART$part %>%
  recode(
    remain_mono = "Monophagous",
    remain_olig = "Oligophagous",
    remain_poly = "Polyphagous",
    decr = "Decrease",
    incr = "Increase"
  ) -> dat_PART$part

dat_PART$part <- factor(
  dat_PART$part,
  levels = c("Polyphagous", "Oligophagous", "Monophagous", "Decrease", "Increase")
)

# graph
ggplot(dat_PART, aes(x = value)) +
  stat_halfeye(aes(fill = part, color = part),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  theme_classic(base_line_size = 0.3) +
  theme(
    axis.title = element_text(color = "black", size = 9),
    axis.text = element_text(color = "black", size = 8.5),
    legend.title = element_text(color = "black", size = 9),
    legend.text = element_text(color = "black", size = 8.5),
    legend.key.height = unit(3.6, "mm"),
    legend.key.width = unit(5, "mm"),
    legend.key.spacing.y = unit(1, "mm")
  ) +
  labs(
    y = "Density",
    x = expression("Log"[10] ~ sigma^2),
    fill = "Branch partition",
    color = "Branch partition"
  ) +
  scale_fill_brewer(
    palette = "Set2",
    breaks = c(
      "Increase",
      "Decrease",
      "Monophagous",
      "Oligophagous",
      "Polyphagous"
    )
  ) +
  scale_color_brewer(
    palette = "Set2",
    breaks = c(
      "Increase",
      "Decrease",
      "Monophagous",
      "Oligophagous",
      "Polyphagous"
    )
  ) -> pp1

# pMCMC: proportion of the posterior sample in which the difference in rate between
# partitions crossed 0 using a critical level of 0.05
part_BT %>%
  mutate(
    remain_mono = Sigma2_1,
    remain_olig = (Sigma2_1 * Olig_Branch),
    remain_poly = (Sigma2_1 * Poly_Branch),
    incr = (Sigma2_1 * Incr_Branch),
    decr = (Sigma2_1 * Decr_Branch)
  ) -> dd1

# pMCMC (two-tailed)
z <- mean(dd1$remain_poly > dd1$remain_mono)
min(z, 1 - z) * 2


## Three Partitions ----

# In this partition, branches are assigned according to the state of their ancestral nodes.

# read in BayesTraits output
part3_BT <- read_table("bt/output/res-3p-log.txt")

# convergence (ESS)
part3_BT |>
  select(-c(Iteration, Lh, Tree_No)) |>
  as.mcmc() |>
  effectiveSize()

# data for plotting
part3_BT %>%
  mutate(
    Polyphagous = log10(Sigma2_1),
    Monophagous = log10(Sigma2_1 * ancM_Branch),
    Oligophagous = log10(Sigma2_1 * ancO_Branch)
  ) %>%
  select(Iteration, Monophagous, Oligophagous, Polyphagous) %>%
  pivot_longer(!Iteration, names_to = "part", values_to = "value") -> dat_PART2

# dat_PART2$part <- factor(dat_PART2$part, levels = c("Increase", "Decrease", "Unchanged"))
# dat_PART2$part <- factor(dat_PART2$part, levels = c("Unchanged", "Decrease", "Increase"))

# graph
ggplot(dat_PART2, aes(x = value)) +
  stat_halfeye(aes(fill = part, color = part),
    alpha = 0.35,
    point_interval = "median_qi",
    point_alpha = 1,
    interval_alpha = 1,
    point_size = 2,
    position = position_dodgejust(),
    scale = 0.9
  ) +
  theme_classic(base_line_size = 0.3) +
  theme(
    axis.title = element_text(color = "black", size = 9),
    axis.text = element_text(color = "black", size = 8.5),
    legend.title = element_text(color = "black", size = 9),
    legend.text = element_text(color = "black", size = 8.5),
    legend.key.height = unit(3.6, "mm"),
    legend.key.width = unit(5, "mm"),
    legend.key.spacing.y = unit(1, "mm")
  ) +
  labs(
    y = "Density",
    x = expression("Log"[10] ~ sigma^2),
    fill = "Branch partition",
    color = "Branch partition"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") -> p1

# pMCMC tags
part3_BT %>%
  mutate(
    Polyphagous = (Sigma2_1),
    Monophagous = (Sigma2_1 * ancM_Branch),
    Oligophagous = (Sigma2_1 * ancO_Branch)
  ) -> dd1

# two-tailed pMCMC
z <- mean(dd1$Monophagous > dd1$Polyphagous)
min(z, 1 - z) * 2

p1 +
  annotate(geom = "text", label = "a", x = -3.97, y = -0.43, size = 3) +
  annotate(geom = "text", label = "ab", x = -3.8, y = -0.10, size = 3) +
  annotate(geom = "text", label = "b", x = -3.73, y = 0.24, size = 3) -> fig1b


## Posterior Summary ----

# Posterior summary, including ESS

# Five Partitions
# effective sample size
part_BT |>
  select(-c(Iteration, Lh, Tree_No)) |>
  as.mcmc() |>
  effectiveSize() -> ess_5p

# summary table
part_BT |>
  mutate(
    log10_rate_M = log10(Sigma2_1),
    log10_rate_O = log10(Sigma2_1 * Olig_Branch),
    log10_rate_P = log10(Sigma2_1 * Poly_Branch),
    log10_rate_I = log10(Sigma2_1 * Incr_Branch),
    log10_rate_D = log10(Sigma2_1 * Decr_Branch)
  ) |>
  select(
    Iteration,
    Alpha_1,
    log10_rate_M,
    log10_rate_O,
    log10_rate_P,
    log10_rate_I,
    log10_rate_D,
    No_RJ_Local_Branch,
    No_RJ_Local_Node
  ) |>
  pivot_longer(!Iteration, names_to = "part", values_to = "value") |>
  group_by(part) |>
  summarise(
    med = median(value, na.rm = TRUE),
    l95 = quantile(value, 0.025, na.rm = TRUE),
    u95 = quantile(value, 0.975, na.rm = TRUE)
  ) |>
  # parameter summary
  mutate(
    ESS = c(
      ess_5p["Alpha_1"],
      ess_5p["No_RJ_Local_Branch"],
      ess_5p["No_RJ_Local_Node"],
      ess_5p["Decr_Branch"],
      ess_5p["Incr_Branch"],
      ess_5p["Sigma2_1"],
      ess_5p["Olig_Branch"],
      ess_5p["Poly_Branch"]
    )
  ) |>
  rename("parameter" = "part") |>
  mutate(Branch_part = "Five partitions") -> summ_5p

# Three Partitions
# effective sample size
part3_BT |>
  select(-c(Iteration, Lh, Tree_No)) |>
  as.mcmc() |>
  effectiveSize() -> ess_3p

# summary table
part3_BT |>
  mutate(
    log10_rate_P = log10(Sigma2_1),
    log10_rate_M = log10(Sigma2_1 * ancM_Branch),
    log10_rate_O = log10(Sigma2_1 * ancO_Branch)
  ) |>
  select(
    Iteration,
    Alpha_1,
    log10_rate_M,
    log10_rate_O,
    log10_rate_P,
    No_RJ_Local_Branch,
    No_RJ_Local_Node
  ) |>
  pivot_longer(!Iteration, names_to = "part", values_to = "value") |>
  group_by(part) |>
  summarise(
    med = median(value, na.rm = TRUE),
    l95 = quantile(value, 0.025, na.rm = TRUE),
    u95 = quantile(value, 0.975, na.rm = TRUE)
  ) |>
  # parameter summary
  mutate(
    ESS = c(
      ess_3p["Alpha_1"],
      ess_3p["No_RJ_Local_Branch"],
      ess_3p["No_RJ_Local_Node"],
      ess_3p["ancM_Branch"],
      ess_3p["ancO_Branch"],
      ess_3p["Sigma2_1"]
    )
  ) |>
  rename("parameter" = "part") |>
  mutate(Branch_part = "Three partitions") -> summ_3p

rm(ess_3p, ess_5p)

# Combine partitions
bind_rows(summ_5p, summ_3p) |>
  select(Branch_part, parameter, med, l95, u95, ESS) |>
  flextable() |>
  theme_vanilla() |>
  merge_v(~Branch_part) |>
  colformat_double(digits = 2) |>
  colformat_double(j = "ESS", digits = 0) # |> save_as_docx(path = "out/bt-table.docx")
