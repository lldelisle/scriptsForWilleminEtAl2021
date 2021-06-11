options(stringsAsFactors = F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("scales")
safelyLoadAPackageInCRANorBioconductor("ggpmisc") # To display equations

output <- commandArgs(TRUE)[1]

pretty.names <- c("Genome wide", "All TADs", "Btg1 TAD",
                  paste0("TAD", c(1:3, 6:8)))
names(pretty.names) <- c("genome_wide", "tads", "welcoming_tad",
                         paste0("tad_", c(1:3, 6:8)))
my.colors <- c("red", "black", hue_pal()(8)[c(8, 2:7)])
names(my.colors) <- pretty.names
decay.fn <- paste0("decay_", names(pretty.names), ".txt")
names(decay.fn) <- pretty.names
decay.df.list <- lapply(pretty.names, function(pn){
  temp.df <- read.delim(decay.fn[pn])
  temp.df$region <- pn
  return(temp.df)
})
decay.df <- do.call(rbind, decay.df.list)
decay.df$region <- factor(decay.df$region, levels = pretty.names)

# Formula used to plot
my.formula <- y ~ x

ggplot(decay.df, aes(x = Distance, y = Contacts, color = region)) +
  geom_line() +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_smooth(formula = my.formula, method = "lm", aes(weight = Number_bins),
              se = FALSE, lwd = 0.5, lty = 2) +
  stat_poly_eq(formula = my.formula, # Put the equations
               output.type = "numeric",
               aes(weight = Number_bins,
                   label = paste0("atop(NA, atop(alpha*\' = ", lapply(stat(coef.ls), function(v){round(v[[2, "Estimate"]], 2)}), "\',\'R\u00b2 = ", round(stat(r.squared), 2), "\'))")), # Use atop to have one on top of the other
               parse = TRUE,
               label.x = c(rep(0.2, 4),
                           rep(0.85, 5)), # Specify positions
               label.y = c(seq(0.5, 0.1, length.out = 4),
                           seq(0.95, 0.45, length.out = 5)),
               hjust = 0.5, vjust = 0.5,
               size = 6) +
  expand_limits(x=5e6) +
  scale_color_manual(values = my.colors) +
  theme_classic() +
  labs(color = "")

ggsave(output, width = 5, height = 5)

lm.welcoming <- lm(log10(Contacts) ~ log10(Distance), weight = Number_bins, subset(decay.df, region == "Btg1 TAD" & Distance > 0))
cat(lm.welcoming$coefficients[2], file = "powerlaw_decay_on_welcoming_weight.txt")
