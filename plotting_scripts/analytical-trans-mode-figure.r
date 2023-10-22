require("ggplot2")
require("ggdark")
require("ggthemes")
require("gridExtra")
require("reshape")

# read in data from analytical model
adat <- read.csv("dat/a_dat.csv")

# melt z and xi columns
no_melty = colnames(adat)[c(2:3,5:14)]
meltadat = melt(adat,no_melty)

#
# labelling for facet_grid parser
#

# label which is which (ab~abiot/biotic)
meltadat$ab <- factor(
    meltadat$variable,
    labels = c("z̄","ξ")
    )

# label ℓ values
meltadat$ell <- factor(
    meltadat$ℓ,
    labels = c("ℓ = 0","ℓ = 0.25","ℓ = 0.5","ℓ = 0.75","ℓ = 1.0")
    )

# plots different κ values as different lines
meltadat$κ <- as.character(adat$κ)

ap <- ggplot(meltadat) +
    geom_line(aes(y = value, x = t, color = κ)) +
    facet_grid(vars(ab), vars(ell)) +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Host Generation") +
    ylab("Phenotypic Units") +
    theme(axis.title = element_text(size=14, face="bold"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent',color='transparent'), 
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size=14),
        strip.background =element_rect(fill="transparent"),
        strip.text.y.right = element_text(angle = 0, size=16),
        strip.text.x.top = element_text(size=16))

ggsave("transmission-modes.png", ap, width = 9, height = 5, bg='transparent')
