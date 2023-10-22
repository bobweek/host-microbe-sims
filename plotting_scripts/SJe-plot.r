require("ggplot2")
require("ggdark")
require("ggthemes")
require("gridExtra")
require("reshape")
require("scales")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

clrs = gg_color_hue(5)

cmnthm_slides = theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'),
    plot.title = element_text(hjust = 0.5, size = 20),
    legend.position="none",
    axis.title = element_text(size=16))

cmnthm_ms = theme(
    axis.title = element_text(size=16))


fnl = read.csv("dat/SJₑ/fnl_rndpar_dat.csv")

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

ggplot() +
    geom_point(
        data=subdf,
        aes(S,dm),
        alpha=0.5) + 
    # geom_segment(
    #     aes(x=10,y=0,xend=1000,yend=0),
    #     color=clrs[4],
    #     lwd=1) +
    # geom_smooth(
    #     data=fnl,
    #     aes(Ne,dm),
    #     method=lm,
    #     fullrange=TRUE,
    #     # formula = y ~ x + I(x^2),
    #     color=clrs[1]) +
    coord_cartesian(ylim=c(-0.25,0.25)) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight() +
    xlab("Microbiome Richness, S") +
    ylab("Cohen's d") + 
    cmnthm_slides
    # cmnthm_ms


ggplot() +
    geom_point(
        data=subdf,
        aes(Je,dm),
        alpha=0.5) + 
    # geom_segment(
    #     aes(x=10,y=0,xend=1000,yend=0),
    #     color=clrs[4],
    #     lwd=1) +
    # geom_smooth(
    #     data=fnl,
    #     aes(Ne,dm),
    #     method=lm,
    #     fullrange=TRUE,
    #     # formula = y ~ x + I(x^2),
    #     color=clrs[1]) +
    coord_cartesian(ylim=c(-0.1,0.1)) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight() +
    xlab("Microbiome Size, Jₑ") +
    ylab("Cohen's d") + 
    cmnthm_slides
    