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



#
# Nₑ
#


fnl = read.csv("dat/posNₑ/fnl_rndpar_dat.csv")

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

Nepl = ggplot() +
    geom_point(
        data=subdf,
        aes(Ne,dm),
        alpha=0.05) + 
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
    dark_mode(theme_fivethirtyeight()) +
    xlab("Host Population Size, Nₑ") +
    ylab("Cohen's d") + 
    cmnthm_slides
    # cmnthm_ms



#
# S
#


fnl = read.csv("dat/posS/fnl_rndpar_dat.csv")

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

Spl = ggplot() + 
    geom_point(
        data=subdf,
        aes(S,dm),
        alpha=0.05) +
    # geom_segment(
    #     aes(x=10,y=0,xend=1000,yend=0),
    #     color=clrs[4],
    #     lwd=1) +
    # geom_smooth(
    #     data=fnl,
    #     aes(S,dm),
    #     method=lm,
    #     # formula = y ~ x + I(x^2),
    #     color=clrs[1]) +
    coord_cartesian(ylim=c(-0.25,0.25)) +
    scale_x_continuous(trans='log10') +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Microbiome Richness, S") +
    ylab("") + 
    cmnthm_slides
    # cmnthm_ms



#
# κ
#


fnl = read.csv("dat/posκ/fnl_rndpar_dat.csv")

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

κpl = ggplot() +
    geom_point(data=subdf,
        aes(κ,dm),
        alpha=0.05) + 
    # geom_segment(
    #     aes(x=0,y=0,xend=1,yend=0),
    #     color=clrs[4],
    #     lwd=1) +
    # geom_smooth(data=fnl,
    #     aes(κ,dm),
    #     method=lm,
    #     fullrange=TRUE,
    #     # formula = y ~ x + I(x^2),
    #     color=clrs[1]) +
    coord_cartesian(ylim=c(-0.25,0.25)) +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Collective Inheritance, κ\nℓ=1-κ") +
    ylab("Cohen's d") + 
    cmnthm_slides
    # cmnthm_ms


#
# ℓκ
#


fnl = read.csv("dat/posκℓ/fnl_rndpar_dat.csv")

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

κℓpl = ggplot() +
    geom_point(data=subdf,
        aes(κℓ,dm),
        alpha=0.05) +
    # geom_segment(
    #     aes(x=0,y=0,xend=1,yend=0),
    #     color=clrs[4],
    #     lwd=1) +
    # geom_smooth(data=fnl,
    #     aes(κℓ,dm),
    #     method=lm,
    #     fullrange=TRUE,
    #     # formula = y ~ x + I(x^2),
    #     color=clrs[1]) +
    coord_cartesian(ylim=c(-0.25,0.25)) +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Microbial Inheritance, κ+ℓ\n") +
    ylab("") + 
    cmnthm_slides
    # cmnthm_ms

dpl = grid.arrange(Nepl,Spl,κpl,κℓpl,nrow=2);

ggsave("d_dm.svg",dpl,width=12,height=6,bg='transparent')
ggsave("d_dm.png",dpl,width=8,height=8)

#
# β
#

Sfnl = read.csv("dat/posS/fnl_rndpar_dat.csv")
Nfnl = read.csv("dat/posNₑ/fnl_rndpar_dat.csv")
κfnl = read.csv("dat/posκ/fnl_rndpar_dat.csv")
κℓfnl = read.csv("dat/posκℓ/fnl_rndpar_dat.csv")

fnl = rbind(Nfnl,Sfnl,κfnl,κℓfnl)

rm(Sfnl,Nfnl,κfnl,κℓfnl)

subdf = fnl[sample(nrow(fnl), 1000),]

# subdf = fnl

ggplot() +
    # geom_density2d(data=subdf,aes(β,dm)) +
    geom_point(data=subdf,
        aes(β,dm),
        alpha=0.05) +
    # geom_smooth(data=fnl,
    #     aes(β,dm),
    #     method=lm,
    #     # fullrange=TRUE,
    #     # formula = y ~ x + I(x^2),
    #     color=clrs[1]) +
    # coord_cartesian(
    #     xlim=c(-0.0001,0.0001),
    #     ylim=c(-0.25,0.25)) +
    theme_fivethirtyeight() +
    xlab("Selection Gradient, β") +
    ylab("") + 
    cmnthm_slides
# cmnthm_ms
head(fnl)
ggplot() +
    # geom_density2d(data=subdf,aes(β,dm)) +
    geom_point(data=subdf,
        aes(β,dm),
        alpha=0.05) +
    # geom_smooth(data=fnl,
    #     aes(β,dm),
    #     method=lm,
    #     # fullrange=TRUE,
    #     # formula = y ~ x + I(x^2),
    #     color=clrs[1]) +
    # coord_cartesian(
    #     xlim=c(-0.0001,0.0001),
    #     ylim=c(-0.25,0.25)) +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Selection Gradient, β") +
    ylab("") + 
    cmnthm_slides
# cmnthm_ms