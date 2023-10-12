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


fnl = read.csv("Nₑ/fnl_rndpar_dat.csv")

fnl$abszd = abs(fnl$dz)
fnl$absgd = abs(fnl$dg)
fnl$absmd = abs(fnl$dm)

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

Nepl = ggplot() +
    geom_point(
        data=subdf,
        aes(Ne,dz),
        alpha=0.1) + 
    geom_segment(
        aes(x=10,y=0,xend=1000,yend=0),
        color=clrs[4],
        lwd=1) +
    geom_smooth(
        data=fnl,
        aes(Ne,dz),
        method=lm,
        fullrange=TRUE,
        # formula = y ~ x + I(x^2),
        color=clrs[1]) +
    coord_cartesian(ylim=c(-0.3,0.3)) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight() +
    xlab("Host Population Size, Nₑ") +
    ylab("Cohen's d") + 
    cmnthm_slides
    # cmnthm_ms



#
# S
#


fnl = read.csv("S/fnl_rndpar_dat.csv")

fnl$abszd = abs(fnl$dz)
fnl$absgd = abs(fnl$dg)
fnl$absmd = abs(fnl$dm)

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

Spl = ggplot() + 
    geom_point(
        data=subdf,
        aes(S,dz),
        alpha=0.1) +
    geom_segment(
        aes(x=10,y=0,xend=1000,yend=0),
        color=clrs[4],
        lwd=1) +
    geom_smooth(
        data=fnl,
        aes(S,dz),
        method=lm,
        # formula = y ~ x + I(x^2),
        color=clrs[1]) +
    coord_cartesian(ylim=c(-0.3,0.3)) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight() +
    xlab("Microbiome Richness, S") +
    ylab("") + 
    cmnthm_slides
    # cmnthm_ms



#
# κ
#


fnl = read.csv("κ/fnl_rndpar_dat.csv")

fnl$abszd = abs(fnl$dz)
fnl$absgd = abs(fnl$dg)
fnl$absmd = abs(fnl$dm)

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

κpl = ggplot() +
    geom_point(data=subdf,
        aes(κ,dz),
        alpha=0.1) + 
    geom_segment(
        aes(x=0,y=0,xend=1,yend=0),
        color=clrs[4],
        lwd=1) +
    geom_smooth(data=fnl,
        aes(κ,dz),
        method=lm,
        fullrange=TRUE,
        # formula = y ~ x + I(x^2),
        color=clrs[1]) +
    coord_cartesian(ylim=c(-0.3,0.3)) +
    theme_fivethirtyeight() +
    xlab("Collective Inheritance, κ") +
    ylab("Cohen's d") + 
    cmnthm_slides
    # cmnthm_ms



#
# ℓκ
#


fnl = read.csv("ℓκ/fnl_rndpar_dat.csv")

fnl$abszd = abs(fnl$dz)
fnl$absgd = abs(fnl$dg)
fnl$absmd = abs(fnl$dm)

subdf = fnl[sample(nrow(fnl), 10000),]

# subdf = fnl

κℓpl = ggplot() +
    geom_point(data=subdf,
        aes(κℓ,dz),
        alpha=0.1) +
    geom_segment(
        aes(x=0,y=0,xend=1,yend=0),
        color=clrs[4],
        lwd=1) +
    geom_smooth(data=fnl,
        aes(κℓ,dz),
        method=lm,
        fullrange=TRUE,
        # formula = y ~ x + I(x^2),
        color=clrs[1]) +
    coord_cartesian(ylim=c(-0.3,0.3)) +
    theme_fivethirtyeight() +
    xlab("Microbial Inheritance, κ+ℓ") +
    ylab("") + 
    cmnthm_slides
    # cmnthm_ms



#
# Corr
#


# Sfnl = read.csv("S/fnl_rndpar_dat.csv")
# Nfnl = read.csv("Nₑ/fnl_rndpar_dat.csv")
# κfnl = read.csv("κ/fnl_rndpar_dat.csv")

# fnl = rbind(Sfnl,Nfnl,κfnl)

# rm(Sfnl,Nfnl,κfnl)

# fnl$abszd = abs(fnl$dz)
# fnl$absgd = abs(fnl$dg)
# fnl$absmd = abs(fnl$dm)

# subdf = fnl[sample(nrow(fnl), 10000),]

# # subdf = fnl

# Corrpl = ggplot() +
#     geom_point(
#         data=subdf,
#         aes(Corr,dz),
#         alpha=0.1
#         ) +
#     geom_segment(
#         aes(x=-1,y=0,xend=1,yend=0),
#         color=clrs[4],
#         lwd=1) +
#     geom_smooth(
#         data=fnl,
#         aes(Corr,dz),
#         method=lm,
#         fullrange=TRUE,
#         color=clrs[1]
#         ) +
#     coord_cartesian(ylim=c(-0.1,0.1)) +
#     theme_fivethirtyeight() +
#     xlab("Corr(g,m)") +
#     ylab("") + 
#     cmnthm_slides
#     # cmnthm_ms

dpl = grid.arrange(Nepl,Spl,κpl,κℓpl,nrow=2)

ggsave("d.svg",dpl,width=12,height=12,bg='transparent')


ggsave("d.png",dpl,width=8,height=8)
