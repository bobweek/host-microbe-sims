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

# ℓ

fnl = read.csv("ℓ/fnl_dat.csv")

fnl$abszd = abs(fnl$z.ᵈ)
fnl$absgd = abs(fnl$ḡᵈ)
fnl$absmd = abs(fnl$m.ᵈ)

ℓs = unique(fnl$ℓ)
subdf = data.frame()
for(i in 1:length(ℓs)){
    sbst = subset(fnl,ℓ==ℓs[i])
    sbst = sbst[sample(nrow(sbst), 1000),]
    subdf = rbind(subdf,sbst)
}

ggplot() + 
    geom_point(data=subdf,
        aes(ℓ,z.ᵈ),
        alpha=0.1) +
    geom_smooth(data=fnl,
        aes(ℓ,z.ᵈ),
        method=lm,
        fullrange=TRUE,
        formula = y ~ x + I(x^2),
        color=clrs[1]) +
    geom_abline(intercept=0,slope=1e-1) +
    coord_cartesian(ylim=c(-0.05,0.05)) +
    theme_fivethirtyeight()

ggplot() + 
    geom_point(
        data=subdf,
        aes(ℓ,abszd),
        alpha=0.1) +
    geom_smooth(
        data=fnl,
        aes(ℓ,abszd),
        method=lm,
        fullrange=TRUE,
        formula = y ~ x + I(x^2),
        color=clrs[1]) +
    coord_cartesian(ylim=c(0,0.1)) +
    theme_fivethirtyeight()

# Nₑ

fnl = read.csv("Nₑ/fnl_dat.csv")

fnl$abszd = abs(fnl$z.ᵈ)
fnl$absgd = abs(fnl$ḡᵈ)
fnl$absmd = abs(fnl$m.ᵈ)

s = unique(fnl$Nₑ)
subdf = data.frame()
for(i in 1:length(s)){
    sbst = subset(fnl,Nₑ==s[i])
    sbst = sbst[sample(nrow(sbst), 1000),]
    subdf = rbind(subdf,sbst)
}

ggplot() + 
    geom_point(
        data=subdf,
        aes(Nₑ,z.ᵈ),
        alpha=0.1) +
    geom_smooth(
        data=fnl,
        aes(Nₑ,z.ᵈ),
        method=lm,
        fullrange=TRUE,
        formula = y ~ x + I(x^2),
        color=clrs[1]) +
    coord_cartesian(ylim=c(-0.05,0.05)) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight()

ggplot() + 
    geom_point(
        data=subdf,
        aes(Nₑ,abszd),
        alpha=0.1) +
    geom_smooth(
        data=fnl,
        aes(Nₑ,abszd),
        method=lm,
        formula = y ~ x + I(x^2),
        color=clrs[1]) +
    coord_cartesian(ylim=c(0,0.1)) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight()


# S

fnl = read.csv("S/fnl_dat.csv")

fnl$abszd = abs(fnl$z.ᵈ)
fnl$absgd = abs(fnl$ḡᵈ)
fnl$absmd = abs(fnl$m.ᵈ)

s = unique(fnl$S)
subdf = data.frame()
for(i in 1:length(s)){
    sbst = subset(fnl,S==s[i])
    sbst = sbst[sample(nrow(sbst), 1000),]
    subdf = rbind(subdf,sbst)
}

ggplot() + 
    geom_point(
        data=subdf,
        aes(S,z.ᵈ),
        alpha=0.1) +
    geom_smooth(
        data=fnl,
        aes(S,z.ᵈ),
        method=lm,
        formula = y ~ x + I(x^2),
        color=clrs[1]) +
    ylim(-0.05,0.05) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight()

ggplot() + 
    geom_point(data=subdf,
        aes(S,abszd),
        alpha=0.1) +
    geom_smooth(data=fnl,
        aes(S,abszd),
        method=lm,
        color=clrs[1]) +
    scale_x_continuous(trans='log10') +
    theme_fivethirtyeight()

# Corr

Sfnl = read.csv("S/fnl_dat.csv")
Nfnl = read.csv("Nₑ/fnl_dat.csv")
ℓfnl = read.csv("ℓ/fnl_dat.csv")

fnl = rbind(Sfnl,Nfnl,ℓfnl)

rm(Sfnl,Nfnl,ℓfnl)

fnl$abszd = abs(fnl$z.ᵈ)
fnl$absgd = abs(fnl$ḡᵈ)
fnl$absmd = abs(fnl$m.ᵈ)

subdf = fnl[sample(nrow(fnl), 10000),]

ggplot() +
    geom_point(
        data=subdf,
        aes(Corr,z.ᵈ),
        alpha=0.1
        ) +
    geom_smooth(
        data=fnl,
        aes(Corr,z.ᵈ),
        method=lm,
        formula = y ~ x + I(x^2),
        color=clrs[1]
        ) +
    geom_smooth(
        data=fnl,
        aes(Corr,abszd),
        method=lm,
        formula = y ~ x + I(x^2),
        color=clrs[1]) +
    ylim(-0.05,0.05) +
    theme_fivethirtyeight()

ggplot() +
    geom_point(
        data=subdf,
        aes(Corr,abszd),
        alpha=0.1) +
    geom_smooth(
        data=fnl,
        aes(Corr,abszd),
        method=lm,
        formula = y ~ x + I(x^2),
        color=clrs[1]) +
    ylim(0,0.1) +
    theme_fivethirtyeight()

