require("ggplot2")
require("ggdark")
require("ggthemes")
require("gridExtra")
require("reshape")
library("scales")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

clrs = gg_color_hue(5)

# 1 ~ ℓ = 0
# 2 ~ ℓ = 1/2
# 3 ~ ℓ = 1

par1 <- read.csv("dat/κ/par_001.csv")
par2 <- read.csv("dat/κ/par_005.csv")
par3 <- read.csv("dat/κ/par_009.csv")

pop1 <- read.csv("dat/κ/pop_dat_001.csv")
pop2 <- read.csv("dat/κ/pop_dat_005.csv")
pop3 <- read.csv("dat/κ/pop_dat_009.csv")

# pop1$corr = pop1$GM/sqrt(pop1$G*pop1$M)
# pop2$corr = pop2$GM/sqrt(pop2$G*pop2$M)
# pop3$corr = pop3$GM/sqrt(pop3$G*pop3$M)

minmin = min(pop1$m,pop2$m,pop3$m,
                pop1$g,pop2$g,pop3$g)
# zmin = min(c(pop1$z),c(pop2$z),c(pop3$z))
# minmin = min(mmin,zmin)
zmax = max(pop1$g,pop2$g,pop3$g,
            pop1$m,pop2$m,pop3$m)
Pmax = max(pop1$G,pop2$G,pop3$G,
            pop1$M,pop2$M,pop3$M)

pop1$k <- as.character(pop1$k)
pop2$k <- as.character(pop2$k)
pop3$k <- as.character(pop3$k)

cmnthm = theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'),
    plot.title = element_text(hjust = 0.5, size = 20),
    legend.position="none",
    axis.title = element_text(size=16, face="bold"),
    axis.title.y.left = element_text(angle = 0, vjust = 0.5, size=25))

zp1 = ggplot(pop1, aes(t,z)) +
    geom_line(,alpha=0.25,color=clrs[1]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,zmax) +
    ylab("ℓ = 0") +
    ggtitle("z̄") +
    theme(axis.title.x = element_blank()) +
    cmnthm
zp2 = ggplot(pop2, aes(t,z)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[3]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,zmax) +
    ggtitle("") +
    ylab("ℓ = 1/2") +
    theme(axis.title.x = element_blank()) +
    cmnthm
zp3 = ggplot(pop3, aes(t,z)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[4]) +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Host Generation") +
    ylim(minmin,zmax) +
    ggtitle("") +
    ylab("ℓ = 1") +
    cmnthm

gp1 = ggplot(pop1, aes(t,g)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[1]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,zmax) +
    ggtitle("ℓ = 1, κ = 0") +
    ylab("ḡ") +
    theme(axis.title.x = element_blank()) +
    cmnthm
gp2 = ggplot(pop2, aes(t,g)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[3]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,zmax) +
    ggtitle("ℓ = κ = 1/2") +
    ylab("") +
    theme(axis.title.x = element_blank()) +
    cmnthm
gp3 = ggplot(pop3, aes(t,g)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[4]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,zmax) +
    ylab("") +
    ggtitle("ℓ = 0, κ = 1") +
    theme(axis.title.x = element_blank()) +
    cmnthm

mp1 = ggplot(pop1, aes(t,m)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[1]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,zmax) +
    ylab("m̄") +
    ggtitle("") +
    xlab("Host Generation") +
    cmnthm
mp2 = ggplot(pop2, aes(t,m)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[3]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,zmax) +
    ylab("") +
    ggtitle("") +
    xlab("Host Generation") +
    cmnthm
mp3 = ggplot(pop3, aes(t,m)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[4]) +
    dark_mode(theme_fivethirtyeight()) +
    ylab("") +
    xlab("Host Generation") +
    ylim(minmin,zmax) +
    ggtitle("") +
    cmnthm

zp = grid.arrange(
    gp1,gp2,gp3,
    mp1,mp2,mp3,    
    ncol=3)

ggsave("zp-κ.svg",zp,width=12,height=9,bg='transparent')


Pp1 = ggplot(pop1, aes(t,P)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[1]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,Pmax) +
    ylab("ℓ = 0") +
    ggtitle("P") +
    theme(axis.title.x = element_blank()) +
    cmnthm
Pp2 = ggplot(pop2, aes(t,P)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[3]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,Pmax) +
    ggtitle("") +
    ylab("ℓ = κ = 1/2") +
    theme(axis.title.x = element_blank()) +
    cmnthm
Pp3 = ggplot(pop3, aes(t,P)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[4]) +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Host Generation") +
    ylim(minmin,Pmax) +
    ggtitle("") +
    ylab("ℓ = 1, κ = 0") +
    cmnthm

Gp1 = ggplot(pop1, aes(t,G)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[1]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,Pmax) +
    ggtitle("ℓ = 1, κ = 0") +
    ylab("G") +
    theme(axis.title.x = element_blank()) +
    cmnthm
Gp2 = ggplot(pop2, aes(t,G)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[3]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,Pmax) +
    ggtitle("ℓ = κ = 1/2") +
    theme(axis.title.x = element_blank()) +
    ylab("") +
    cmnthm
Gp3 = ggplot(pop3, aes(t,G)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[4]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,Pmax) +
    ylab("") +
    theme(axis.title.x = element_blank()) +
    ggtitle("ℓ = 0, κ = 1") +
    cmnthm

Mp1 = ggplot(pop1, aes(t,M)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[1]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,Pmax) +
    ylab("M") +
    ggtitle("") +
    theme(axis.title.x = element_blank()) +
    cmnthm
Mp2 = ggplot(pop2, aes(t,M)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[3]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(minmin,Pmax) +
    ylab("") +
    ggtitle("") +
    theme(axis.title.x = element_blank()) +
    cmnthm
Mp3 = ggplot(pop3, aes(t,M)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[4]) +
    dark_mode(theme_fivethirtyeight()) +
    ylab("") +
    theme(axis.title.x = element_blank()) +
    ylim(minmin,Pmax) +
    ggtitle("") +
    cmnthm

GMp1 = ggplot(pop1, aes(t,corr)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[1]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(-1,1) +
    ylab("Corr(g,m)") +
    ggtitle("") +
    xlab("Host Generation") +
    cmnthm +
    theme(axis.title.y.left = element_text(angle = 90, vjust = 1.0, size=25))
GMp2 = ggplot(pop2, aes(t,corr)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[3]) +
    dark_mode(theme_fivethirtyeight()) +
    ylim(-1,1) +
    ylab("") +
    ggtitle("") +
    xlab("Host Generation") +
    cmnthm
GMp3 = ggplot(pop3, aes(t,corr)) +
    geom_line(aes(group=k),alpha=0.25,color=clrs[4]) +
    dark_mode(theme_fivethirtyeight()) +
    xlab("Host Generation") +
    ylim(-1,1) +
    ylab("") +
    ggtitle("") +
    cmnthm

# ,
#         plot.margin = unit(c(0,0,0,0), "pt")

Vp = grid.arrange(
    Gp1,Gp2,Gp3,
    Mp1,Mp2,Mp3,
    GMp1,GMp2,GMp3,
    ncol=3)

ggsave("Vp-κ.svg",Vp,width=12,height=9,bg='transparent')
