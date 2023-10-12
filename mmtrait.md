---
marp: true
theme: uncover
class:
    - lead
    - invert
paginate: false
# backgroundColor: #fff
# backgroundImage: url('https://marp.app/assets/hero-background.svg')
math: mathjax
---

<style>
{
   font-size: 0.7rem;
}
</style>

![bg left:40% 100%](metaorganism.jpg)

## **The Evolution of <br> Microbiome-Mediated Traits**

### Bob Week

#### 🐘 @bweek@ecoevo.social

Postdoc, University of Oregon
_Advised by:
Brendan Bohannan,
Peter Ralph, Bill Cresko_

---

## What is a Microbiome-Mediated Trait?

A quantitative character of a host organism with variance explained by variation of the host microbiome

* $P=\textcolor{orange}{G}+\textcolor{lightgreen}{M}$
* $\textcolor{lightgreen}{M}>0$

---

## Host  Trait Architecture
<br>

$$z=\textcolor{orange}{g}+\textcolor{lightgreen}{m}$$
<br>

* $\textcolor{orange}{g}\sim$ sum of additive effects across host loci

  * _additive_ genetic value
<br>

* $\textcolor{lightgreen}{m}\sim$ sum of additive effects across microbe species

  * _additive_ microbial value

---

## Host Life Cycle

![bg right:50% 80%](life-cycle.svg)

---

## Host Trait Evolution
<br>

![](trait-evo.svg)

* Want: $\Delta\bar z = \bar z'-\bar z$

* Need: $\bar m'$

---

## Microbial Inheritance

* $\textcolor{cyan}{\ell}=$ Lineal inheritance
  
* $1-\textcolor{cyan}{\ell}=$ Environmental acquisition

* $\textcolor{lightgreen}{\varepsilon}=$ Env-Env transmission

* $1-\textcolor{lightgreen}{\varepsilon}=$ Host shedding

* Collective inheritance:

  * $\textcolor{orange}{\kappa}=(1-\textcolor{cyan}{\ell})(1-\textcolor{lightgreen}{\varepsilon})$
<br>

* $\textcolor{cyan}{\ell}+\textcolor{orange}{\kappa}\leq1$

![bg right:45% 95%](acquisition-transmission.svg)

---

### Microbial Value Inheritance
<br>

* $\xi, \ \xi' =$ Environmental microbial values
<br>

* $m'=\textcolor{cyan}{\ell} m+(1-\textcolor{cyan}{\ell})(1-\textcolor{lightgreen}{\varepsilon})\bar m^*+(1-\textcolor{cyan}{\ell})\textcolor{lightgreen}{\varepsilon}\xi$
<br>

* $\bar m'=(\textcolor{cyan}{\ell}+\textcolor{orange}{\kappa})\bar m^*+(1-\textcolor{cyan}{\ell}-\textcolor{orange}{\kappa})\xi$

![bg right:40% 60%](mprime.svg)

---

### The Analytical Model
<br>

Assuming normally distributed $g,m,z$...
<br>

* $\Delta\bar z = \textcolor{orange}{G}\beta+\textcolor{lightgreen}{(\ell+\kappa)M}\beta+\textcolor{cyan}{(1-\ell-\kappa)(\xi-\bar m)}$
<br>

  * $\textcolor{orange}{G}=Var(\textcolor{orange}{g})$

  * $\beta=$ selection gradient

  * $\textcolor{lightgreen}{M}=Var(\textcolor{lightgreen}{m})$
<br>

* $\Delta\xi=\kappa(\textcolor{lightgreen}{M}\beta+\textcolor{cyan}{\bar m-\xi})/(1-\ell)$

---

<!-- ### Host Trait Evolution<br>Across Transmission Modes

$\small G=0, \ \ M=1, \ \ \beta=1, \ \ \bar z_0=\xi_0=0$
![](transmission-modes.svg)

--- -->

### Simulation Model
<br>

* Finite host population size $N_e$
<br>

* Explicit microbiome dynamics

  * $m=\sum_{i=1}^S\alpha_in_i$

  * $S \ =$ species richness

---

### Microbiome Dynamics

* Hubbell's neutral model

  * Fixed community size: $\sum_in_i=J_e$

  * $\vec n(t+1)\sim$ Multinomial$(J_e,\vec n(t))$

  * $L$ iterations per host generation

* Source of host trait variation

  * Microbiome analog to mutation

![bg left:40% 70%](micr-dyn.png)

---

### Simulation Setup

* $G_0=M_0=0$

* $L=30$

* $10$ replicates per parameter combination

  * Across $N_e, S, \ell, \kappa$

---

##### Mean Trait Evolution Across $N_e$

$\ell=1$

$\kappa=0$

$S=100$

![bg right:75% 90%](zp-N.svg)

---

###### Trait Variance Evolution Across $N_e$

$\ell=1$

$\kappa=0$

$S=100$

![bg right:75% 90%](Vp-N.svg)

---

![bg left:75% 90%](zp-S.svg)

##### Mean Trait Evolution Across $S$

$\ell=1$

$\kappa=0$

$N_e=1000$

---

![bg left:75% 90%](Vp-S.svg)

###### Trait Variance Evolution Across $S$

$\ell=1$

$\kappa=0$

$N_e=1000$

---

##### Mean Trait Evolution Across $\ell,\kappa$

$S=100$

$N_e=1000$

![bg right:75% 90%](zp-ℓ.svg)

---

###### Trait Variance Evolution Across $\ell,\kappa$

$S=100$

$N_e=1000$

![bg right:75% 90%](Vp-ℓ.svg)

---

### TODO

How well does the analytical model predict the simulation model over a single host generation?

---

### Fin