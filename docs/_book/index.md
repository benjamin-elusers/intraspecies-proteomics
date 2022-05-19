--- 
title: "Intraspecies proteome expression"
project: "Integrating timescales in protein evolution from distant species to populations"
author: "Benjamin Dubreuil"
github: "benjamin-elusers"
mail: "benjamin.dubreuil@weizmann.ac.il"
created: "24/04/2022"
date: "19-May-2022 15:09 "
home: ""
site: bookdown::bookdown_site
documentclass: book
bibliography:
- book.bib
- packages.bib
description: |
  This book is a report for analyzing a proteomic experiment between 8 isolated
  strains from Saccharomyces cerevisiae.
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---


# Introduction

## About

This project results from a collaboration between the groups of E. Levy (Weizmmann, IL)
and J. Schacherer (CNRS/Univ. Strasbourg). 

## Abstract

Cells and their proteomes constantly evolve and adapt to survive in a wide
variety of ecological niches. While phenotypic diversity arises in populations
on relatively short timescales (1 to 100s thousand years), longer timescales (1
to 100s million years) are involved in divergence between species. The
signatures of both evolutionary processes are carved in each and every single
protein sequence through accumulation of mutations, thereby reshaping cell
machineries, including protein complexes, signaling pathways, and metabolic
pathways for example.

Interestingly, over long timescales, different proteins accumulate mutations at
markedly different rates. For example, orthologous proteins in S. cerevisiae and
S. pombe share *42%* (+/- 14%) sequence identity on average, but some orthologs
exhibit higher conservation. For example, actins share *89.7%* sequence identity,
illustrating that some proteins tolerate mutations more than others. One
biophysical property correlates with sequence divergence more than any other and
that is protein abundance. The more abundant a protein, the more conserved its
sequence. While the abundance-conservation correlation is well established,
its mechanistic origin is not well understood. Toxicity associated with
mutation-induced misfolding has been suggested as a possible origin, but recent
works, including ours (Dubreuil et al. 2019) show it is unlikely
the main driving mechanism.

To resolve mechanisms shaping protein evolution, we propose to integrate
analyses of sequence evolution across timescales. On the one hand, evolution
across distant species shows wide dynamic range of sequence conservation, and
reflects how entire proteomes diverge with time. On the other hand, it also
presents two important drawbacks. First, that functional information of one
species (e.g., protein abundance) is implicitly extrapolated to other species,
and such an assumption is frequently incorrect. Second, the order in which
mutations occur can hardly be traced. For example, considering a gene
duplication event that occurred several million years ago, it is impossible to
distinguish mutations that occurred right after the duplication event from those
that arose later. Remarkably, both of these limitations are resolved when
measuring evolution across a population. Indeed, the short evolutionary
distances seen across strains mean that functional information can be
extrapolated with high confidence, and mutations can be situated in a narrow
time-window. Thus, by resolving these limitations, the integration of
evolutionary data across species and strains (Figure \@ref(fig:project-goal))
will provide new insights into constraints that biophysics (in particular
abundance) versus function place on protein evolution.

<div class="figure">
<img src="images/Fungi-SpeciesTree-v7.png" alt="Integrating timescales in protein evolution from distant species to populations" width="100%" />
<p class="caption">(\#fig:project-goal)Integrating timescales in protein evolution from distant species to populations</p>
</div>


