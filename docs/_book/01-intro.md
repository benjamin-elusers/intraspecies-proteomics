# Proteomics experiment

## Experiment design

For this experiment, we selected eight strains from a population
of 1,011 *S. cerevisiae* isolates, representing the overall ecological, geographical
and genetic diversity (Figure \@ref(fig:phylo-8-strain). 

The eight strains were grown on synthetic defined media (SD). We carefully monitored
the strains growth to harvest cells close to the mid-log phase (OD ~ 0.5). Then,
we proceed to wash the samples in PBS and flash-freeze cell pellets for a 
whole-lysate proteomics profiling.

For each strains, we prepared  biological replicates that originated from two
distinct colonies which and two technical replicates.

<div class="figure" style="text-align: center">
<img src="images/phylogenetic-tree-8-strains.png" alt="Phylogenetic tree of 1011 *cerevisiae* isolates highlighting the 8 strains used for RNASeq/RiboSeq/proteomics exepriment" width="40%" />
<p class="caption">(\#fig:phylo-8-strain)Phylogenetic tree of 1011 *cerevisiae* isolates highlighting the 8 strains used for RNASeq/RiboSeq/proteomics exepriment</p>
</div>

Altogether there were 32 samples of yeast cell pellets including 4 replicates 
(two biological and two technical) that we sent to the Weizmann proteomics unit.
([INPCM](https://g-incpm.weizmann.ac.il/units/deBottonProteinProfiling/how-does-it-work%3F)).

The samples submitted correspond to 4mL cultures with OD ranging from 0.4-0.8. 

**OD of samples (table)**


## Sample preparation

The cell pellets were subjected to lysis and in solution tryptic digestion using
the S-Trap method (by Protifi) followed by a solid phase extraction cleaning step (Oasis HLB).

## Liquid chromatography mass spectrometry

The resulting peptides were analyzed using nanoflow liquid chromatography 
(nanoAcquity) coupled to high resolution, high mass accuracy mass spectrometry 
(Thermo Exploris 480). Each sample was analyzed on the instrument separately in 
a random order in discovery mode.

## Peptide identification and quantification

Raw data was processed with MaxQuant v1.6.6.0. The data were searched with the 
Andromeda search engine against a database containing protein sequences of 
*Saccharomyces cerevisiae* as downloaded from Uniprot.org, and appended with
common lab protein contaminants. 

The following modifications were defined for the search: 
Fixed modification- cysteine carbamidomethylation. 
Variable modifications- methionine oxidation and protein N-terminal acetylation. 

The quantitative comparisons were calculated using Perseus v1.6.0.7. 
Decoy hits were filtered out and only proteins that were detected
in at least two replicates of at least one experimental group were kept.
