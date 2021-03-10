# Simulations

This folder includes a bash script for simulating small-to-large TR expansions in the FMR1 CGG motif, and a R script for plotting results. Simulations are performed using [VISOR](https://github.com/davidebolo1993/VISOR). To run the simulate.sh bash script, one must first replace the included path to the clusteranalysis python with a proper one.

``` bash
sed -i 's#/home/davide/tools/TREADMILL/scripts/clusteranalysis.py#/path/to/TREADMILL/scripts/clusteranalysis.py#g' simulate.sh
```

Simulations can then be run and results can be plotted with the included R script.

``` bash
bash simulate.sh
Rscript plotsim.R overview.tsv tpr.tsv
```
