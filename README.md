# energy-allocation
## The environment selects: Modeling energy allocation in microbial communities under dynamic environments 

Leonor Guedes da Silva*, Sergio Tomas-Martinez*, Mark C.M. van Loosdrecht, S. Aljoscha Wahl

(*) Equal contribution


Department of Biotechnology, Delft University of Technology, The Netherlands

Corresponding Author: Leonor Guedes da Silva, LeonorGuedesdaSilva@gmail.com

### This model was adapted from the following previous work:
Elucidating temporal resource allocation and diurnal dynamics in phototrophic metabolism using conditional FBA

Marco Rugen, Alexander Bockmayr, and Ralf Steuerb

Sci Rep. 2015; 5: 15247.

Published online 2015 Oct 26. doi: 10.1038/srep15247

PMCID: PMC4620596

PMID: 26496972

#### Download original model here:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4620596/bin/srep15247-s2.zip

## How to...

1. Download or Clone "energy-allocation" into your computer
2. Open MATLAB (this was tested in MATLAB 2018a for macOS Mojave)
    >> load PAO.mat
    >> run_cFBA(model)
3. Once the simulation is done, results can be checked with:
    >> load results/res_cFBA_model_PAO.mat      %(be careful: '/' and '\' depends on OS)
    >> plotFlux(res)
    >> plotFlux(res, "PHB_S")
    >> plotCompound(res)
    >> plotCompound(res, "PHB")