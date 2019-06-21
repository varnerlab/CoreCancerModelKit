# Introduction
The ``CoreCancerModelKit`` package is a [Julia](https://julialang.org/downloads/) implementation
of the Core Cancer Metabolic Network published by Palsson and coworkers:

[Zielinski et al., (2017) Systems biology analysis of drivers underlying hallmarks of cancer cell metabolism. Sci Reports, 7:41241, doi: 10.1038/srep41241](https://rdcu.be/Olwc)

We've ported parts of the [COBRA](https://opencobra.github.io/cobratoolbox/stable/) toolbox to [Julia](https://julialang.org/downloads/) and added/reformulated the constraints in the flux balance analysis problem.

### Model edits
To edit the core model, we developed the [CBModelTools.jl](https://github.com/varnerlab/CBModelTools) package to extract the [MATLAB](https://www.mathworks.com/products/matlab.html) [COBRA](https://opencobra.github.io/cobratoolbox/stable/) binary file into a human readable/editable text format, updated this file with new reactions, and additional information such as [KEGG metabolite IDs](https://www.genome.jp/kegg/compound/), [EC numbers](https://en.wikipedia.org/wiki/Enzyme_Commission_number) for the reactions, etc and then re-encoded the updated version as a [COBRA](https://opencobra.github.io/cobratoolbox/stable/) binary file. To help with this, we also developed the [KEGG.jl](https://github.com/varnerlab/Kegg) package to automatically query the [KEGG API](https://www.kegg.jp/kegg/docs/keggapi.html) from [Julia](https://julialang.org/downloads/).   

## Installation and requirements
``CoreCancerModelKit.jl`` can be installed in the ``package mode`` of Julia.
Start of the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/index.html) and enter the ``package mode`` using the ``]`` key (to get back press the ``backspace`` or ``^C`` keys). Then, at the prompt enter:

    (v1.1) pkg> add https://github.com/varnerlab/CoreCancerModelKit.git

This will install the ``CoreCancerModelKit.jl`` package and other all required packages.
``CoreCancerModelKit.jl`` requires Julia 1.x and above.

## Table of contents
* [Example calculation scripts](/examples/README.md)

## Funding
The work described was supported by the [Center on the Physics of Cancer Metabolism at Cornell University](https://psoc.engineering.cornell.edu) through Award Number 1U54CA210184-01 from the [National Cancer Institute](https://www.cancer.gov). The content is solely the responsibility of the authors and does not necessarily
represent the official views of the [National Cancer Institute](https://www.cancer.gov) or the [National Institutes of Health](https://www.nih.gov).  
