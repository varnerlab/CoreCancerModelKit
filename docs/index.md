## Introduction
The ``CoreCancerModelKit`` package is a [Julia](https://julialang.org/downloads/) implementation
of the Core Cancer Metabolic Network published by Palsson and coworkers:

[Zielinski et al., (2017) Systems biology analysis of drivers underlying hallmarks of cancer cell metabolism. Sci Reports, 7:41241, doi: 10.1038/srep41241](https://rdcu.be/Olwc)

We've ported the parts of the [COBRA](https://opencobra.github.io/cobratoolbox/stable/) toolbox to [Julia](https://julialang.org/downloads/), added a few constraints, reformulated some of the constraints and added additional
biological data to the flux estimation problem.

### Requirements
To use the ``CoreCancerModelKit``package requires [Julia](https://julialang.org/downloads/) v1.1 or above.
The ``CoreCancerModelKit``package requires several additional packages, which are typically downloaded/installed automatically when ``CoreCancerModelKit`` is installed. See the [installation instructions](/installation/README.md)
for details.

### Table of contents
* [Installation](/installation/README.md)
* [Example calculation scripts](/examples/README.md)

## Funding
The work described was supported by the [Center on the Physics of Cancer Metabolism at Cornell University](https://psoc.engineering.cornell.edu) through Award Number 1U54CA210184-01 from the [National Cancer Institute](https://www.cancer.gov). The content is solely the responsibility of the authors and does not necessarily
represent the official views of the [National Cancer Institute](https://www.cancer.gov) or the [National Institutes of Health](https://www.nih.gov).  
