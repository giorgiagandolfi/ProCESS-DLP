# ProCESS-DLP
Repository for testing and development of ProCESS for DLP+ sequencing.

## Set-up
In order to perform a DLP-like ProCESS simulation, you need to install a specific ProCESS version from a specific branch.

```r
devtools::install_github("caravagnalab/ProCESS", ref="on-the-fly_mutations")
```

Or either using the docker image with singularity:

```bash
singularity pull docker://giorgiagandolfi97/process_on_the_fly:v2
```

## Simulation types

### DLP-Sim1

Low number of cells ~240.
