# HPC workflow

The computationally intensive analyses were performed on a SLURM-based high-performance computing (HPC) environment. As HPC systems differ in available software modules, scheduler configuration, paths, and resource allocation policies, these scripts should be regarded as templates and will require adaptation before use on a different system.

The original analyses used:

-   SLURM workload manager
-   project-local R installation
-   renv-managed R package environment
-   JAGS for MCMC sampling

Users should review and adapt:

-   module load commands
-   account and partition settings
-   paths to R installations
-   resource requests (CPU, memory, runtime)
-   scheduler-specific directives
-   exact model formulation (in `hpc/scripts`)
-   seeds and MCMC settings (in `hpc/submit.sh`)

Original analysis environment used for the publication:

``` bash
# Load required HPC modules
module load USS/2022
module load gcc/9.4.0-pe5.34
module load libxt/1.1.5-pe5.34
module load libx11/1.7.0-pe5.34
module load readline/8.1
module load libiconv/1.16
module load zlib/1.2.12
module load bzip2/1.0.8
module load pcre2/10.39
module load curl/7.83.0-pe5.34
module load gettext/0.21
module load cairo/1.16.0-py3.9-pe5.34
module load pango/1.50.7-py3.9-pe5.34
module load freetype/2.11.1
module load libjpeg-turbo/2.1.3-pe5.34
module load libpng/1.6.37
module load libtiff/4.3.0
module load libtirpc/1.2.6-pe5.34
module load libdeflate/1.10
module load zstd/1.5.2
module load icu4c/67.1-py3.9-pe5.34
module load ghostscript/9.54.0-py3.9-pe5.34
module load jdk/19
module load cmake/3.23.1-pe5.34
module load udunits/2.2.28
module load proj/6.3.2-pe5.34
module load sqlite/3.38.5
module load gdal/3.0.4-py3.9-pe5.34
module load geos/3.9.1-py3.9-pe5.34
module load openblas/0.3.20-pe5.34
module load jags/4.3.2-openblas
module load gsl/2.7.1

# Additional environment variables were configured for package
# compilation (e.g. GDAL, GEOS, PROJ, GSL and UDUNITS libraries).
```
