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
