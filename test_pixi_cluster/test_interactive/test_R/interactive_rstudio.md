# Interactive RStudio session from a pixi env

Steps to start an interactive RStudio session on the cluster, using a pixi environment.

1. Request an interactive node:

   ```bash
   srun -p interactive -J rstudio --cpus-per-task=8 --mem=128GB --pty bash
   ```

2. Once on the compute node, activate the pixi environment:

   ```bash
   pixi shell
   ```

3. Start RStudio:

   ```bash
   rstudio.sh start
   ```
