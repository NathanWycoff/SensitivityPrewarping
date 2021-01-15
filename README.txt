This respository contains R scripts used to produce each figure in the article. In each case, the script will assume that your working directory is the directory which contains this file. In each case, the output image may be found in the images directory after running the script from the terminal like so:

Figure 1) Rscript R/ridge_func.R
Figure 2) This illustration was created using GIMP. 
Figure 3) Rscript R/neighbor_change.R

Figures 4, 5 and 6 give our main numerical results.
The script R/bakeoff.R produces our numerical results, and accepts command line arguments specifying which test functions to try, as well as which prewarping methods.
After running this script with one of the below configurations, a file will be created in ""./data/runs/"" .
Run the concomitant file in R/plotfactory to produce the plot, making sure to modify the "load" statement to point to the newly created file in ""./data/runs/"" .

Figure 4) To run the comparisons:
- Rscript R/bakeoff.R -p de -d bcd
    To plot the results:
- Rscript R/plotfactory/cnc_temp_boxplots.R

Figure 5) To run the comparisons:
- Rscript R/bakeoff.R -p abc -d bcd
    To plot the results:
- Rscript R/plotfactory/lowDtest.R

Figure 6) Figure 6 involves a dataset of size 1.7G, which unfortunately cannot be uploaded to Manuscript Central. 

To create the associated pareto plots see "R/plotfactory/pareto.R" .

You will need to install the following R packages to use run these scripts:
install.packages(c("laGP", "lhs", "activegp", "hetGP", "FNN", "akima"))
