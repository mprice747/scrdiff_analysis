# scrdiff_analysis
Simulation Analysis for Stationary Point Constrained Inference via Diffeomorphisms. Each file can be run line by line.

## Files: 
**applied_data_analysis.R** - Applies method to ERP data. Data included in the package. <br>
**example.R** - Runs method on example dataset. Identical to that of the vignette from the *scrdiff* R package. <br>
**images.R** - Code used create images in manuscript. <br>
**Sim_1.R** - Simulation study for the dataset with one stationary point. Records metrics for both the Diffeomorphism and DGP method. <br>
**Sim_2.R** - Simulation study for the dataset with two stationary points. Only runs the Diffeomorphism Method. <br>

## Data: 
All data used is located in the *Data* folder. The ERP data is located in *Data/Applied_Data_Analysis* folder, the data simulated from the model with one stationary point is located in *Data/Sim_1* and the data simulated from the model with two stationary points is located in *Data/Sim_2*. For the simulations, the "X" points are located in the files titled *X_matrix* and the "Y" points are located in the files titled *Y_matrix*. Finally, the number at the end of these files indicate the sample size. 



## Simulation Metrics:
We perform simulations for datasets of size 50, 100, 200 and 300. 
- For the first simulation, we record the MAP estimates, posterior means, and posterior standard deviations of the stationary point for both methods. The results are located in the *Sim_Results/Sim_1/Metrics* folder. We also record the 90, 95, and 99 HDP Intervals for the posterior of the stationary point, which are located in the *Sim_Results/Sim_1/HDP_Ints* folder. 
- For the second simulation, we record the MAP estimates, posterior means and marginal posterior standard deviations for each stationary point. The results are located in the *Sim_Results/Sim_2/Metrics*. The 90, 95, 99 and Bonferonni Corrected HDP Intervals for each stationary point were also recorded, and is located in *Sim_Results/Sim_2/HDP_Ints* folder. 