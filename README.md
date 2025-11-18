# scrdiff_analysis
Simulation Analysis for Stationary Point Constrained Inference via Diffeomorphisms. 

## Files: 
*applied_data_analysis.R* - Applies method to ERP data. Data included in the package. 
*example.R* - Runs method on example dataset. Identical to that of the vignette from the *scrdiff* R package.
*images.R* - Code used create images in manuscript. 
*Sim_1.R* - Simulation study for the dataset with one stationary point. Records metrics for both the Diffeomorphism and DGP method. 
*Sim_2.R* - Simulation study for the dataset with two stationary points. Only runs the Diffeomorphism Method.



## Simulation Metrics:
We perform simulations for datasets of size 50, 100, 200 and 300. 
- For the first simulation, we record the MAP estimates, posterior means, and posterior standard deviations of the stationary point for both methods. The results are located in the *Sim_Results/Sim_1/Metrics* folder. We also record the 90, 95, and 99 HDP Intervals for the posterior of the stationary point, which are located in the *Sim_Results/Sim_1/HDP_Ints* folder. 
- For the second simulation, we record the MAP estimates, posterior means and marginal posterior standard deviations for each stationary point. The results are located in the *Sim_Results/Sim_2/Metrics*. The 90, 95, 99 and Bonferonni Corrected HDP Intervals for each stationary point was also recorded, and is located in *Sim_Results/Sim_2/HDP_Ints* folder. 