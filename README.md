# SingleMoleculeDataAnalysis
A repository for code related to the Global Fit data analysis used and produced in the laboratory of Julie S. Biteen at the University of Michigan and described in Rowland and Biteen, "Measuring molecular motions inside single cells with improved analysis of
single-particle trajectories," Chemical Physics Letters (2017) http://dx.doi.org/10.1016/j.cplett.2017.02.05

Code included:

CPDGlobal.m : Matlab program that uses tracking files (see header for data types) to estimate diffusion parameters via fitting to empirical cumulative probability distributions of squared step size.

cpdFunFinder.m : Matlab program that is auxiliary to CPDGlobal.m. This code designs the particular fitting function to be used and can be edited to include diffusion models of one's choice.

simpleDiffusion.m : Matlab program that generates trajectories of diffusing particles. The output has the data type that is expected by CPDGlobal.m.
