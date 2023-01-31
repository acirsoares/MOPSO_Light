# MOPSO_Light  (Multi-Objective Particle Swarm Optimization - Light)

1) The MOPSO_Light project implements the MOPSO_Light algorithm: a stochastic multi-objective optimization method. 

It is based on the Particle Swarm Optimization algorithm adapted for multi-objective problems.(Ref#01)

 	Main-Features (implemented in FORTRAN):

	+  Optimality Criteria: ε-dominância

	+  Global leadersheap: Fixed Leader (up 3 to 8 iterations). 

	+  Particle Velocity Calculation (Vi):  follows the original expression from Eberhart and Kennedy [REF#2] with the addition of the inertial factor $(ω)$ suggested by Chatterjee e Siarry [REF#3] to enhance the global assessment capability at the beginning and the local assessment capability at final iterations
	
	+ Independent Pareto Front-object to storarge the global leaders : 


# How to install

clone the git project : 

	git clone git@github.com:acirsoares/MOPSO_Light.git


# Example usage

change the directory :

	cd MOPSO_Light/Examples/Ex_01_MOPSOL

execute the make file :

	make

run the .a file :

	./test1.a

run the python script to plot the results in files MOPSOL_test1_$i.csv ($i = 01, 02, 03, ... , 10)

	python3 test1_2D_Plot.py 


# Change log

version_01 (initial version)  

# License and author info

! @licence GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)
! @date Oct 12, 2022
! @author Acir M. Soares Jr. <acir@ufsj.edu.br>


# Appendix 

The velocity "Vi" formula for the "K+1" iteration follows the original expression from Eberhart and Kennedy [REF#2] with the addition of the inertial factor $(ω)$ suggested by Chatterjee e Siarry [REF#3] to enhance the global assessment capability at the beginning and the local assessment capability at final iterations:

$Vi(k+1)=ω.Vi(k)+C_1.φ_1[Xi_pbest(k)−Xi_(k)] + C_2.φ_2[Xi_gbest(k)−Xi(k)]$ (1)

where $ ω = ω_min + (ω_max-ω_min)(maxIT-k)**q /(maxIT)**q $ ;
$(ω_min=0.4d0)$                   REF#[3] Chatterjee and Siarry
$(ω_max=0.9d0)$                   REF#[3] Chatterjee and Siarry
$(q=1.2d0)$                       REF#[3] Chatterjee and Siarry
$(C1=2.05d0)$                     REF#[2] Eberhart and Kennedy
$(C2=2.05d0)$                     REF#[2] Eberhart and Kennedy
$ φ_1 , φ_2 = random number       REF#[2] Eberhart and Kennedy
Xi_pbest and Xi_gbest correspond to the local and global leadership respectively.
____________________________________________________________________________________

# REFERENCES                                                                        

REF#[1] Proceeding Series of the Brazilian Society of Computational and Applied   
          Mathematics, v. 7, n. 1, 2020.                                          

REF#[2] R. Eberhart and J. Kennedy. A new optimizer using particle swarm theory.  
          In Proceedings of the Sixth International Symposium on Micro Machine    
          and Human Science, 39-43, 1995. DOI: 10.1109/MHS.1995.494215.           

REF#[3] A. Chatterjee and P. Siarry. Nonlinear inertia weight variation for       
          dynamic adaptation in particle swarm optimization. Computers &          
          operations research, 33:859-871, 2006. DOI:10.1016/j.cor.2004.08.012    

