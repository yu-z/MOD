# MOD

## What is MOD?

* MOD stands for Graphics Processing Units Molecular Dynamics. It is a new molecular dynamics (MD) code implemented fully on graphics processing units (GPUs). This code is highly efficient. For details, see Ref. [1].

* Force evaluation for many-body potentials has been significantly accelerated by using GPUs. Our efficient and flexible GPU implementation of the force evaluation for many-body potentials relies on a set of simple expressions for force, virial stress, and heat current derived in Ref. [2]. Detailed algorithms for efficient CUDA-implementation have been presented in Ref. [1]. We have implemented the following many-body potentials in GPUMD:
   * The EAM-type potential with some analytical forms
   * The Tersoff (1989) potential with single or double atom types
   * The REBO potential for Mo-S systems (2009)
   * The Stillinger-Weber (1985) potential
   
* Apart from being highly efficient, another unique feature of GPUMD is that it has useful utilities to study heat transport. The current version of GPUMD can calculate the following quantities related to heat transport:
   * It can calculate the phonon density of states (DOS) from the velocity autocorrelation (VAC).
   * It can calculate the equilibrium heat current autocorrelation (HAC), whose time integral gives the running thermal conductivity   according to the Green-Kubo relation. As stressed in Ref. [2], the heat current as implemented in LAMMPS does not apply to many-body  potentials and significantly underestimates the thermal conductivity in 2D materials described by many-body potentials. GPUMD also contains the thermal conductivity decomposition method as introduced in Ref. [3], which is essential for 2D materials.
   * It can calculate the thermal conductivity of a system of finite length or the thermal boundary resistance (Kapitza resistance) of an interface or similar structures using nonequilibrium MD (NEMD) methods. The spectral decomposition method as described in Ref. [3] has also been implemented.
     
* GPUMD was firstly used for heat transport simulations only but we are now making it more and more general. However, the functionalities in GPUMD are still very limited. We are working on implementing (1) more potential models (including mixed potentials), (2) more integrators (including external conditions), and (3) more measurements.

* Note that we will not implement any function for building simulation models. Users of GPUMD are supposed to be able to build simulation models by their own.

## Compile GPUMD

* To compile GPUMD, one just needs to go to the "src" directory and type "make". one may want to first do "make clean". When the compilation finishes, an executable named "gpumd" will be generated in the "src" directory. See the manual in the "doc" directory for details.

## Run GPUMD

* Take the first example in the "examples" directory for example, one can use the following command to run it:
  * src/gpumd < examples/input1.txt
* Please read the manual to study the examples. These examples should get you started. 
  
## Authors:

* Yu Zhang (Aalto University)
  * zhangy9@sustc.edu.cn
* Ville Vierimaa (Aalto University)  
* Mikko Ervasti (Aalto University)  
* Ari Harju (Aalto University) 

## Feedbacks:

There is a comprehensive manual in the "doc" directory. You can e-mail the first author (Zheyong Fan) if you find errors in the manual or bugs in the source code, or have any suggestions/questions about the manual and code. Thank you!

## Citations

If you use MOD in your published work, we kindly ask you to cite the following paper which describes the central algorithms used in MOD:
* [1] xxx


