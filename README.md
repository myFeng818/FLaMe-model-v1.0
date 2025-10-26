The Fluxes of Lake Methane (FLaMe-v1.0) model is written in the programming language Fortran (For) and runs well on PC or a clustering using a gfortran Compliler. 
The use of the cross platform Visual Studio Code, which can be downloaded for free from https://code.visualstudio.com/download, is recommended by the authors. 
A minimum random access memory (RAM) of 500 MB is required. 
A standard model simulation with 10 years can be finished in about several minutes on a regular PC, and the tests were performed successfully with Windows, Mac, and Linux operating systems. 
If users have problems or difficulties in compiling and running the model, please feel free to contact us for help (Maoyuan.feng@ulb.be or maoyuanfeng93@gmail.com).

1. A package named “FLaMe”, containing the source code of FLaMe as well as the model project, is available as supplementary material. The package contains 2 sub-folders:
(1) Subfolder named “INPUT_OUTPUT”, which contains two subfolders of “input” and “output” for the model inputs and outputs, respectively.
(2) Subfolder named “CODE”, which contains 36 subroutines that are classified into two groups:
	A) Subroutines inherited from Canadian Small Lake Model (CSLM):
- CLASSI.f
- CLASSL.f
- density.f
- DIASURFZ.f
- DRCOEF.f
- DRCOEFL.f
- EQNST.f
- FLXSURFZ.f
- FREECONV.f
- LKTRANS.f
- MIXLYR.f
- RUNLAKE.f (main for a single run)
- SCREENRH.f
- SLDIAG.f
- SNINFL.f
- SNOADD.f
- SNOALBA.f
- SNOALBW.f
- SNOVAP.f
- TLSPOST.f
- TLSPREP.f
- TMELT.f
- TSOLVE.f
- TSOLVL.f
- XIT.f

(B) Biogeochemical modules, newly developed in this study:
- area_layer.f
- CH4_module.f
- Clabile_module.f
- geometry_bgc.f
- initialisation.f
- O2_module.f
- oxidation.f
- search_zebmin.f
- transport.f
- Writing.f

2. Model compilation and run

	(1) Before compiling and running the FLaMe model, the users need to modify the parameters in RUNLAKE.f, which is the main of the model.
	The parameters need to be modified:
	GEN_DIR and FILE_DIR: parameters for the directory where the codes are deposited.
    SET: the scenario the users may use (the users may have multiple scenarios to test the model, so this parameter is used to set the scenario directory)
    x_per_m: the number of printing the output per month, we have two choice for this parameter. i.e., 1 for monthly output and 4 for weekly output.
    NLAT: Number of grid cells
    NMOS: Number of representative lakes in each grid cell

	(2) In the “CODE” subfolder, we provide a makefile to compile the FLaMe model, with the following commands:
    rm *.o: remove all the previously compiled files
    make -f Makefile_check: Compile the FLaMe model with a gfortran compiler
    ./RUNLAKE_check: Run the model, and you will find the new results in the “output” subfolder.


