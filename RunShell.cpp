/*
 *  RunShell.cpp
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 7/27/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

/*


 */
#include <filesystem>
namespace fs = std::filesystem;


#include "Main.H"
#include "Vector.H"
#include <cstring>
#include <string>
#include "Errors.H"
#include "Matrix.H"
#include "MatlabFileHandle.H"
#include "BinaryFileHandle.H"
#include "TextFileHandle.H"
// #include "RefCountedPointer.H"  // obsolete RKLibrary header
#include "SpecialFunctions.H"
// #include "LapackWrapper.H"      // obsolete: LU routines are defined in TinyMatrix.H
#include "Tensor2D.H"
#include "TinyVector.H"
#include "TinyMatrix.H"
#include "NonEuclideanShell.H"
#include "mathexpr.h"

#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>



/* forward function declarations (for optimizer) */
double my_f  (const gsl_vector *, void *);
void   my_df (const gsl_vector *, void *, gsl_vector *);
void   my_fdf(const gsl_vector *, void *, double *, gsl_vector *);

/* Forward declaration of input functions */
TinyMatrix<double,2> inputFunctionAbar     (double, double);
TinyMatrix<double,2> inputFunctionBbar     (double, double);
TinyMatrix<double,2> inputFunctionCbar     (double, double);
TinyVector<double,3> inputFunctionPos0     (double, double);
double               inputFunctionThickness(double, double);
double               inputFunctionLambda   (double, double);
double               inputFunctionMu       (double, double);

/* Global strings of the formulas */
char* abar11_formula;
char* abar12_formula;
char* abar21_formula;
char* abar22_formula;
char* bbar11_formula;
char* bbar12_formula;
char* bbar21_formula;
char* bbar22_formula;

char* cbar11_formula;
char* cbar12_formula;
char* cbar21_formula;
char* cbar22_formula;

char* thickness_formula;
char* Y_formula;
char* nu_formula;
char* X0_formula;
char* Y0_formula;
char* Z0_formula;
// ────────────────────────────────────────────────────────────
// Target‐connection Γ̄ᵏᵢⱼ formulas (k,i,j ∈ {1,2})
char* gammabar111_formula;
char* gammabar112_formula;
char* gammabar121_formula;
char* gammabar122_formula;
char* gammabar211_formula;
char* gammabar212_formula;
char* gammabar221_formula;
char* gammabar222_formula;
// ----------------------------------------------------------------------------
 TinyMatrix<TinyMatrix<double,2>,2> inputFunctionGammaBar(double, double);

// Build Γ̄ᵏᵢⱼ from the eight user strings:

TinyMatrix<TinyMatrix<double,2>,2>
inputFunctionGammaBar(double u, double v)
{
    RVar uvar("u",&u), vvar("v",&v);
    RVar* vars[2] = { &uvar, &vvar };

    ROperation op000(gammabar111_formula,2,vars),
               op001(gammabar112_formula,2,vars),
               op010(gammabar121_formula,2,vars),
               op011(gammabar122_formula,2,vars),
               op100(gammabar211_formula,2,vars),
               op101(gammabar212_formula,2,vars),
               op110(gammabar221_formula,2,vars),
               op111(gammabar222_formula,2,vars);

    TinyMatrix<TinyMatrix<double,2>,2> G;
    G(0,0)(0,0) = op000.Val();  G(0,0)(0,1) = op001.Val();
    G(0,1)(0,0) = op010.Val();  G(0,1)(0,1) = op011.Val();
    G(1,0)(0,0) = op100.Val();  G(1,0)(0,1) = op101.Val();
    G(1,1)(0,0) = op110.Val();  G(1,1)(0,1) = op111.Val();
    return G;
}

int main() {

	/* Declare parameter */
	int NumberOfLoops      = 1000;
	int HowOftenToPrint    = 5000;
	int restart            = 1;
	double ThicknessAdjust = 1.0;
	double MetricAdjust    = 1.0;
	double EtaAdjust = 1.0;
	


	std::string  nodeOutputFileName;
	std::string  faceOutputFileName;
	std::string  EFGLMNOutputFileName;
	std::string  restartFileName;
	std::string  energyFileName;
	std::string  saveFileName;
	std::string  verticesFileName;
	std::string  facesFileName;
	std::string  inputFormula;

	// std::cout << "Compiled 13/06/2025" << std::endl;
	/* input parameters */
	// std::cout << "Enter the vertices file name " << std::endl;
	std::cin >> verticesFileName;	// std::cout << "  " << verticesFileName << std::endl;
	// std::cout << "Enter the faces file name " << std::endl;
	std::cin >> facesFileName;		// std::cout << "  " << facesFileName << std::endl;
	// std::cout << "Enter Number of loops (time steps or optimization)" << std::endl;
	std::cin >> NumberOfLoops;		// std::cout << "  " << NumberOfLoops << std::endl;
	// std::cout << "Enter how often to save state " << std::endl;
	std::cin >> HowOftenToPrint;	// std::cout << "  " << HowOftenToPrint << std::endl;

	// std::cout << "Enter formulas in (u,v) for the elements of the target metric " << std::endl;
	std::cin >> inputFormula;
	abar11_formula = new char[inputFormula.length() + 1];
	strcpy(abar11_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	abar12_formula = new char[inputFormula.length() + 1];
	strcpy(abar12_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	abar21_formula = new char[inputFormula.length() + 1];
	strcpy(abar21_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	abar22_formula = new char[inputFormula.length() + 1];
	strcpy(abar22_formula, inputFormula.c_str());
	// std::cout << "  (  " << abar11_formula << "  ,  "  << abar12_formula << "  )  "  << std::endl;
	// std::cout << "  (  " << abar21_formula << "  ,  "  << abar22_formula << "  )  "  << std::endl;
	// std::cout << "Enter formulas in (u,v) for the elements of the target curvature " << std::endl;
	std::cin >> inputFormula;
	bbar11_formula = new char[inputFormula.length() + 1];
	strcpy(bbar11_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	bbar12_formula = new char[inputFormula.length() + 1];
	strcpy(bbar12_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	bbar21_formula = new char[inputFormula.length() + 1];
	strcpy(bbar21_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	bbar22_formula = new char[inputFormula.length() + 1];
	strcpy(bbar22_formula, inputFormula.c_str());
	// std::cout << "  (  " << bbar11_formula << "  ,  "  << bbar12_formula << "  )  "  << std::endl;
	// std::cout << "  (  " << bbar21_formula << "  ,  "  << bbar22_formula << "  )  "  << std::endl;
	// Read formulas for the reference third fundamental form cbar
	std::cin >> inputFormula;
	cbar11_formula = new char[inputFormula.length() + 1];
	strcpy(cbar11_formula, inputFormula.c_str());

	std::cin >> inputFormula;
	cbar12_formula = new char[inputFormula.length() + 1];
	strcpy(cbar12_formula, inputFormula.c_str());

	std::cin >> inputFormula;
	cbar21_formula = new char[inputFormula.length() + 1];
	strcpy(cbar21_formula, inputFormula.c_str());

	std::cin >> inputFormula;
	cbar22_formula = new char[inputFormula.length() + 1];
	strcpy(cbar22_formula, inputFormula.c_str());
	// std::cout << "Enter formula in (u,v) for the thickness " << std::endl;
	std::cin >> inputFormula;
	thickness_formula = new char[inputFormula.length() + 1];
	strcpy(thickness_formula, inputFormula.c_str());
	std::cout << "[DEBUG] thickness_formula = \"" << thickness_formula << "\"" << std::endl;
	// std::cout << "  " << thickness_formula << std::endl;
	// std::cout << "Enter formula in (u,v) for Young's Modulus " << std::endl;
	std::cin >> inputFormula;
	Y_formula = new char[inputFormula.length() + 1];
	strcpy(Y_formula, inputFormula.c_str());
	// std::cout << "  " << Y_formula << std::endl;
	// std::cout << "Enter formula in (u,v) for Poisson's ratio " << std::endl;
	std::cin >> inputFormula;
	nu_formula = new char[inputFormula.length() + 1];
	strcpy(nu_formula, inputFormula.c_str());
	// std::cout << "  " << nu_formula << std::endl;

	// std::cout << "Enter formulas in (u,v) for the initial state " << std::endl;
	std::cin >> inputFormula;
	X0_formula = new char[inputFormula.length() + 1];
	strcpy(X0_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	Y0_formula = new char[inputFormula.length() + 1];
	strcpy(Y0_formula, inputFormula.c_str());
	std::cin >> inputFormula;
	Z0_formula = new char[inputFormula.length() + 1];
	strcpy(Z0_formula, inputFormula.c_str());
	// std::cout << "  " << X0_formula << std::endl << "  " << Y0_formula << std::endl << "  " << Z0_formula << std::endl;
	// ────────────────────────────────────────────────────────────
	// Read the eight Γ̄ᵏᵢⱼ formulas
	// --- read eta scalar value ---
	double eta;
	std::cin >> eta;
	std::cout << "[DEBUG] Input eta = " << eta << std::endl;

    // --- first, read the connection‐energy Lamé parameters ---
    double lambdaG, muG;
    // std::cout << "Enter scalar λ_G for connection energy: " << std::endl;
    std::cin >> lambdaG;
    // std::cout << "Enter scalar μ_G for connection energy: " << std::endl;
    std::cin >> muG;

    // --- now read the eight Γ̄ᵏᵢⱼ formulas in lex order (k,i,j = 1..2) ---
    // std::cout << "=== DEBUG: entering gamma-read block ===" << std::endl;

    // std::cout << "Enter Γ̄^1_{11}(u,v): "<< std::endl;
    std::cin >> inputFormula;
    gammabar111_formula = new char[inputFormula.length()+1];
    strcpy(gammabar111_formula, inputFormula.c_str());

    // std::cout << "Enter Γ̄^1_{12}(u,v): "<< std::endl;
    std::cin >> inputFormula;
    gammabar112_formula = new char[inputFormula.length()+1];
    strcpy(gammabar112_formula, inputFormula.c_str());

    // std::cout << "Enter Γ̄^1_{21}(u,v): " << std::endl;
    std::cin >> inputFormula;
    gammabar121_formula = new char[inputFormula.length()+1];
    strcpy(gammabar121_formula, inputFormula.c_str());

    // std::cout << "Enter Γ̄^1_{22}(u,v): " << std::endl;
    std::cin >> inputFormula;
    gammabar122_formula = new char[inputFormula.length()+1];
    strcpy(gammabar122_formula, inputFormula.c_str());

    // std::cout << "Enter Γ̄^2_{11}(u,v): " << std::endl;
    std::cin >> inputFormula;
    gammabar211_formula = new char[inputFormula.length()+1];
    strcpy(gammabar211_formula, inputFormula.c_str());

    // std::cout << "Enter Γ̄^2_{12}(u,v): " << std::endl;
    std::cin >> inputFormula;
    gammabar212_formula = new char[inputFormula.length()+1];
    strcpy(gammabar212_formula, inputFormula.c_str());

    // std::cout << "Enter Γ̄^2_{21}(u,v): " << std::endl;
    std::cin >> inputFormula;
    gammabar221_formula = new char[inputFormula.length()+1];
    strcpy(gammabar221_formula, inputFormula.c_str());

    // std::cout << "Enter Γ̄^2_{22}(u,v): " << std::endl;
    std::cin >> inputFormula;
    gammabar222_formula = new char[inputFormula.length()+1];
    strcpy(gammabar222_formula, inputFormula.c_str());


	// std::cout << "Enter Initial thickness adjustment parameter " << std::endl;
	std::cin >> ThicknessAdjust;		// std::cout << "  " << ThicknessAdjust << std::endl;
	// std::cout << "Enter Initial metric adjustment parameter " << std::endl;
	std::cin >> MetricAdjust;			// std::cout << "  " << MetricAdjust << std::endl;

	// std::cout << "Enter Initial eta adjustment parameter " << std::endl;
	std::cin >> EtaAdjust;				// std::cout << "  " << EtaAdjust << std::endl;


	// std::cout << "Enter whether to restart (0=No, 1=restart) " << std::endl;
	std::cin >> restart;

	if (restart != 0)
	{
		// std::cout << "Enter the restart file name " << std::endl;
		std::cin >> restartFileName;
	}

	/* Setting the file names */
	// nodeOutputFileName  = verticesFileName + ".dat";
	// faceOutputFileName  = facesFileName + ".dat";
	// EFGLMNOutputFileName  = facesFileName + ".EFGLMN";
	// energyFileName      = verticesFileName + ".energy";
	// saveFileName        = verticesFileName + ".save";
	/* Setting the file names (put them alongside the input files) */
	{
	fs::path vPath(verticesFileName);
	fs::path fPath(facesFileName);
	fs::path dir = vPath.parent_path();          // e.g. "mesh/"
	std::string vStem = vPath.stem().string();   // e.g. "verts"
	std::string fStem = fPath.stem().string();   // e.g. "faces"

	nodeOutputFileName   = (dir / (vStem + ".dat")).string();
	energyFileName       = (dir / (vStem + ".energy")).string();
	saveFileName         = (dir / (vStem + ".save")).string();

	faceOutputFileName   = (dir / (fStem + ".dat")).string();
	EFGLMNOutputFileName = (dir / (fStem + ".EFGLMN")).string();
	}


	/* Create the file handles */
	TextFileHandle nodeOutputFileHandle(nodeOutputFileName,FileHandle::OPEN_WR);
	TextFileHandle faceOutputFileHandle(faceOutputFileName,FileHandle::OPEN_WR);
	TextFileHandle EFGLMNOutputFileHandle(EFGLMNOutputFileName,FileHandle::OPEN_WR);
	BinaryFileHandle saveFileHandle(saveFileName,FileHandle::OPEN_WR);
	TextFileHandle   energyFileHandle(energyFileName,FileHandle::OPEN_WR);
	// for(double u=0; u<=1; u+=0.5)
  	// for(double v=0; v<=1; v+=0.5)
    // std::cout << "thickness("<<u<<","<<v<<") = " << inputFunctionThickness(u,v) << "\n";

	/* Construct the NonEuclideanShell */
	NonEuclideanShell lattice(verticesFileName,facesFileName);

	/* If restart then do it (set the initial state to what's in file) */
	if (restart==1)
	{
		BinaryFileHandle restartFileHandle(restartFileName,FileHandle::OPEN_RD);
		lattice.restart(&restartFileHandle);
		// lattice.checkNodePositions();

	}
	else
	{
		lattice.defaultInitialization();
	}


	/* Set various parameters of the NonEuclideanShell */
	lattice.setVerbosity(1);
	// A trivial zero‐reference connection
    // // build Γ̄ from the user‐supplied formulas:
    // auto inputFunctionGammaBar = [](double u, double v) {
    //   RVar uvar("u",&u), vvar("v",&v);
    //   RVar* vars[2] = {&uvar,&vvar};
    //   ROperation op000(gammabar111_formula,2,vars),
    //              op001(gammabar112_formula,2,vars),
    //              op010(gammabar121_formula,2,vars),
    //              op011(gammabar122_formula,2,vars),
    //              op100(gammabar211_formula,2,vars),
    //              op101(gammabar212_formula,2,vars),
    //              op110(gammabar221_formula,2,vars),
    //              op111(gammabar222_formula,2,vars);

    //   TinyMatrix<TinyMatrix<double,2>,2> G;
    //   // row 0, col 0
    //   G(0,0)(0,0) = op000.Val(); G(0,0)(0,1) = op001.Val();
    //   G(0,1)(0,0) = op010.Val(); G(0,1)(0,1) = op011.Val();
    //   G(1,0)(0,0) = op100.Val(); G(1,0)(0,1) = op101.Val();
    //   G(1,1)(0,0) = op110.Val(); G(1,1)(0,1) = op111.Val();
    //   return G;
    // };
    // lambdaG, muG were read just above




	// --- Debug: Print per-face energy densities ---
	// const Vector<Face*>& faces = lattice.getFaces();
	// for (int i = 0; i < faces.length(); i++) {
	// 	std::cout << "Face " << i
	// 			<< " Stretch = " << faces(i)->stretchingEnergyContentDensity()
	// 			<< ", Bend = " << faces(i)->bendingEnergyContentDensity()
	// 			<< ", Conn = " << faces(i)->connectionEnergyContentDensity()
	// 			<< std::endl;
	// }




    lattice.setParameters(&inputFunctionThickness,
						  &inputFunctionLambda,
						  &inputFunctionMu,
						  &inputFunctionAbar,
						  &inputFunctionBbar,
						  &inputFunctionCbar,
                          &inputFunctionPos0,
                          &inputFunctionGammaBar,
						  eta,
                          lambdaG,
                          muG);
	lattice.setAdjust(ThicknessAdjust,
                  MetricAdjust,
                  EtaAdjust);



	// lattice.checkNodePositions(); // <-- This should ALSO print nothing

	// lattice.setAdjust(ThicknessAdjust, MetricAdjust);


	/* temporary: test the gradient */
	// lattice.testGradient();
	// exit(1);

	/* ******************************************************************* */
	/* OPTIMIZATION *** OPTIMIZATION *** OPTIMIZATION *** OPTIMIZATION *** */
	/* ******************************************************************* */

	/* Calculate the initial energy and output it */
	std::cout.precision(15);
	std::cout << "\tThe initial bending energy is " << lattice.bendingEnergy() * pow(ThicknessAdjust * MetricAdjust, 0) << std::endl;
	std::cout << "\tThe initial stretching energy is " << lattice.stretchingEnergy() * pow(ThicknessAdjust * MetricAdjust, 0) << std::endl;
	std::cout << "\tThe initial connection energy is " << lattice.connectionEnergy() * pow(ThicknessAdjust * MetricAdjust, 0) << std::endl;



	/* The size of the optimization problem */
	int size = lattice.SizeOfOptimizationProblem();

	/* define the minimizer: we work with BFGS */
	// const gsl_multimin_fdfminimizer_type* method = gsl_multimin_fdfminimizer_vector_bfgs;
	const gsl_multimin_fdfminimizer_type* method = gsl_multimin_fdfminimizer_conjugate_fr;
	// const gsl_multimin_fdfminimizer_type* method = gsl_multimin_fdfminimizer_vector_bfgs2;

	// const gsl_multimin_fdfminimizer_type* method = gsl_multimin_fdfminimizer_steepest_descent;
	gsl_multimin_fdfminimizer*      optimizer = gsl_multimin_fdfminimizer_alloc (method, size);

	/* Define the function to be minimized */
	gsl_multimin_function_fdf my_func;
	my_func.n      = size;
	my_func.f      = &my_f;
	my_func.df     = &my_df;
	my_func.fdf    = &my_fdf;
	my_func.params = static_cast<void *>(&lattice);

	/* construct the initial state */
	gsl_vector *IC;
	IC = gsl_vector_alloc(size);
	lattice.getPositionVector(IC);

	/* set paramters of optimizer */
	// gsl_multimin_fdfminimizer_set(optimizer, &my_func, IC, 0.01, 1e-6);
	
	/* set paramters of optimizer */
	gsl_multimin_fdfminimizer_set(optimizer, &my_func, IC, 0.01, 1e-6);	
	


	/* the optimization loop */
	int status;
	int print_counter=0;
	for (int iter=0; iter<NumberOfLoops; iter++)
	{
		/* set the adjustment parameters */
		double thicknessAdjustSign = 1.0;
		
		if (ThicknessAdjust < 1.0)	{thicknessAdjustSign = -1.0;}
		// double adjustParamThickness	= 1.0 + (ThicknessAdjust - 1.0) * pow((1.0 - 5.0 * iter / NumberOfLoops),3);
		double adjustParamThickness	= ThicknessAdjust + (1.0 - ThicknessAdjust) * pow(5.0 * iter / NumberOfLoops,3);
		/*double adjustParamThickness	= ThicknessAdjust + (1.0 - ThicknessAdjust) * 0.5*(1 + sin(3.141592*(5.0 * iter / NumberOfLoops) - 0.727)/abs(sin(3.141592*(5.0 * iter / NumberOfLoops) - 0.727))); */
		// double adjustParamThickness	= pow(ThicknessAdjust + (1.0 - ThicknessAdjust) * iter / NumberOfLoops,-2);//adjust wgal then calibrate
		if (adjustParamThickness*thicknessAdjustSign < thicknessAdjustSign)	{adjustParamThickness = 1.0;}
		double adjustParamMetric = 1.0 + (MetricAdjust - 1.0) * (1.0 - 1.0 * iter / NumberOfLoops);
		if (adjustParamMetric < 1)		{adjustParamMetric = 1.0;}
		double adjustParamEta = 1.0 + (EtaAdjust - 1.0) * pow((1.0 - 5.0 * iter / NumberOfLoops), 3);
		// if (iter == 0) {
		// 	std::cout << "[DEBUG] Initial eta = " << eta << std::endl;
		// 	std::cout << "[DEBUG] Initial adjustParamEta = " << adjustParamEta << std::endl;
		// }
		lattice.setAdjust(adjustParamThickness, adjustParamMetric, adjustParamEta);


		
		// double adjustParamThickness = ThicknessAdjust + (1.0 - ThicknessAdjust) * pow(5.0 * iter / NumberOfLoops,3);
		// double adjustParamMetric = 1.0 + (MetricAdjust - 1.0) * (1.0 - 1.0 * iter / NumberOfLoops);
		// if (adjustParamMetric < 1.0) adjustParamMetric = 1.0;
		// lattice.setAdjust(adjustParamThickness, adjustParamMetric);
		
		// if (iter % HowOftenToPrint == 0) {
    	// 	std::cout << "iter " << iter
        //       << ": adjustThickness=" << adjustParamThickness
        //       << ", adjustMetric=" << adjustParamMetric
        //       << ", stretching=" << lattice.stretchingEnergy()
        //       << ", bending=" << lattice.bendingEnergy()
        //       << ", connection=" << lattice.connectionEnergy()
        //       << ", total=" << (lattice.stretchingEnergy() + lattice.bendingEnergy() + lattice.connectionEnergy())
        //       << std::endl;
		// }



		// Temporary fixed adjustment parameters
		// double adjustParamThickness = 1.0;
		// double adjustParamMetric = 1.0;
		/* perform an iteration */
		status = gsl_multimin_fdfminimizer_iterate(optimizer);
	    // if (iter % 100 == 0 ) {
        // std::cout << "Iter " << iter
        //           << " E = " << lattice.energy()
        //           << " stretch = " << lattice.stretchingEnergy()
        //           << " bend = " << lattice.bendingEnergy()
        //           << " conn = " << lattice.connectionEnergy()
        //           << std::endl;
        // std::cout << "Gradient norm: " << gsl_blas_dnrm2(optimizer->gradient) << std::endl;
    	// }

		// if (status == GSL_ENOPROG) {
		// 	std::cout << "Line search stagnated at iter " << iter << " — exiting\n";
		// 	break;
		// } else if (status) {
		// 	std::cerr << "Optimization stopped: " << gsl_strerror(status) << std::endl;
		// 	break;
		// }


		// if (!status) {
		// // update the shell to the new x
		// double* currentX = gsl_vector_ptr(optimizer->x, 0);
		// lattice.setPositionVector(currentX);
		
		// }
		

		if (status)
		{
			if (((ThicknessAdjust == 1.0) && (MetricAdjust == 1.0)) || ((5.0 * iter / NumberOfLoops) > 1.0)) /*Allow exit only after 1/5 of maximum number of iterations (if adjusted)*/
			{
				// std::cout << "Exit minimizer with status " << status << std::endl;
				std::cout << "Exit minimizer with status " << status
				<< " (" << gsl_strerror(status) << ")\n";

				break;
			}
			
		}
		    // — pull the optimizer’s current x-vector back into the lattice —
		{
			double* currentX = gsl_vector_ptr(optimizer->x, 0);
			lattice.setPositionVector(currentX);
		}

		


		/* check the size of the gradient */
		status = gsl_multimin_test_gradient(optimizer->gradient, 1e-6);
		if (status == GSL_SUCCESS)
		{
			if (((ThicknessAdjust == 1.0) && (MetricAdjust == 1.0)) ||((5.0 * iter / NumberOfLoops) > 1.0)) /*Allow exit only after 1/5 of maximum number of iterations (if adjusted)*/
			{
				std::cout << "Minimum found" << std::endl;
				break;
			}
		}
		

		if (iter%HowOftenToPrint==0)
		{
		/* print the current energy and thickness */
			char *space = (char*)(" ");
			double Es = lattice.stretchingEnergy() * pow(adjustParamThickness * adjustParamMetric, 0);
			double Eb = lattice.bendingEnergy() * pow(adjustParamThickness * adjustParamMetric, 0);
			double Eg = lattice.connectionEnergy()  * pow(adjustParamThickness * adjustParamMetric, 0);
			// double Es = lattice.stretchingEnergy();
			// double Eb = lattice.bendingEnergy();
			// double Eg = lattice.connectionEnergy();
    		double E  = Es + Eb + Eg;
			energyFileHandle.write(iter);
			energyFileHandle.write(space);
			energyFileHandle.write(E);
			energyFileHandle.write(space);
			energyFileHandle.write(Es);
			energyFileHandle.write(space);
			energyFileHandle.write(Eb);
			energyFileHandle.write(space);
			energyFileHandle.write(Eg);
    		energyFileHandle.newline();
			energyFileHandle.write(adjustParamMetric);
			energyFileHandle.newline();

		/* Dump the surface to the openGL file */
			double Efinal = lattice.stretchingEnergy() + lattice.bendingEnergy() + lattice.connectionEnergy();
			int PercentDone = 100*iter/NumberOfLoops;
			std::cout.precision(15);
			std::cout << PercentDone << "% Done: " << "The current energy is " << Efinal << std::endl;
			std::cout << " stretch = " << lattice.stretchingEnergy()
				<< " bend = " << lattice.bendingEnergy()
				<< " conn = " << lattice.connectionEnergy()
				<< std::endl;
			lattice.DumpStateTextFormat(&nodeOutputFileHandle, &faceOutputFileHandle, ++print_counter);
			lattice.DumpFormsTextFormat(&EFGLMNOutputFileHandle);
		}
	// 		if (iter % 100 == 0) {
    // 	std::cout << "x[0] = " << gsl_vector_get(optimizer->x, 0) << std::endl;
	// }

	}



	/* Set the state of the system with final state */
	double* ptr = gsl_vector_ptr(optimizer->x, 0);
	lattice.setPositionVector(ptr);

	/* free the memory */
	gsl_multimin_fdfminimizer_free(optimizer);
	gsl_vector_free(IC);

		/* Calculate the final energy and output it */
	std::cout.precision(15);
	std::cout << "The final energy is " << lattice.energy() << std::endl;

	char *space = (char*)(" ");
	double Es = lattice.stretchingEnergy();
	double Eb = lattice.bendingEnergy();
	double Eg = lattice.connectionEnergy();
	double E = Es + Eb + Eg;
	energyFileHandle.write(0);
	energyFileHandle.write(space);
	energyFileHandle.write(E);
	energyFileHandle.write(space);
	energyFileHandle.write(Es);
	energyFileHandle.write(space);
	energyFileHandle.write(Eb);
	energyFileHandle.write(space);
	energyFileHandle.write(Eg);
	energyFileHandle.write(space);
	energyFileHandle.write(1.0);
	energyFileHandle.newline();

    // Save the results
    // lattice.DumpStateBinaryFormat(&saveFileHandle);
    lattice.DumpFormsTextFormat(&EFGLMNOutputFileHandle);
    lattice.DumpStateTextFormat(&nodeOutputFileHandle, &faceOutputFileHandle, ++print_counter);

    // --- DEBUG: verify that the final dump actually ran ---
    // std::cout << "=== DEBUG: final dump complete ===\n";

    /* Close the file handles */
    nodeOutputFileHandle.close();
    faceOutputFileHandle.close();
    saveFileHandle.close();
    energyFileHandle.close();
    EFGLMNOutputFileHandle.close();  

	return 0;
}


/* ============================================================================== */
/* FUNCTIONS FOR OPTIMIZER                                                        */
/* ============================================================================== */
double my_f(const gsl_vector *a_x, void *a_lattice)
{
	NonEuclideanShell *lattice_ptr = static_cast<NonEuclideanShell*>(a_lattice);
	return lattice_ptr->getEnergy(a_x);
}

void my_df(const gsl_vector *a_x, void *a_lattice, gsl_vector *a_df)
{
	NonEuclideanShell *lattice_ptr = static_cast<NonEuclideanShell*>(a_lattice);
	lattice_ptr->getEnergyGradient(a_x, a_df);
}

void my_fdf(const gsl_vector *a_x, void *a_lattice, double *a_f, gsl_vector *a_df)
{
	NonEuclideanShell *lattice_ptr = static_cast<NonEuclideanShell*>(a_lattice);
	lattice_ptr->getEnergyAndEnergyGradient(a_x, a_f, a_df);
}

/* ============================================================================== */
/* INPUT FUNCTIONS (lattice parameters)                                           */
/* ============================================================================== */
TinyMatrix<double,2>  inputFunctionAbar(double u, double v)
{

	RVar uvar ("u", &u);
	RVar vvar ("v", &v);
	RVar* vararray[2];
	vararray[0]=&uvar;
	vararray[1]=&vvar;
	ROperation op11(abar11_formula, 2, vararray);
	ROperation op12(abar12_formula, 2, vararray);
	ROperation op21(abar21_formula, 2, vararray);
	ROperation op22(abar22_formula, 2, vararray);

	TinyMatrix<double,2> ret;
	ret(0,0) = op11.Val();
	ret(0,1) = op12.Val();
	ret(1,0) = op21.Val();
	ret(1,1) = op22.Val();
	return ret;
}

TinyMatrix<double,2> inputFunctionBbar(double u, double v)
{
	RVar uvar ("u", &u);
	RVar vvar ("v", &v);
	RVar* vararray[2];
	vararray[0]=&uvar;
	vararray[1]=&vvar;
	ROperation op11(bbar11_formula, 2, vararray);
	ROperation op12(bbar12_formula, 2, vararray);
	ROperation op21(bbar21_formula, 2, vararray);
	ROperation op22(bbar22_formula, 2, vararray);

	TinyMatrix<double,2> ret;
	ret(0,0) = op11.Val();
	ret(0,1) = op12.Val();
	ret(1,0) = op21.Val();
	ret(1,1) = op22.Val();
	return ret;
}

TinyMatrix<double,2> inputFunctionCbar(double u, double v)
{
    RVar uvar ("u", &u);
    RVar vvar ("v", &v);

    RVar* vararray[2];
    vararray[0] = &uvar;
    vararray[1] = &vvar;

    ROperation op11(cbar11_formula, 2, vararray);
    ROperation op12(cbar12_formula, 2, vararray);
    ROperation op21(cbar21_formula, 2, vararray);
    ROperation op22(cbar22_formula, 2, vararray);

    TinyMatrix<double,2> ret;

    ret(0,0) = op11.Val();
    ret(0,1) = op12.Val();
    ret(1,0) = op21.Val();
    ret(1,1) = op22.Val();

    return ret;
}

double inputFunctionThickness(double u, double v)
{
	RVar uvar ("u", &u);
	RVar vvar ("v", &v);
	RVar* vararray[2];
	vararray[0]=&uvar;
	vararray[1]=&vvar;
	ROperation op(thickness_formula, 2, vararray);

	return op.Val();
}

double inputFunctionLambda(double u, double v)
{
	RVar uvar ("u", &u);
	RVar vvar ("v", &v);
	RVar* vararray[2];
	vararray[0]=&uvar;
	vararray[1]=&vvar;
	ROperation opY(Y_formula, 2, vararray);
	ROperation opnu(nu_formula, 2, vararray);
	double Y  = opY.Val();
	double nu = opnu.Val();

	return Y * nu / (1 - nu * nu) / 8;
}

double inputFunctionMu (double u, double v)
{
	RVar uvar ("u", &u);
	RVar vvar ("v", &v);
	RVar* vararray[2];
	vararray[0]=&uvar;
	vararray[1]=&vvar;
	ROperation opY(Y_formula, 2, vararray);
	ROperation opnu(nu_formula, 2, vararray);
	double Y  = opY.Val();
	double nu = opnu.Val();

	return Y / (nu + 1) / 8;
}

TinyVector<double,3>  inputFunctionPos0(double u, double v)
{

	RVar uvar ("u", &u);
	RVar vvar ("v", &v);
	RVar* vararray[2];
	vararray[0]=&uvar;
	vararray[1]=&vvar;
	ROperation op1(X0_formula, 2, vararray);
	ROperation op2(Y0_formula, 2, vararray);
	ROperation op3(Z0_formula, 2, vararray);

	TinyVector<double,3> ret;
	ret(0) = op1.Val();
	ret(1) = op2.Val();
	ret(2) = op3.Val();
	return ret;
}