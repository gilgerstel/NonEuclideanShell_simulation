/*
 *  NonEuclideanShell.H
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/8/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

#ifndef _NONEUCLIDEANSHELL_H_
#define _NONEUCLIDEANSHELL_H_

#include "Main.H"
#include "TinyVector.H"
#include "TinyMatrix.H"
#include "Matrix.H"
#include "Errors.H"
#include "MatlabFileHandle.H"
#include "BinaryFileHandle.H"
#include "TextFileHandle.H"
#include "gsl/gsl_vector.h"
#include <string>
#include <vector> 
#include <set>


/*
 class Node

 A node consists of:
 a TinyVector<double,2> holding the coordinates (u,v).
 a TinyVector<double,3> holding the location in R3.
 a TinyVector<double,3> holding the force exerted on the node.
 The Node knows nothing about the lattice.

*/

class Node
	{
	public:

		/* Default constructor */
		Node() :
			m_coordinates(),
			m_position(),
			m_force(),
			m_fixed(),
			m_offset()
		{}

		/* Copy constructor */
		Node(const Node &a_rhs) :
			m_coordinates(a_rhs.m_coordinates),
			m_position   (a_rhs.m_position),
			m_force      (a_rhs.m_force),
			m_fixed		 (a_rhs.m_fixed),
			m_offset	 (a_rhs.m_offset)
		{}

		/* Assignment */
		void operator=(const Node &a_rhs)
		{
			m_coordinates = a_rhs.m_coordinates;
			m_position    = a_rhs.m_position;
			m_force       = a_rhs.m_force;
			m_fixed		  = a_rhs.m_fixed;
			m_offset	  = a_rhs.m_offset;
		}

		/* References to data (const and non-const) */
		TinyVector<double,3>& position()          {return m_position;}
		TinyVector<double,3>& force()             {return m_force;}
		TinyVector<double,2>& coordinates()       {return m_coordinates;}
		TinyVector<double,3>& offset()			  {return m_offset;}
		double&               position(int i)     {return m_position(i);}
		double&               force(int i)        {return m_force(i);}
		double&               coordinates(int i)  {return m_coordinates(i);}
		int&				  fixed()			  {return m_fixed;}
		double&               offset(int i)        {return m_offset(i);}


		const TinyVector<double,3>& position()          const {return m_position;}
		const TinyVector<double,3>& force()             const {return m_force;}
		const TinyVector<double,2>& coordinates()       const {return m_coordinates;}
		const TinyVector<double,3>& offset()             const {return m_offset;}
		double                      position(int i)     const {return m_position(i);}
		double                      force(int i)        const {return m_force(i);}
		double                      coordinates(int i)  const {return m_coordinates(i);}
		int	                        fixed()				const {return m_fixed;}
		double                      offset(int i)        const {return m_offset(i);}

		/* A pointed to the position vector */
		double* getPointer()        {return m_position.getPointer();}

	private:

		TinyVector<double,2>   m_coordinates;
		TinyVector<double,3>   m_position;
		TinyVector<double,3>   m_force;
		int					   m_fixed;
		TinyVector<double,3>   m_offset;
	};


/*
 class Face

 A face consists of:
 a TinyVector<Node*,6> holding pointes to 6 nodes. The first 3 are the vertices of
	the triangle (in CCW order). The next 3 are the next-nearest-neighbors.
 a TinyVector<Face*,3> holding pointers to the 3 neighboring Faces.
 a TinyVector<double,2> holding the coordinates of the center of the Face.
 a TinyVector<double,2> holding the coordinates difference between Nodes 1 and 2
 a TinyVector<double,2> holding the coordinates difference between Nodes 2 and 3
 a TinyVector<double,2> holding the coordinates difference between Nodes 3 and 1
 3 TinyVector<double,2> holding the coordinates difference between the Face center and its neighbors
 the area of the Face (cooridnates)
 The parameters Lambda, Mu, and the Thickness.
 2 TinyMatrix<double,2> for abar and bbar

 A face cannot be constructed empty, nor be assigned.

*/

class Face
	{
	public:

		/* Forbid empty constructor */
		Face();

		/* Constructor with pointers to Nodes */
		/* All the other members are assigned separately */
		Face(const TinyVector<Node*,6> &a_nodes);

		/* Forbid copy constructor */
		Face(const Face &a_rhs);

		/* Forbid assignemnt */
		void operator=(const Face &);

		/* Assign neighboring faces */
		void assignNeighbors(const TinyVector<Face*,3> &a_faces);
		/* Initialize */
		/* Gets the parameters of the Face */
		/* Calculates there the auxiliary matrices */
        void initialize(double a_thickness,
                        double a_lambda,
                        double a_mu,
                        const TinyMatrix<double,2> a_abar,
                        const TinyMatrix<double,2> a_bbar,
						const TinyMatrix<double,2> a_cbar,
						double a_cPrefactor,
						double a_K1111,
						double a_K1112,
						double a_K1122,
						double a_K1212,
						double a_K1222,
						double a_K2222,
                        const TinyMatrix< TinyMatrix<double,2>,2 >& a_gammabar,
						double a_eta,
                        double a_lambdaG,
                        double a_muG);




		/* Modify adjustment parameter */
		// void setAdjust(double a_adjust1, double a_adjust2) {m_adjust1 = a_adjust1; m_adjust2 = a_adjust2;}
		void setAdjust(double a_adjust1, double a_adjust2, double a_adjust3 = 1.0) {
			m_adjust1 = a_adjust1;
			m_adjust2 = a_adjust2;
			m_adjust3 = a_adjust3;
		}

		/* return the position */
		TinyVector<double,3> position() const;

		/* return the coordinates */
		TinyVector<double,2>  coordinates() const {return m_coordinates;}
		TinyVector<double,2>& coordinates()       {return m_coordinates;}
		double  coordinates(int i) const {return m_coordinates(i);}
		double& coordinates(int i)       {return m_coordinates(i);}

		/* Calculate the unit normal */
		TinyVector<double,3> calculateUnitNormal() const;

       /* Calculate energy densities */
		double stretchingEnergyContentDensity() const;
		double bendingEnergyContentDensity() const;
		double connectionEnergyContentDensity() const;
		double thirdFundamentalEnergyContentDensity() const;

		double stretchingEnergy() const;
		// {return pow(m_adjust2,-6) * stretchingEnergyContentDensity() * m_area * m_thickness * m_adjust1;}
		double bendingEnergy() const;
		//  const {return pow(m_adjust2,-6) * bendingEnergyContentDensity() * m_area * pow(m_thickness * m_adjust1,3);}
		double connectionEnergy() const;
		double thirdFundamentalEnergy() const;
		
		//  {return pow(m_adjust2,-6) * connectionEnergyContentDensity() * m_area * pow(m_thickness * m_adjust1,3);}
		double energy() const {
			return stretchingEnergy()
				+ bendingEnergy()
				+ connectionEnergy()
				+ thirdFundamentalEnergy();
		}

		/* Calculate fundamental forms */
		TinyVector<double,3> EFG() const;
		TinyVector<double,3> LMN() const;
		TinyVector<double,3> ABC() const;

        // /* γ-energy helpers */
        // TinyMatrix<double,2> computeMetric() const;
        // std::pair<TinyMatrix<double,2>,TinyMatrix<double,2>> computeMetricDerivatives() const;
        // TinyMatrix<TinyMatrix<double,2>,2> computeConnection() const;

		    /* γ-energy helpers */
		TinyMatrix<double,2> computeMetric() const;
		std::pair<TinyMatrix<double,2>,TinyMatrix<double,2>> computeMetricDerivatives() const;
		TinyMatrix<TinyMatrix<double,2>,2> computeConnection() const;

	

		    // in class Face { … }

		// access one of the three neighbor pointers
		Face* neighbor(int e) const { return m_faces(e); }

		// read-only access to your stored reference γ-bar and weights
		const TinyMatrix< TinyMatrix<double,2>,2 >& gammabar() const { return m_gammabar; }
		const TinyMatrix<double,2>& cbar() const { return m_cbar; }
		double cPrefactor() const { return m_cPrefactor; }
		double K1111() const { return m_K1111; }
		double K1112() const { return m_K1112; }
		double K1122() const { return m_K1122; }
		double K1212() const { return m_K1212; }
		double K1222() const { return m_K1222; }
		double K2222() const { return m_K2222; }
		double lambdaG() const { return m_lambdaG; }
		double muG() const     { return m_muG; }

        /* Calculate forces (energy gradient) */
        void setForce();

        /* Output the area */
        double area() const { return m_area; }
		double getLambdaG() const { return m_lambdaG; }
		double getMuG() const { return m_muG; }
		Node* node(int i) const { return m_nodes(i); }    // 0,1,2 are the triangle’s vertices



	private:
		
		TinyVector<Node*,6> m_nodes;
		TinyVector<Face*,3> m_faces;

		TinyVector<double,2> m_coordinates;
		TinyVector<double,2> m_edge12;
		TinyVector<double,2> m_edge23;
		TinyVector<double,2> m_edge31;
		TinyVector<double,2> m_nodeconnector1;
		TinyVector<double,2> m_nodeconnector2;
		TinyVector<double,2> m_nodeconnector3;
		TinyVector<double,2> m_nodeconnector4;
		TinyVector<double,2> m_nodeconnector5;
		TinyVector<double,2> m_nodeconnector6;
		double               m_area;
		double			     m_lambda;
		double               m_mu;
		double               m_thickness;
		double               m_adjust1;
		double               m_adjust2;
		double               m_adjust3;
		// Lamé parameters for the γ-energy term
		double m_eta;

        double m_lambdaG;
        double m_muG;

		TinyMatrix<double,3> m_Amatrix;
		TinyVector<int,3> 	 m_Apivot;
		TinyMatrix<double,6> m_Bmatrix;
		TinyVector<int,6>	 m_Bpivot;
		TinyMatrix<double,2> m_abar;
		TinyMatrix<double,2> m_bbar;
		TinyMatrix<double,2> m_cbar;
		// Effective third-fundamental-form energy parameters
		double m_cPrefactor;
		double m_K1111;
		double m_K1112;
		double m_K1122;
		double m_K1212;
		double m_K1222;
		double m_K2222;
		TinyMatrix<double,2> m_invabar;
		TinyMatrix< TinyMatrix<double,2>, 2 > m_gammabar;



	};


/*
 class NonEuclideanShell

 A NonEuclideanShell consists of
 a Vector<Node*>
 a Vector<Face*>
 a BoundaryConditionsType (currently inactive)

 */


class NonEuclideanShell
	{
	public:
		void checkNodePositions() const;
	    double meanEdgeLength() const;
		/// Return all mesh‐adjacent vertex indices for vertex v
		std::vector<int> getVertexNeighbors(int v) const;


		const Vector<Face*>& getFaces() const;

        void restart(BinaryFileHandle* a_fh);
        void getEnergyGradientFull(const gsl_vector* a_state, gsl_vector* a_gradient);
        void getEnergyAndEnergyGradientFull(const gsl_vector* a_state, double* a_energy, gsl_vector* a_gradient);

		/* Constructor */
		NonEuclideanShell(const std::string &nodesFileName,
						  const std::string &facesFileName);

		/* Destructor (cleanup) */
		~NonEuclideanShell();

		/* Set the verbosity */
		void setVerbosity(int a_verbosity) {m_verbosity = a_verbosity;}

		/* Set the parameters */
//       void setParameters(double (*)(double, double),
//						   double (*)(double, double),
//						   double (*)(double, double),
//						   TinyMatrix<double,2> (*)(double, double),
//						   TinyMatrix<double,2> (*)(double, double),
//						   TinyVector<double,3> (*)(double, double),
//						   TinyMatrix<TinyMatrix<double,2>,2> (*gammabarFunc)(double u, double v),
//                           double lambdaG,
//                           double muG);
        /** Set the reference forms and γ-term weights */
        // void setParameters(
        //     double (*f_thickness)(double,double),
        //     double (*f_lambda)(double,double),
        //     double (*f_mu)(double,double),
        //     TinyMatrix<double,2> (*f_abar)(double,double),
        //     TinyMatrix<double,2> (*f_bbar)(double,double),
        //     TinyVector<double,3> (*f_Pos0)(double,double),
        //     TinyMatrix<TinyMatrix<double,2>,2> (*f_gammabar)(double,double),
        //     double lambdaG,
        //     double muG
        // );
		void setParameters(
		double (*f_thickness)(double,double),
		double (*f_lambda)(double,double),
		double (*f_mu)(double,double),
		TinyMatrix<double,2> (*f_abar)(double,double),
		TinyMatrix<double,2> (*f_bbar)(double,double),
		TinyMatrix<double,2> (*f_cbar)(double,double),

		double (*f_cPrefactor)(double,double),
		double (*f_K1111)(double,double),
		double (*f_K1112)(double,double),
		double (*f_K1122)(double,double),
		double (*f_K1212)(double,double),
		double (*f_K1222)(double,double),
		double (*f_K2222)(double,double),

		TinyVector<double,3> (*f_Pos0)(double,double),
		TinyMatrix<TinyMatrix<double,2>,2> (*f_gammabar)(double,double),
		double eta,
		double lambdaG,
		double muG
	);



        /* set the adjustment parameters */
		// void setAdjust(double a_adjust1, double a_adjust2);
	void setAdjust(double a_adjust1, double a_adjust2, double a_adjust3 = 1.0);


		/* Defaults initialization of positions */
		void defaultInitialization();


    void initializeForce();
    void setForce();

    double stretchingEnergy() const;
    double bendingEnergy()    const;
    double connectionEnergy() const;
	double thirdFundamentalEnergy() const;
    double energy() const;

    int SizeOfOptimizationProblem() const;
    void setPositionVector(const double*);
    void setPositionVector(const gsl_vector*);
    void getPositionVector(gsl_vector*);

    double getEnergy(const gsl_vector*);
    void   getEnergyGradient(const gsl_vector*, gsl_vector*);
    void   getEnergyAndEnergyGradient(const gsl_vector*, double*, gsl_vector*);
    void   testGradient();

    void DumpStateBinaryFormat(BinaryFileHandle*);
    void DumpStateTextFormat(TextFileHandle*, TextFileHandle*, int);
    void DumpFormsTextFormat(TextFileHandle*);

private:
    Vector<Node*>                    m_nodes;
    Vector<Face*>                    m_faces;
    int                              m_verbosity;
	double 							 m_eta;
    bool                             m_includeConn = false;
};

#endif // _NONEUCLIDEANSHELL_H_
