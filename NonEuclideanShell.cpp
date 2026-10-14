 /*  NonEuclideanShell.cpp
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/9/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

#include "NonEuclideanShell.H"
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <map>
// #include "LapackWrapper.H"   // obsolete: luDecompose/luSolve are defined in TinyMatrix.H
#include <algorithm>   // for std::clamp
#include <limits>




/* ==============================================================================  */
/* Face Face Face Face Face Face Face Face Face Face Face Face Face Face Face Face */
/* ==============================================================================  */
/* Constructor with pointers to Nodes */
Face::Face(const TinyVector<Node*,6> &a_nodes) :
m_nodes(a_nodes),
m_faces(),
m_coordinates(),
m_edge12(),
m_edge23(),
m_edge31(),
m_nodeconnector1(),
m_nodeconnector2(),
m_nodeconnector3(),
m_nodeconnector4(),
m_nodeconnector5(),
m_nodeconnector6(),
m_Amatrix(),
m_Apivot(),
m_Bmatrix(),
m_Bpivot(),
m_abar(),
m_bbar(),
m_cbar(),
m_invabar(),
m_adjust1(1.0),
m_adjust2(1.0),
m_adjust3(1.0),
// **new members** — defaulted to zero/identity as appropriate
m_gammabar(),      // default‐constructed 2×2 of 2×2 (all zeros)
m_lambdaG(0.0),
m_muG(0.0)
{
	/* Check that the first 3 Node* are not NULL */
	// assert(m_nodes(0)!=NULL && m_nodes(1)!=NULL && m_nodes(2)!=NULL);

	/* Set the Face coordinates to the average of the Nodes */
	m_coordinates = (m_nodes(0)->coordinates() + m_nodes(1)->coordinates() + m_nodes(2)->coordinates());
	m_coordinates.scale(1.0/3.0);

	/* Set the delta-Nodes */
	m_edge12      = m_nodes(1)->coordinates() - m_nodes(0)->coordinates();
	m_edge31      = m_nodes(0)->coordinates() - m_nodes(2)->coordinates();
	m_edge23      = m_nodes(2)->coordinates() - m_nodes(1)->coordinates();
}

/* ============================================================================== */
/* Assign pointers to neighboring faces */
void Face::assignNeighbors(const TinyVector<Face*,3> &a_faces)
{
	m_faces = a_faces;
}

/* ============================================================================== */
/* Initialization */
void Face::initialize(double a_thickness,
                      double a_lambda,
                      double a_mu,
                      const TinyMatrix<double,2> a_abar,
                      const TinyMatrix<double,2> a_bbar,
                      const TinyMatrix<double,2> a_cbar,
                      const TinyMatrix< TinyMatrix<double,2>,2 >& a_gammabar,
                      double a_eta,
                      double a_lambdaG,
                      double a_muG)
{
	/* Assign the input to member data */
	m_thickness = a_thickness;
	m_lambda    = a_lambda;
	m_mu        = a_mu;
	m_abar      = a_abar;
	m_bbar      = a_bbar;
	m_cbar      = a_cbar;
    m_gammabar  = a_gammabar;
	m_eta       = a_eta;

    m_lambdaG   = a_lambdaG;
    m_muG       = a_muG;

	/* Calculate the inverse of abar */
	m_invabar  = m_abar.inverse();

	/* Calculate the (reference) area */
	m_area      = 0.5 * sqrt(m_abar.det()) * fabs(m_edge12(0)*m_edge31(1) - m_edge12(1)*m_edge31(0));

	/* Connectors to neighboring nodes */
	m_nodeconnector1 = m_nodes(0)->coordinates() - m_coordinates;
	m_nodeconnector2 = m_nodes(1)->coordinates() - m_coordinates;
	m_nodeconnector3 = m_nodes(2)->coordinates() - m_coordinates;
	m_nodeconnector4(0) = 0.023;
	m_nodeconnector4(1) = 0.321;
	m_nodeconnector5(0) = 0.432;
	m_nodeconnector5(1) = 0.923;
	m_nodeconnector6(0) = 0.754;
	m_nodeconnector6(1) = 0.147;
	if (m_nodes(3) != NULL) m_nodeconnector4 = m_nodes(3)->coordinates() - m_coordinates;
	if (m_nodes(4) != NULL) m_nodeconnector5 = m_nodes(4)->coordinates() - m_coordinates;
	if (m_nodes(5) != NULL) m_nodeconnector6 = m_nodes(5)->coordinates() - m_coordinates;

	/* Construct the A matrix */
	m_Amatrix(0,0) =   m_edge23(0)*m_edge23(0);
	m_Amatrix(0,1) = 2*m_edge23(0)*m_edge23(1);
	m_Amatrix(0,2) =   m_edge23(1)*m_edge23(1);
	m_Amatrix(1,0) =   m_edge31(0)*m_edge31(0);
	m_Amatrix(1,1) = 2*m_edge31(0)*m_edge31(1);
	m_Amatrix(1,2) =   m_edge31(1)*m_edge31(1);
	m_Amatrix(2,0) =   m_edge12(0)*m_edge12(0);
	m_Amatrix(2,1) = 2*m_edge12(0)*m_edge12(1);
	m_Amatrix(2,2) =   m_edge12(1)*m_edge12(1);
	luDecompose(m_Amatrix, m_Apivot);
	double detA = m_Amatrix.det();
	if (fabs(detA) < 1e-12) {
		std::cerr << "Warning: m_Amatrix nearly singular!" << std::endl;
}


	/* Construct the B matrix */
	m_Bmatrix(0,0) = 1.0;
	m_Bmatrix(0,1) = m_nodeconnector1(0);
	m_Bmatrix(0,2) = m_nodeconnector1(1);
	m_Bmatrix(0,3) = 0.5*m_nodeconnector1(0)*m_nodeconnector1(0);
	m_Bmatrix(0,4) = m_nodeconnector1(0)*m_nodeconnector1(1);
	m_Bmatrix(0,5) = 0.5*m_nodeconnector1(1)*m_nodeconnector1(1);
	m_Bmatrix(1,0) = 1.0;
	m_Bmatrix(1,1) = m_nodeconnector2(0);
	m_Bmatrix(1,2) = m_nodeconnector2(1);
	m_Bmatrix(1,3) = 0.5*m_nodeconnector2(0)*m_nodeconnector2(0);
	m_Bmatrix(1,4) = m_nodeconnector2(0)*m_nodeconnector2(1);
	m_Bmatrix(1,5) = 0.5*m_nodeconnector2(1)*m_nodeconnector2(1);
	m_Bmatrix(2,0) = 1.0;
	m_Bmatrix(2,1) = m_nodeconnector3(0);
	m_Bmatrix(2,2) = m_nodeconnector3(1);
	m_Bmatrix(2,3) = 0.5*m_nodeconnector3(0)*m_nodeconnector3(0);
	m_Bmatrix(2,4) = m_nodeconnector3(0)*m_nodeconnector3(1);
	m_Bmatrix(2,5) = 0.5*m_nodeconnector3(1)*m_nodeconnector3(1);
	m_Bmatrix(3,0) = 1.0;
	m_Bmatrix(3,1) = m_nodeconnector4(0);
	m_Bmatrix(3,2) = m_nodeconnector4(1);
	m_Bmatrix(3,3) = 0.5*m_nodeconnector4(0)*m_nodeconnector4(0);
	m_Bmatrix(3,4) = m_nodeconnector4(0)*m_nodeconnector4(1);
	m_Bmatrix(3,5) = 0.5*m_nodeconnector4(1)*m_nodeconnector4(1);
	m_Bmatrix(4,0) = 1.0;
	m_Bmatrix(4,1) = m_nodeconnector5(0);
	m_Bmatrix(4,2) = m_nodeconnector5(1);
	m_Bmatrix(4,3) = 0.5*m_nodeconnector5(0)*m_nodeconnector5(0);
	m_Bmatrix(4,4) = m_nodeconnector5(0)*m_nodeconnector5(1);
	m_Bmatrix(4,5) = 0.5*m_nodeconnector5(1)*m_nodeconnector5(1);
	m_Bmatrix(5,0) = 1.0;
	m_Bmatrix(5,1) = m_nodeconnector6(0);
	m_Bmatrix(5,2) = m_nodeconnector6(1);
	m_Bmatrix(5,3) = 0.5*m_nodeconnector6(0)*m_nodeconnector6(0);
	m_Bmatrix(5,4) = m_nodeconnector6(0)*m_nodeconnector6(1);
	m_Bmatrix(5,5) = 0.5*m_nodeconnector6(1)*m_nodeconnector6(1);
	luDecompose(m_Bmatrix, m_Bpivot);
}

/* ============================================================================== */
/* Calculate the position of the Face */
TinyVector<double,3> Face::position() const
{
	TinyVector<double,3> ret = (m_nodes(0)->position() + m_nodes(1)->position() + m_nodes(2)->position());
	ret.scale(1.0/3.0);
	return ret;
}

/* ============================================================================== */
/* Calculate the unit normal */
TinyVector<double,3> Face::calculateUnitNormal() const
{
	TinyVector<double,3> r1 = m_nodes(0)->position();
	TinyVector<double,3> r2 = m_nodes(1)->position();
	TinyVector<double,3> r3 = m_nodes(2)->position();
	TinyVector<double,3> ret = CrossProduct(r2-r1, r3-r2);
	ret.scale(1.0/ret.norm());
	return ret;
}
/* ============================================================================== */
/* Debugging  */



// — stretch —
double Face::stretchingEnergy() const {
    double ρ      = stretchingEnergyContentDensity();
    // double factor = std::pow(m_adjust2, -6)
    //               * m_area
    //               * m_thickness
    //               * m_adjust1;
	double invAdjust2_6 = 1.0 / (m_adjust2 * m_adjust2 * m_adjust2 *
		m_adjust2 * m_adjust2 * m_adjust2
	);

	// double base = m_thickness * m_area * m_adjust1;

	double factor = invAdjust2_6 * m_area * m_thickness * m_adjust1;
	// std::cout << "[DEBUG Face] ρ=" << ρ
    //           << " area="    << m_area
    //           << " t="       << m_thickness
    //           << " a1="      << m_adjust1
    //           << " a2="      << m_adjust2
    //           << " factor="  << factor
    //           << "\n";
    return ρ * factor;
}

// — bend —
double Face::bendingEnergy() const {
    double β      = bendingEnergyContentDensity();
    // double factor = std::pow(m_adjust2, -6)
    //               * m_area
    //               * std::pow(m_thickness * m_adjust1, 3);
	double invAdjust2_6 = 1.0 / (m_adjust2 * m_adjust2 * m_adjust2 *
		m_adjust2 * m_adjust2 * m_adjust2
	);

	double base = m_thickness  * m_adjust1;

	double factor = invAdjust2_6 * m_area * base * base * base;
    // std::cout << "[DEBUG Face] β=" << β
    //           << " area="       << m_area
    //           << " t³a1³="      << std::pow(m_thickness * m_adjust1,3)
    //           << " a2="         << m_adjust2
    //           << " factor="     << factor
    //           << "\n";
    return β * factor;
}

// — connect —
double Face::connectionEnergy() const {
	double γ = connectionEnergyContentDensity();

	// Compute η_eff = η * adjust3
	double eta_eff = m_eta * m_adjust3;

	// Compute factor: η_eff² * h * a1 / a2⁶
	double invAdjust2_6 = 1.0 / (m_adjust2 * m_adjust2 * m_adjust2 *
	                             m_adjust2 * m_adjust2 * m_adjust2);

	double factor = eta_eff * eta_eff
	              * m_thickness * m_adjust1
	              * invAdjust2_6 * m_area;

	return γ * factor;
}
/* ============================================================================== */
/* Calculate stretching energy density */
double Face::stretchingEnergyContentDensity() const
{
	/* Calculate the 2D metric a */
	TinyVector<double,3> r1 = m_nodes(0)->position();
	TinyVector<double,3> r2 = m_nodes(1)->position();
	TinyVector<double,3> r3 = m_nodes(2)->position();
	TinyVector<double,3> dr12 = r2 - r1;
	TinyVector<double,3> dr23 = r3 - r2;
	TinyVector<double,3> dr31 = r1 - r3;
	TinyVector<double,3> l2;
	l2(0) = innerProduct(dr23,dr23);
	l2(1) = innerProduct(dr31,dr31);
	l2(2) = innerProduct(dr12,dr12);
	TinyVector<double,3> acomp = luSolve(m_Amatrix,m_Apivot,l2);
	TinyMatrix<double,2> a;
	a(0,0) = acomp(0);
	a(0,1) = acomp(1);
	a(1,0) = acomp(1);
	a(1,1) = acomp(2);

	/* Adjust abar with adjustment parameter */
	TinyMatrix<double,2> abarAdj = m_abar;
	abarAdj.scale(m_adjust2);
	abarAdj(0,0) += 1.0 - m_adjust2;
	abarAdj(1,1) += 1.0 - m_adjust2;
	TinyMatrix<double,2> invabarAdj = abarAdj.inverse();

	/* Calculate inv(abar)(a - abar) */
	TinyMatrix<double,2> tmp  = invabarAdj*(a - abarAdj);

	// return (m_lambda*tmp.trace()*tmp.trace() + m_mu*(tmp*tmp).trace());
	double E = (m_lambda * tmp.trace() * tmp.trace() + m_mu * (tmp*tmp).trace());
// std::cout << "m_lambda = " << m_lambda
//           << ", m_mu = " << m_mu
//           << ", thickness = " << m_thickness
//           << ", area = " << m_area
//           << ", adjust1 = " << m_adjust1
//           << ", adjust2 = " << m_adjust2
//           << ", E = " << E
//           << std::endl;
	// std::cout
    //   << "[DEBUG] m_lambda=" << m_lambda
    //   << ", m_mu="     << m_mu
    //   << ", E="        << E
    //   << "\n";
	return E;
}

/* ============================================================================== */
/* Calculate bending energy density */
// double Face::bendingEnergyContentDensity() const
// {
// 	/* Save the position of the Face */
// 	TinyVector<double,3> my_position = position();
// 	TinyVector<double,3> unitnormal = calculateUnitNormal();
// 	TinyVector<double,6> rhs;
// 	for (int nodeindex=0; nodeindex<6; nodeindex++)
// 	{
// 		if (m_nodes(nodeindex) != NULL)
// 		{
// 			TinyVector<double,3> dr = m_nodes(nodeindex)->position() - my_position;
// 			rhs(nodeindex) = innerProduct(dr,unitnormal);
// 		}
// 		else
// 		{
// 			const TinyVector<double,2>  *ptr = NULL;
// 			if (nodeindex==3)      ptr = &m_nodeconnector4;
// 			else if (nodeindex==4) ptr = &m_nodeconnector5;
// 			else                   ptr = &m_nodeconnector6;

// 			rhs(nodeindex) = m_bbar(0,0) * (*ptr)(0) * (*ptr)(0) +
// 				2*m_bbar(0,1)* (*ptr)(0) * (*ptr)(1) +
// 				m_bbar(1,1)* (*ptr)(1) * (*ptr)(1);
// 		}
// 	}

// 	TinyVector<double,6> bcomp = luSolve(m_Bmatrix, m_Bpivot, rhs);
// 	TinyMatrix<double,2> b;
// 	b(0,0) = bcomp(3);
// 	b(0,1) = bcomp(4);
// 	b(1,0) = bcomp(4);
// 	b(1,1) = bcomp(5);

// 	/* Adjust abar with adjustment parameter */
// 	TinyMatrix<double,2> abarAdj = m_abar;
// 	abarAdj.scale(m_adjust2);
// 	abarAdj(0,0) += 1.0 - m_adjust2;
// 	abarAdj(1,1) += 1.0 - m_adjust2;
// 	TinyMatrix<double,2> invabarAdj = abarAdj.inverse();

// 	/* Calculate inv(abar)(b - bbar) */
// 	TinyMatrix<double,2> tmp  = invabarAdj*(b - m_bbar);

// 	return (m_lambda*tmp.trace()*tmp.trace() + m_mu*(tmp*tmp).trace()) / 3;
// }


double Face::bendingEnergyContentDensity() const
{
    // get L,M,N
    auto LMN = this->LMN();
    TinyMatrix<double,2> b;
    b(0,0)= LMN(0);
    b(0,1)= LMN(1);
    b(1,0)= LMN(1);
    b(1,1)= LMN(2);

    // adjust reference curvature m_bbar
    TinyMatrix<double,2> abarAdj = m_abar;
    abarAdj.scale(m_adjust2);
    abarAdj(0,0) += 1.0 - m_adjust2;
    abarAdj(1,1) += 1.0 - m_adjust2;
    TinyMatrix<double,2> invA = abarAdj.inverse();

    TinyMatrix<double,2> tmp = invA * (b - m_bbar);
    return (m_lambda*tmp.trace()*tmp.trace() + m_mu*(tmp*tmp).trace())/3.0;
}

/* ============================================================================== */
/* Calculate connection energy density */
#include <iostream>    // for std::cout
double Face::connectionEnergyContentDensity() const {
    if (m_eta == 0.0)
        return 0.0;

    // Adjust reference metric abar to avoid numerical issues
    TinyMatrix<double,2> abarAdj = m_abar;
    abarAdj.scale(m_adjust2);
    abarAdj(0,0) += 1.0 - m_adjust2;
    abarAdj(1,1) += 1.0 - m_adjust2;

    // Elasticity tensor B^{μναβ} built from adjusted abar and Lamé parameters
    auto B = [&](int mu, int nu, int alpha, int beta) {
        return m_lambdaG * abarAdj(mu, nu) * abarAdj(alpha, beta)
             + m_muG * (abarAdj(mu, alpha) * abarAdj(nu, beta)
                      + abarAdj(mu, beta)  * abarAdj(nu, alpha));
    };

    // Compute actual and reference Christoffel symbols
    TinyMatrix<TinyMatrix<double,2>,2> Gamma     = computeConnection();
    const TinyMatrix<TinyMatrix<double,2>,2>& GammaBar = m_gammabar;

    // Perform the contraction:
    // W_conn = B^{μναβ} * (Γ^λ_{μν} - 𝛤̄^λ_{μν}) * (Γ^ρ_{αβ} - 𝛤̄^ρ_{αβ}) * 𝑎̄_{λρ}
	// TODO change
    // double result = 0.0;
    // for (int mu = 0; mu < 2; ++mu)
    // for (int nu = 0; nu < 2; ++nu)
    // for (int alpha = 0; alpha < 2; ++alpha)
    // for (int beta = 0; beta < 2; ++beta)
    // for (int lambda = 0; lambda < 2; ++lambda)
    // for (int rho = 0; rho < 2; ++rho) {
    //     double dGamma1 = Gamma(lambda, mu)(0, nu) - GammaBar(lambda, mu)(0, nu);
    //     double dGamma2 = Gamma(rho, alpha)(0, beta) - GammaBar(rho, alpha)(0, beta);
    //     result += B(mu, nu, alpha, beta) * dGamma1 * dGamma2 * abarAdj(lambda, rho);
    // }
	
	double result = 0.0;
	bool printedNeg = false;

	for (int μ = 0; μ < 2; ++μ) {
	for (int ν = 0; ν < 2; ++ν) {
		for (int α = 0; α < 2; ++α) {
		for (int β = 0; β < 2; ++β) {
			for (int λ = 0; λ < 2; ++λ) {
			for (int ρ = 0; ρ < 2; ++ρ) {
				double d1      = Gamma(λ,μ)(0,ν) - GammaBar(λ,μ)(0,ν);
				double d2      = Gamma(ρ,α)(0,β) - GammaBar(ρ,α)(0,β);
				double contrib = B(μ,ν,α,β) * d1 * d2 * abarAdj(λ,ρ);

				// if this single term is negative, print it once
				// if (!printedNeg && contrib < 0.0) {
				// std::cerr
				// 	<< "[DBG] NEG contrib at ("
				// 	<< μ<<","<<ν<<","<<α<<","<<β<<","<<λ<<","<<ρ
				// 	<< ") = " << contrib << "\n";
				// printedNeg = true;
				

				result += contrib;
			}
			}
		}
		}
	}
	}
	// after the loops:
	// std::cerr << "[DBG] full raw_result = " << result << "\n";

	// only print final if still negative
	if (result < 0.0) {
	std::cerr << "[DBG] FINAL raw_result = " << result << "\n";
	}





    if (result < 0.0) {
        std::cerr
          << "[DBG] NEGATIVE raw_result = " << result
          << ", m_thickness=" << m_thickness
          << ", m_eta="      << m_eta
          << ", m_adjust3="  << m_adjust3
          << "\n";
    }
    // Final result: h * η² * adjust3 * contraction

	// double factor  = 0.5 * m_thickness * m_eta * m_eta * m_adjust3 * result;

    return result;
}
/* ============================================================================== */
/* Calculate the first fundamental form */
TinyVector<double,3> Face::EFG() const
{
	/* Calculate the 2D metric a */
	TinyVector<double,3> r1 = m_nodes(0)->position();
	TinyVector<double,3> r2 = m_nodes(1)->position();
	TinyVector<double,3> r3 = m_nodes(2)->position();
	TinyVector<double,3> dr12 = r2 - r1;
	TinyVector<double,3> dr23 = r3 - r2;
	TinyVector<double,3> dr31 = r1 - r3;
	TinyVector<double,3> l2;
	l2(0) = innerProduct(dr23,dr23);
	l2(1) = innerProduct(dr31,dr31);
	l2(2) = innerProduct(dr12,dr12);
	TinyVector<double,3> acomp = luSolve(m_Amatrix,m_Apivot,l2);
	// std::cout << "Face " << this << " a = " << acomp << std::endl;

	return acomp;

}


/* ============================================================================== */
/* Helper: build the 2×2 metric a from EFG */
TinyMatrix<double,2> Face::computeMetric() const
{
    // EFG() returns {E, F, G} for this face
    TinyVector<double,3> ef = EFG();
    TinyMatrix<double,2> a;
    a(0,0) = ef(0);   // E
    a(0,1) = ef(1);   // F
    a(1,0) = ef(1);   // F
    a(1,1) = ef(2);   // G
    return a;
}
/* ============================================================================== */
/* computeMetricDerivatives */
// std::pair<TinyMatrix<double,2>, TinyMatrix<double,2>>
// Face::computeMetricDerivatives() const
// {
//     // ——— DEBUG: print connection stiffness params instead of da/du, da/dv ———
//     // std::cout
//     //   << "DEBUG Face @ (" << m_coordinates(0) << ", " << m_coordinates(1) << "): "
//     //   << "lambdaG = " << m_lambdaG
//     //   << ", muG = "     << m_muG
//     //   << std::endl;
//     // —————————————————————————————————————————————————————————————————————————

//     // 1) make sure neighbors exist
	
//     // assert(m_faces(0) && m_faces(1) && m_faces(2));
// 	{
//     // if any neighbor is missing, just return zeros
//     if (!m_faces(0) || !m_faces(1) || !m_faces(2)) {
// 		std::cerr << "Warning: Missing neighbor face(s) detected!" << std::endl;

//         return { TinyMatrix<double,2>(), TinyMatrix<double,2>() };
//     }
// 	}

//     // 2) pull the metric a on the three neighbor faces
//     TinyMatrix<double,2> N0 = m_faces(0)->computeMetric();
//     TinyMatrix<double,2> N1 = m_faces(1)->computeMetric();
//     TinyMatrix<double,2> N2 = m_faces(2)->computeMetric();

//     // 3) compute parametric differences Δu and Δv
//     double du = m_faces(0)->coordinates()(0)
//               - m_faces(1)->coordinates()(0);
//     double dv = m_faces(0)->coordinates()(1)
//               - m_faces(2)->coordinates()(1);
// 	if (fabs(du) < 1e-12 || fabs(dv) < 1e-12) {
// 		// handle boundary or degenerate case…
// 		return { TinyMatrix<double,2>(), TinyMatrix<double,2>() };
// 	}

//     // 4) finite‐difference the metric
//     TinyMatrix<double,2> da_du = (N0 - N1) * (1.0/du);
//     TinyMatrix<double,2> da_dv = (N0 - N2) * (1.0/dv);
// 	// std::cout << "du = " << du << ", dv = " << dv << std::endl;
// 	// if(fabs(du) < 1e-8 || fabs(dv) < 1e-8)
// 	// 	std::cerr << "Warning: Very small du/dv detected!" << std::endl;

//     return std::make_pair(da_du, da_dv);
// }

// /* computeMetricDerivatives - works well but cooses 2 neighbors*/
// std::pair<TinyMatrix<double,2>, TinyMatrix<double,2>>
// Face::computeMetricDerivatives() const
// {
//     const double tol = 1e-8;

//     // 1) pull neighbor metrics (if they exist)
//     TinyMatrix<double,2> N0, N1, N2;
//     if (m_faces(0)) N0 = m_faces(0)->computeMetric();
//     if (m_faces(1)) N1 = m_faces(1)->computeMetric();
//     if (m_faces(2)) N2 = m_faces(2)->computeMetric();

//     // 2) this face’s own metric
//     TinyMatrix<double,2> a0 = computeMetric();

//     // 3) uv coords of this face and its neighbors
//     double u0  = m_coordinates(0),      v0  = m_coordinates(1);
//     double u0p = m_faces(0) ? m_faces(0)->coordinates()(0) : 0.0;
//     double u1p = m_faces(1) ? m_faces(1)->coordinates()(0) : 0.0;
//     double v0p = m_faces(0) ? m_faces(0)->coordinates()(1) : 0.0;
//     double v2p = m_faces(2) ? m_faces(2)->coordinates()(1) : 0.0;

//     // 4) allocate outputs
//     TinyMatrix<double,2> da_du, da_dv;

//     // —— ∂a/∂u using neighbors 0 & 1 —— 
//     if (m_faces(0) && m_faces(1)) {
//         double du = u0p - u1p;
//         if (fabs(du) > tol)
//             da_du = (N0 - N1) * (1.0/du);
//         else
//             da_du = TinyMatrix<double,2>();
//     }
//     else if (m_faces(0)) {
//         double du = u0p - u0;
//         if (fabs(du) > tol)
//             da_du = (N0 - a0) * (1.0/du);
//         else
//             da_du = TinyMatrix<double,2>();
//     }
//     else if (m_faces(1)) {
//         double du = u0 - u1p;
//         if (fabs(du) > tol)
//             da_du = (a0 - N1) * (1.0/du);
//         else
//             da_du = TinyMatrix<double,2>();
//     }
//     else {
//         da_du = TinyMatrix<double,2>();  // no neighbors ⇒ zero
//     }

//     // —— ∂a/∂v using neighbors 0 & 2 —— 
//     if (m_faces(0) && m_faces(2)) {
//         double dv = v0p - v2p;
//         if (fabs(dv) > tol)
//             da_dv = (N0 - N2) * (1.0/dv);
//         else
//             da_dv = TinyMatrix<double,2>();
//     }
//     else if (m_faces(0)) {
//         double dv = v0p - v0;
//         if (fabs(dv) > tol)
//             da_dv = (N0 - a0) * (1.0/dv);
//         else
//             da_dv = TinyMatrix<double,2>();
//     }
//     else if (m_faces(2)) {
//         double dv = v0 - v2p;
//         if (fabs(dv) > tol)
//             da_dv = (a0 - N2) * (1.0/dv);
//         else
//             da_dv = TinyMatrix<double,2>();
//     }
//     else {
//         da_dv = TinyMatrix<double,2>();  // no neighbors ⇒ zero
//     }

//     // 5) return the two partial derivatives
//     return { da_du, da_dv };
// }
// std::pair<TinyMatrix<double,2>, TinyMatrix<double,2>>
// Face::computeMetricDerivatives() const
// {
//     // Collect available neighbors and their UVs
//     struct Neighbor {
//         TinyMatrix<double,2> metric;
//         double du, dv;
//         double weight;
//     };
//     std::vector<Neighbor> neighbors;

//     const TinyMatrix<double,2> a0 = computeMetric();
//     const double u0 = m_coordinates(0), v0 = m_coordinates(1);

//     for (int ni = 0; ni < 3; ++ni) {
//         if (!m_faces(ni)) continue;
//         TinyMatrix<double,2> aj = m_faces(ni)->computeMetric();
//         double uj = m_faces(ni)->coordinates()(0);
//         double vj = m_faces(ni)->coordinates()(1);

//         double du = uj - u0;
//         double dv = vj - v0;
//         double dist2 = du*du + dv*dv;
//         if (dist2 < 1e-12) continue;

//         // Weight: inverse squared distance (safe, regularizes for closer neighbors)
//         double w = 1.0 / dist2;

//         neighbors.push_back({aj, du, dv, w});
//     }

//     // If no neighbors, return zero
//     if (neighbors.empty())
//         return {TinyMatrix<double,2>(), TinyMatrix<double,2>()};

//     // Least-squares: Solve W*Δa = W*[du, dv] * [da_du; da_dv]
//     // For each metric component, solve independently (since TinyMatrix is small)

//     // Each metric is symmetric 2x2, we can unroll as [a00, a01, a11]
//     // We'll fit each component separately

//     TinyMatrix<double,2> da_du, da_dv;
//     for (int m=0; m<2; ++m) for (int n=0; n<=m; ++n) { // (0,0), (0,1), (1,1)
//         // Set up least-squares
//         double Sww = 0, Swu = 0, Swv = 0, Swuu = 0, Swvv = 0, Swuv = 0;
//         double Swa = 0, Swau = 0, Swav = 0;

//         for (const auto& N : neighbors) {
//             double da = N.metric(m,n) - a0(m,n);
//             double w = N.weight;
//             Sww  += w;
//             Swu  += w * N.du;
//             Swv  += w * N.dv;
//             Swuu += w * N.du * N.du;
//             Swvv += w * N.dv * N.dv;
//             Swuv += w * N.du * N.dv;
//             Swa  += w * da;
//             Swau += w * da * N.du;
//             Swav += w * da * N.dv;
//         }

//         // Build 2x2 system:
//         // [Swuu  Swuv][X] = [Swau]
//         // [Swuv  Swvv][Y]   [Swav]
//         double det = Swuu*Swvv - Swuv*Swuv;
//         double x = 0, y = 0;
//         if (fabs(det) > 1e-14) {
//             x = ( Swau*Swvv - Swav*Swuv) / det;
//             y = (Swuu*Swav - Swuv*Swau) / det;
//         }
//         // Place into the output derivatives
//         da_du(m,n) = x;
//         da_dv(m,n) = y;
//         if (m != n) {
//             da_du(n,m) = x;
//             da_dv(n,m) = y;
//         }
//     }
//     return {da_du, da_dv};
// }

std::pair<TinyMatrix<double,2>, TinyMatrix<double,2>>
Face::computeMetricDerivatives() const
{
    TinyVector<double,2> uv0;
    uv0(0) = m_coordinates(0);
    uv0(1) = m_coordinates(1);

    TinyMatrix<double,2> a0 = computeMetric();
    // std::cout << "[DEBUG] a = " << a0 << std::endl;


    struct Entry {
        TinyVector<double,2> duv;  // normalized UV direction
        TinyMatrix<double,2> Da;   // directional derivative of metric
    };

    std::vector<Entry> entries;

    for (int i = 0; i < 3; ++i) {
        Face* f = m_faces(i);
        if (!f) continue;

        // UV direction
        TinyVector<double,2> uv_i;
        uv_i(0) = f->coordinates()(0);
        uv_i(1) = f->coordinates()(1);

        TinyVector<double,2> duv = uv_i - uv0;
        double norm = std::sqrt(duv(0)*duv(0) + duv(1)*duv(1));
        if (norm < 1e-6) continue;

        // Normalize manually
        TinyVector<double,2> duv_norm;
        duv_norm(0) = duv(0) / norm;
        duv_norm(1) = duv(1) / norm;

        // Distance in 3D
        TinyVector<double,3> p0 = position();
        TinyVector<double,3> pi = f->position();
        TinyVector<double,3> dp = pi - p0;
        double dist = std::sqrt(dp(0)*dp(0) + dp(1)*dp(1) + dp(2)*dp(2));
        if (dist < 1e-6) continue;

        // Directional derivative estimate
        TinyMatrix<double,2> Da = f->computeMetric() - a0;
        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c)
                Da(r, c) = Da(r, c) / dist;

        entries.push_back({duv_norm, Da});
    }

    if (entries.size() < 2)
        return {TinyMatrix<double,2>(), TinyMatrix<double,2>()};

    // Weighted average: weights normalized
    double weight_sum = 0.0;
    for (const auto& e : entries)
        weight_sum += e.duv(0)*e.duv(0) + e.duv(1)*e.duv(1);

    TinyMatrix<double,2> da_du, da_dv;
    da_du = TinyMatrix<double,2>();
    da_dv = TinyMatrix<double,2>();

    for (const auto& e : entries) {
        double u = e.duv(0);
        double v = e.duv(1);
        double w = (u*u + v*v) / weight_sum;

        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c) {
                da_du(r, c) += w * u * e.Da(r, c);
                da_dv(r, c) += w * v * e.Da(r, c);
            }
    }
	// std::cout << "[DEBUG] ∂a/∂u = " << da_du << "\n";
	// std::cout << "[DEBUG] ∂a/∂v = " << da_dv << "\n";


    return {da_du, da_dv};
}


/* ============================================================================== */
/* Calculate the second fundamental form */
// TinyVector<double,3> Face::LMN() const
// {
//     // collect all non‐NULL nodes
//     std::vector<int> valid;
//     for(int i=0; i<6; ++i)
//         if(m_nodes(i)) valid.push_back(i);

//     if(valid.size()<3)
//         return TinyVector<double,3>(); // not enough to fit

//     // face center and normal
//     auto  P0 = position();
//     auto  n0 = calculateUnitNormal();

//     // build least squares: fit z ≈ L x^2 + 2M x y + N y^2
//     const int Np = valid.size();
//     std::vector<std::array<double,3>> B(Np);
//     std::vector<double>                 rhs(Np);
//     for(int k=0; k<Np; ++k){
//         int idx = valid[k];
//         TinyVector<double,2> c = (idx<3
//             ? m_nodes(idx)->coordinates() - m_coordinates
//             : (idx==3? m_nodeconnector4
//               : idx==4? m_nodeconnector5
//                        : m_nodeconnector6));
//         B[k][0] = c(0)*c(0);
//         B[k][1] = c(0)*c(1);
//         B[k][2] = c(1)*c(1);
//         rhs[k]  = innerProduct(m_nodes(idx)->position() - P0, n0);
//     }

//     // form normal equations BTB x = BTrhs
//     TinyMatrix<double,3> BTB;
//     TinyVector<double,3> BTrhs;
//     for(int i=0;i<3;++i){
//       BTrhs(i)=0; 
//       for(int j=0;j<3;++j) BTB(i,j)=0;
//       for(int k=0;k<Np;++k){
//         BTrhs(i) += B[k][i]*rhs[k];
//         for(int j=i;j<3;++j){
//           BTB(i,j) += B[k][i]*B[k][j];
//           BTB(j,i)  = BTB(i,j);
//         }
//       }
//     }

//     TinyVector<int,3> piv;
//     auto BTBcopy = BTB;
//     luDecompose(BTBcopy,piv);
//     return luSolve(BTBcopy,piv,BTrhs);
// }

/* ============================================================================== */
/* ============================================================================== */
/* Calculate the second fundamental form */
/* ============================================================================== */
/* Calculate the second fundamental form, handling full vs. partial neighbors */
TinyVector<double,3> Face::LMN() const
{
    // 1) check if *all* 6 m_nodes are non‐NULL ⇒ use your original 6×6 solve
    bool haveAll = true;
    for(int i = 0; i < 6; ++i){
        if (!m_nodes(i)) { haveAll = false; break; }
    }
    if (haveAll) {
        // — exactly your original code —
        TinyVector<double,3> P0 = position();
        TinyVector<double,3> n0 = calculateUnitNormal();
        TinyVector<double,6> rhs;
        for (int i = 0; i < 6; ++i) {
            TinyVector<double,3> dr = m_nodes(i)->position() - P0;
            rhs(i) = innerProduct(dr, n0);
        }
        TinyVector<double,6> bcomp = luSolve(m_Bmatrix, m_Bpivot, rhs);
        TinyVector<double,3> Lmn;
        Lmn(0) = bcomp(3);
        Lmn(1) = bcomp(4);
        Lmn(2) = bcomp(5);
        return Lmn;
    }

    // 2) otherwise: fall back to a 3×3 least‐squares on [x², 2xy, y²]
    struct Sample { double x,y,h; };
    std::vector<Sample> samp;
    auto P0 = position();
    auto n0 = calculateUnitNormal();

    // collect only the real nodes
    for(int i=0; i<6; ++i){
        if (!m_nodes(i)) continue;
        TinyVector<double,2> c = (i<3
          ? m_nodes(i)->coordinates() - m_coordinates
          : (i==3 ? m_nodeconnector4
           : i==4 ? m_nodeconnector5
                  : m_nodeconnector6));
        double h = innerProduct(m_nodes(i)->position() - P0, n0);
        samp.push_back({ c(0), c(1), h });
    }
    // need at least 3 points to fit 3 unknowns (L,2M,N)
    if (samp.size() < 3)
        return TinyVector<double,3>();  

    // build the 3×3 normal equations A·c = b
    TinyMatrix<double,3> A;  A.setToZero();
    TinyVector<double,3> b;  b.setToZero();

    for(auto &s : samp){
        double B0 = s.x*s.x;
        double B1 =     2*s.x*s.y;
        double B2 = s.y*s.y;
        // accumulate A
        A(0,0) += B0*B0;  A(0,1) += B0*B1;  A(0,2) += B0*B2;
                          A(1,1) += B1*B1;  A(1,2) += B1*B2;
                                             A(2,2) += B2*B2;
        // accumulate rhs
        b(0) += B0 * s.h;
        b(1) += B1 * s.h;
        b(2) += B2 * s.h;
    }
    // symmetrize
    A(1,0) = A(0,1);
    A(2,0) = A(0,2);
    A(2,1) = A(1,2);

    // solve 3×3
    TinyVector<int,3> piv;
    auto Ac = A;
    luDecompose(Ac, piv);
    TinyVector<double,3> c3 = luSolve(Ac, piv, b);

    // extract L, M, N
    TinyVector<double,3> Lmn;
    Lmn(0) = c3(0);        // L
    Lmn(1) = c3(1) * 0.5;  // 2M → M
    Lmn(2) = c3(2);        // N
    return Lmn;
}

/* ============================================================================== */
/* Calculate the third fundamental form c = b a^{-1} b */
TinyVector<double,3> Face::ABC() const
{
    // First fundamental form
    TinyVector<double,3> ef = EFG();

    TinyMatrix<double,2> a;
    a(0,0) = ef(0);
    a(0,1) = ef(1);
    a(1,0) = ef(1);
    a(1,1) = ef(2);

    // Second fundamental form
    TinyVector<double,3> lmn = LMN();

    TinyMatrix<double,2> b;
    b(0,0) = lmn(0);
    b(0,1) = lmn(1);
    b(1,0) = lmn(1);
    b(1,1) = lmn(2);

    // c = b a^{-1} b
    TinyMatrix<double,2> c = b * a.inverse() * b;

    TinyVector<double,3> abc;
    abc(0) = c(0,0);
    abc(1) = c(0,1);
    abc(2) = c(1,1);

    return abc;
}
/* ============================================================================== */


/* ==============================================================================  */
/* Calculate the Gamma terms */
/*
   Compute the Levi–Civita connection Γᵏ_{ij} on this face from the current metric a and its derivatives.
   Uses:
     - a = computeMetric(u,v)
     - (da_du, da_dv) = computeMetricDerivatives(u,v)
   Formula:
     Γᵏ_{ij} = ½ a^{kℓ} ( ∂ᵢ a_{ℓj} + ∂ⱼ a_{ℓi} - ∂_ℓ a_{ij} )
*/
TinyMatrix<TinyMatrix<double,2>,2> Face::computeConnection() const {
	
    // 1) compute metric and its inverse
    TinyVector<double,3> ef = EFG();
    TinyMatrix<double,2> a;
    a(0,0) = ef(0);  a(0,1) = ef(1);
    a(1,0) = ef(1);  a(1,1) = ef(2);
    TinyMatrix<double,2> inva = a.inverse();

    // 2) get the two derivative‐matrices
    TinyMatrix<double,2> da_du, da_dv;
    std::tie(da_du, da_dv) = computeMetricDerivatives();

    // 3) extract inva entries into locals
    const double i00 = inva(0,0),
                 i01 = inva(0,1),
                 i11 = inva(1,1);

    // 4) alias da_du / da_dv to avoid the (i==0? …) branches inside the loop
    const TinyMatrix<double,2>* d[2] = { &da_du, &da_dv };

    TinyMatrix<TinyMatrix<double,2>,2> Gamma;

    // now only three loops: k,i,j
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < 2; ++i) {
        const TinyMatrix<double,2>& di = *d[i];
        for (int j = 0; j < 2; ++j) {
          const TinyMatrix<double,2>& dj = *d[j];

          // unroll l=0 and l=1 by hand:
          //   term(l) = ∂ᵢa_{ℓj} + ∂ⱼa_{ℓi} – ∂_ℓa_{ij}
          double term0 = di(0,j) + dj(0,i) - (*d[0])(i,j);
          double term1 = di(1,j) + dj(1,i) - (*d[1])(i,j);

          // sum = ∑ₗ a^{kℓ} * term(l)
          double sum = (k==0 ? (i00*term0 + i01*term1)
                             : (i01*term0 + i11*term1));

          Gamma(k,i)(0,j) = 0.5 * sum;
        }
      }
    }

    return Gamma;
}
// Note: the original code had a commented‐out version of this function
// that used a different indexing scheme for Gamma, which was slow.
// The current version uses the correct indexing as per the formula above.
// TinyMatrix<TinyMatrix<double,2>,2> Face::computeConnection() const {
//     TinyVector<double,3> ef = EFG();  // First fundamental form

//     TinyMatrix<double,2> a;
//     a(0,0) = ef(0);  a(0,1) = ef(1);
//     a(1,0) = ef(1);  a(1,1) = ef(2);
//     TinyMatrix<double,2> inva = a.inverse();

//     TinyMatrix<double,2> da_du, da_dv;
//     std::tie(da_du, da_dv) = computeMetricDerivatives();  // Per face

//     TinyMatrix<TinyMatrix<double,2>,2> Gamma;

//     for (int k = 0; k < 2; ++k) {
//         for (int i = 0; i < 2; ++i) {
//             for (int j = 0; j < 2; ++j) {
//                 double sum = 0.0;
//                 for (int l = 0; l < 2; ++l) {
//                     const TinyMatrix<double,2>& d_i = (i == 0 ? da_du : da_dv);
//                     const TinyMatrix<double,2>& d_j = (j == 0 ? da_du : da_dv);
//                     const TinyMatrix<double,2>& d_l = (l == 0 ? da_du : da_dv);

//                     double term = d_i(l,j) + d_j(l,i) - d_l(i,j);
//                     sum += inva(k,l) * term;
//                 }
//                 Gamma(k,i)(0,j) = 0.5 * sum;  // ✅ specify both indices
// 				// Gamma(k,i)(j,0) = 0.5 * sum; 

//             }
//         }
//     }

//     return Gamma;
// }
/* ============================================================================== */
/* Calculate the energy gradient */
void Face::setForce() {
    double ep = 1.e-6;
    for (int i=0; i<6; i++) {
        if (m_nodes(i) != NULL) {
            for (int comp=0; comp<3; comp++) {
                m_nodes(i)->position(comp) += ep;
                double Eplus = energy();
                m_nodes(i)->position(comp) -= 2*ep;
                double Eminus = energy();
                double grad = 0.5*(Eplus-Eminus)/ep;

                // Diagnostic print statement here:
                // std::cout << "Gradient component [" << i << "][" << comp << "] = " << grad << std::endl;

                m_nodes(i)->force(comp) += grad;
                m_nodes(i)->position(comp) += ep;
            }
        }
    }
}

/* ============================================================================== */
/* NonEuclideanShell NonEuclideanShell NonEuclideanShell NonEuclideanShell    */
/* ============================================================================== */
/* Constructor */
NonEuclideanShell::NonEuclideanShell(const std::string &a_nodesFileName,
									 const std::string &a_facesFileName) :
m_nodes(),
m_faces(),
m_verbosity(2)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::NonEuclideanShell()");

	int		numberNodes, numberFaces;
	int		n1, n2, n3, n4, n5, n6, f1, f2, f3;
	int		index;
	double	u,v;
	int		isfixed;

	/* Take care of Nodes */
	TextFileHandle nodesFileHandle(a_nodesFileName,FileHandle::OPEN_RD);
	nodesFileHandle.read(numberNodes);
	m_nodes = Vector<Node*>(numberNodes);

	if (m_verbosity>1)
		std::cout << "NonEuclideanShell::NonEuclideanShell()   Number of Nodes = " << numberNodes << std::endl;
    // ——— DEBUG: collect counts of each isfixed value ———
    std::map<int,int> flagCounts;

	for (int i=0; i<numberNodes; i++)
	{
		nodesFileHandle.read(index);
		nodesFileHandle.read(u);
		nodesFileHandle.read(v);
		nodesFileHandle.read(isfixed);
		m_nodes(i)  = new Node();
        // ——— DEBUG: print each one as it arrives ———
        // std::cerr << "[DEBUG] Read vertex " << index
        //           << " → u=" << u << ", v=" << v
        //           << ", isfixed=" << isfixed << "\n";
        // tally for summary
        flagCounts[isfixed]++;
		m_nodes(i)->coordinates(0) = u;
		m_nodes(i)->coordinates(1) = v;
		m_nodes(i)->fixed() = isfixed;
	}
	nodesFileHandle.close();
	  // ——— DEBUG: print a summary of all flags seen ———
    std::cerr << "[DEBUG] isfixed value summary:\n";
    for (auto &kv : flagCounts) {
        std::cerr << "   " << kv.first << " → " << kv.second << " vertices\n";
    }

	/* Take care of Faces */
	TextFileHandle facesFileHandle(a_facesFileName,FileHandle::OPEN_RD);
	facesFileHandle.read(numberFaces);
	m_faces = Vector<Face*>(numberFaces);

	if (m_verbosity>1)
		std::cout << "NonEuclideanShell::NonEuclideanShell()   Number of Faces = " << numberFaces << std::endl;

	for (int i=0; i<numberFaces; i++)
	{
		facesFileHandle.read(index);
		facesFileHandle.read(n1);
		facesFileHandle.read(n2);
		facesFileHandle.read(n3);
		facesFileHandle.read(n4);
		facesFileHandle.read(n5);
		facesFileHandle.read(n6);
		facesFileHandle.read(f1);
		facesFileHandle.read(f2);
		facesFileHandle.read(f3);

		TinyVector<Node*, 6> nodeVec;
		nodeVec(0) = m_nodes(n1);
		nodeVec(1) = m_nodes(n2);
		nodeVec(2) = m_nodes(n3);
		nodeVec(3) = (n4==-1) ? NULL : m_nodes(n4);
		nodeVec(4) = (n5==-1) ? NULL : m_nodes(n5);
		nodeVec(5) = (n6==-1) ? NULL : m_nodes(n6);

		m_faces(i) = new Face(nodeVec);
	}
	facesFileHandle.close();

	/* Take care of Faces */
	TextFileHandle facesFileHandle2(a_facesFileName,FileHandle::OPEN_RD);
	facesFileHandle2.read(numberFaces);

	for (int i=0; i<numberFaces; i++)
	{
		facesFileHandle2.read(index);
		facesFileHandle2.read(n1);
		facesFileHandle2.read(n2);
		facesFileHandle2.read(n3);
		facesFileHandle2.read(n4);
		facesFileHandle2.read(n5);
		facesFileHandle2.read(n6);
		facesFileHandle2.read(f1);
		facesFileHandle2.read(f2);
		facesFileHandle2.read(f3);

		TinyVector<Face*, 3> faceVec;
		faceVec(0) = (f1==-1) ? NULL : m_faces(f1);
		faceVec(1) = (f2==-1) ? NULL : m_faces(f2);
		faceVec(2) = (f3==-1) ? NULL : m_faces(f3);

		m_faces(i)->assignNeighbors(faceVec);
	}
	facesFileHandle2.close();

}

/* ============================================================================== */
/* Destructor */
NonEuclideanShell::~NonEuclideanShell()
{
	for (int i=0; i<m_nodes.length(); i++) delete m_nodes(i);
	for (int i=0; i<m_faces.length(); i++) delete m_faces(i);
}

/* ============================================================================== */

std::vector<int> NonEuclideanShell::getVertexNeighbors(int v) const {
    std::vector<int> nbrs;
    // the Node* we’re interested in
    Node* target = m_nodes(v);

    // for every face, if 'target' is one of its 3 vertices, add the other two
    for (int fi = 0; fi < m_faces.length(); ++fi) {
        Face* f = m_faces(fi);
        for (int ei = 0; ei < 3; ++ei) {
            if (f->node(ei) == target) {
                // grab the other two corner indices
                for (int k = 1; k <= 2; ++k) {
                    Node* nbrNode = f->node((ei + k) % 3);
                    // find its index in m_nodes
                    for (int ni = 0; ni < m_nodes.length(); ++ni) {
                        if (m_nodes(ni) == nbrNode) {
                            if (std::find(nbrs.begin(), nbrs.end(), ni) == nbrs.end())
                                nbrs.push_back(ni);
                            break;
                        }
                    }
                }
            }
        }
    }
    return nbrs;
}


/* ============================================================================== */
/* Set reference forms and parameters */
void NonEuclideanShell::setParameters(double (*f_thickness)(double, double),
									  double (*f_lambda)(double, double),
									  double (*f_mu)(double, double),
									  TinyMatrix<double,2> (*f_abar)(double, double),
									  TinyMatrix<double,2> (*f_bbar)(double, double),
									  TinyMatrix<double,2> (*f_cbar)(double, double),
									  TinyVector<double,3> (*f_Pos0)(double, double),
									  TinyMatrix<TinyMatrix<double,2>,2> (*f_gammabar)(double, double),
									  double eta,
                                      double lambdaG,
                                      double muG)
{	
	m_eta = eta;

	for (int i=0; i<m_faces.length(); i++)
	{
		double u         = m_faces(i)->coordinates(0);
		double v         = m_faces(i)->coordinates(1);
		double thickness = f_thickness(u,v);
		double lambda    = f_lambda(u,v);
		double mu        = f_mu(u,v);
		TinyMatrix<double,2> abar      = f_abar(u,v);
		TinyMatrix<double,2> bbar      = f_bbar(u,v);
		TinyMatrix<double,2> cbar      = f_cbar(u,v);
        TinyMatrix<TinyMatrix<double,2>,2> gammabar = f_gammabar(u, v);

		m_faces(i)->initialize(thickness, lambda, mu, abar, bbar, cbar, gammabar, eta, lambdaG, muG);
    }

    // --- initialize node positions (unchanged) ---
    for (int i = 0; i < m_nodes.length(); ++i) 
	{
        double u = m_nodes(i)->coordinates()(0);
        double v = m_nodes(i)->coordinates()(1);
        // auto   P = f_Pos0(u, v);
        // m_nodes(i)->position()(0) = P(0);
        // m_nodes(i)->position()(1) = P(1);
        // m_nodes(i)->position()(2) = P(2);
		m_nodes(i)->position(0) = f_Pos0(u,v)(0);
		m_nodes(i)->position(1) = f_Pos0(u,v)(1);
		m_nodes(i)->position(2) = f_Pos0(u,v)(2);
    }

    // --- fix offsets (unchanged) ---
    for (int i = 0; i < m_nodes.length(); ++i) {
        int j = m_nodes(i)->fixed();
        if (j >= 0) {
            m_nodes(i)->offset() =
                m_nodes(i)->position() - m_nodes(j)->position();
        }
    }
}
/* ============================================================================== */
/* set the adjustment parameters */
// void NonEuclideanShell::setAdjust(double a_adjust1, double a_adjust2)
void NonEuclideanShell::setAdjust(double a_adjust1, double a_adjust2, double a_adjust3)

{
	for (int i=0; i<m_faces.length(); i++)
	{
		// m_faces(i)->setAdjust(a_adjust1, a_adjust2);
		m_faces(i)->setAdjust(a_adjust1, a_adjust2, a_adjust3);

	}
}

/* ============================================================================== */
/* Initialization: set initial configuration */
void NonEuclideanShell::defaultInitialization()
{
	for (int i=0; i<m_nodes.length(); i++)
	{
		m_nodes(i)->position(0) = m_nodes(i)->coordinates(0);
		m_nodes(i)->position(1) = m_nodes(i)->coordinates(1);
		m_nodes(i)->position(2) = 0.1*sin(m_nodes(i)->coordinates(0));
	}
}

/* ============================================================================== */
/* Restart from file */
void NonEuclideanShell::restart(BinaryFileHandle* a_fh)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::restart()");

	for (int i=0; i<m_nodes.length(); i++)
	{
		float X,Y,Z;
		a_fh->read(&X);
		a_fh->read(&Y);
		a_fh->read(&Z);
		m_nodes(i)->position(0) = X;
		m_nodes(i)->position(1) = Y;
		m_nodes(i)->position(2) = Z;
	}
}

/* ============================================================================== */
/* Initialization: set forces to zero */
void NonEuclideanShell::initializeForce()
{
	for (int i=0; i<m_nodes.length(); i++)
		m_nodes(i)->force().setToZero();
}
/* ============================================================================== */

double NonEuclideanShell::meanEdgeLength() const {
    int n = m_nodes.length();                                                 // total nodes :contentReference[oaicite:3]{index=3}
    double sumNN = 0.0;

    // For each node, find its closest neighbor
    for (int i = 0; i < n; ++i) {
        auto ri = m_nodes(i)->position();                                      // TinyVector<double,3>
        double best = std::numeric_limits<double>::infinity();

        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            auto rj = m_nodes(j)->position();
            double dx = ri(0) - rj(0);
            double dy = ri(1) - rj(1);
            double dz = ri(2) - rj(2);
            double d  = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d < best) best = d;
        }

        sumNN += best;
    }

    return sumNN / n;
}

/* ============================================================================== */
/* Calculate the stretching energy */
// double NonEuclideanShell::stretchingEnergy() const
// {
// 	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::stretchingEnergy()");
	
// 	double energy=0;
// 	for (int i=0; i<m_faces.length(); i++) {
// 		energy += m_faces(i)->stretchingEnergy();
// 		// std::cout << "Face " << i << " stretch energy: " << m_faces(i)->stretchingEnergy() << std::endl;
// 	}

// 	if (m_verbosity>1)
// 		std::cout << "NonEuclideanShell::stretchingEnergy()    = " << energy << std::endl;
// 	// std::cout << "NonEuclideanShell::stretchingEnergy()    = " << energy << std::endl;

// 	return energy;
// }
double NonEuclideanShell::stretchingEnergy() const {
    double energy = 0;
    for (int i = 0; i < m_faces.length(); ++i) {
        // double density = m_faces(i)->stretchingEnergyContentDensity();
        // double A       = m_faces(i)->area();
        double Eface   = m_faces(i)->stretchingEnergy();
        // std::cout << "[DEBUG] Face " << i
        //           << " ρ=" << density
        //           << "  A=" << A
        //           << " E=" << Eface
        //           << std::endl;
        energy += Eface;
    }
    // std::cout << "[DEBUG] Total stretching = " << energy << std::endl;
    return energy;
}


/* ============================================================================== */
/* Calculate the bending energy */
double NonEuclideanShell::bendingEnergy() const
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::bendingEnergy()");

	double energy=0;
	for (int i=0; i<m_faces.length(); i++)
		energy += m_faces(i)->bendingEnergy();

	if (m_verbosity>1)
		std::cout << "NonEuclideanShell::bendingEnergy()    = " << energy << std::endl;

	return energy;
}
/* ============================================================================== */
/* Calculate the connection (γ) energy */
double NonEuclideanShell::connectionEnergy() const
{
    if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::connectionEnergy()");

    double energy = 0.0;
    for (int i = 0; i < m_faces.length(); ++i)
        energy += m_faces(i)->connectionEnergy();

    if (m_verbosity>1)
        std::cout << "NonEuclideanShell::connectionEnergy()    = " << energy << std::endl;

    return energy;
}
// double NonEuclideanShell::connectionEnergy() const {
//     double sum = 0.0;
//     double maxEnergy = 0.0;
//     int maxIndex = -1;

//     for (int i = 0; i < m_faces.length(); i++) {
//         double e = m_faces(i)->connectionEnergy();
//         sum += e;
//         if (e > maxEnergy) {
//             maxEnergy = e;
//             maxIndex = i;
//         }
//     }

//     std::cout << "[DEBUG] Max connection energy = " << maxEnergy
//               << " at face index " << maxIndex << std::endl;

//     return sum;
// }

/* ============================================================================== */
/* Calculate the total energy */
double NonEuclideanShell::energy() const
{
    double energy;
    if (m_faces.length() > 0 && m_faces(0)->getLambdaG() == 0.0 && m_faces(0)->getMuG() == 0.0) {
    // std::cout << "[DEBUG] Skipping connection energy (lambdaG and muG are zero for all faces)\n";
    energy = stretchingEnergy() + bendingEnergy();
	} else {
    energy = stretchingEnergy() + bendingEnergy() + connectionEnergy();
	}


    if (m_verbosity > 1)
        std::cout << "NonEuclideanShell::energy()    = " << energy << std::endl;

    return energy;
}

/* ============================================================================== */
/* Calculate size of optimization problem */
int NonEuclideanShell::SizeOfOptimizationProblem() const
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::SizeOfOptimizationProblem()");

	int sop = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2) sop+=3;
	}

	return sop;
}

/* ============================================================================== */
/* Calculate the stretching forces */
void NonEuclideanShell::setForce()
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::setForce()");

	for (int i=0; i<m_faces.length(); i++)
		m_faces(i)->setForce();

}

/* ============================================================================== */
/* Dump results to file (binary) */
void NonEuclideanShell::DumpStateBinaryFormat(BinaryFileHandle* a_fh)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::DumpStateBinaryFormat()");

	for (int i=0; i<m_nodes.length(); i++)
	{
		float X = (float)m_nodes(i)->position(0);
		float Y = (float)m_nodes(i)->position(1);
		float Z = (float)m_nodes(i)->position(2);
		a_fh->write(&X);
		a_fh->write(&Y);
		a_fh->write(&Z);
	}
}
/* ============================================================================== */
/* Dump results to file (text) */
void NonEuclideanShell::DumpStateTextFormat(TextFileHandle* a_node, TextFileHandle* a_face, int counter)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::DumpStateTextFormat()");

	for (int i=0; i<m_nodes.length(); i++)
	{
		double X = m_nodes(i)->position(0);
		double Y = m_nodes(i)->position(1);
		double Z = m_nodes(i)->position(2);
		a_node->write(i);
		a_node->tab();
		a_node->write(X);
		a_node->tab();
		a_node->write(Y);
		a_node->tab();
		a_node->write(Z);
		a_node->tab();
		a_node->write(counter);
		a_node->newline();
	}

	for (int i=0; i<m_faces.length(); i++)
	{
		// double Es = m_faces(i)->stretchingEnergyContentDensity();
		// double Eb = m_faces(i)->bendingEnergyContentDensity();
		// double Eg = m_faces(i)->connectionEnergyContentDensity();
		double Es = m_faces(i)->stretchingEnergy();
		double Eb = m_faces(i)->bendingEnergy();
		double Eg = m_faces(i)->connectionEnergy();
		a_face->write(i);
		a_face->tab();
		a_face->write(Es);
		a_face->tab();
		a_face->write(Eb);
		a_face->tab();
		a_face->write(Eg);            // ← write the γ‐energy
        a_face->tab();
		a_face->write(counter);
		a_face->newline();
	}

}
/* ============================================================================== */
/* Dump forms to file (text) */
// void NonEuclideanShell::DumpFormsTextFormat(TextFileHandle* a_face)
// {
// 	for (int i=0; i<m_faces.length(); i++)
// 	{
// 		TinyVector<double,3> EFGv = m_faces(i)->EFG();
// 		TinyVector<double,3> LMNv = m_faces(i)->LMN();
// 		TinyVector<double,2> UVv = m_faces(i)->coordinates();
// 		a_face->write(i);
// 		a_face->tab();
// 		a_face->write(UVv(0));
// 		a_face->tab();
// 		a_face->write(UVv(1));
// 		a_face->tab();
// 		a_face->write(EFGv(0));
// 		a_face->tab();
// 		a_face->write(EFGv(1));
// 		a_face->tab();
// 		a_face->write(EFGv(2));
// 		a_face->tab();
// 		a_face->write(LMNv(0));
// 		a_face->tab();
// 		a_face->write(LMNv(1));
// 		a_face->tab();
// 		a_face->write(LMNv(2));
// 		a_face->newline();
// 	}
// }
/* ============================================================================== */
/* Dump first, second, and third fundamental forms plus energy terms:
   (E,F,G), (L,M,N), (c11,c12,c22), Es, Eb, Eg */
void NonEuclideanShell::DumpFormsTextFormat(TextFileHandle* a_face)
{
    for (int i = 0; i < m_faces.length(); ++i)
    {
        auto   EFGv   = m_faces(i)->EFG();
        auto   LMNv   = m_faces(i)->LMN();
        auto   UVv    = m_faces(i)->coordinates();
		auto   ABCv   = m_faces(i)->ABC();
		auto   Cbarv  = m_faces(i)->cbar();
        // double Es  = m_faces(i)->stretchingEnergyContentDensity();
		// double Eb  = m_faces(i)->bendingEnergyContentDensity();
		// double Eg  = m_faces(i)->connectionEnergyContentDensity();

		double Es  = m_faces(i)->stretchingEnergy();
		double Eb  = m_faces(i)->bendingEnergy();
		double Eg  = m_faces(i)->connectionEnergy();

        a_face->write(i);        a_face->tab();
        a_face->write(UVv(0));   a_face->tab();
        a_face->write(UVv(1));   a_face->tab();

        a_face->write(EFGv(0));  a_face->tab();
        a_face->write(EFGv(1));  a_face->tab();
        a_face->write(EFGv(2));  a_face->tab();

        a_face->write(LMNv(0));  a_face->tab();
        a_face->write(LMNv(1));  a_face->tab();
        a_face->write(LMNv(2));  a_face->tab();

		a_face->write(ABCv(0));  a_face->tab();
		a_face->write(ABCv(1));  a_face->tab();
		a_face->write(ABCv(2));  a_face->tab();

		a_face->write(Cbarv(0,0));  a_face->tab();
		a_face->write(Cbarv(0,1));  a_face->tab();
		a_face->write(Cbarv(1,0));  a_face->tab();
		a_face->write(Cbarv(1,1));  a_face->tab();

        a_face->write(Es);  a_face->tab();
		a_face->write(Eb);  a_face->tab();
		a_face->write(Eg);  
        a_face->newline();
    }
}
/* ============================================================================== */
/* I/O of state and force. Needed for external optimization procedure */
void NonEuclideanShell::setPositionVector(const double *ptr)
{
	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			m_nodes(i)->position(0) = ptr[3*j];
			m_nodes(i)->position(1) = ptr[3*j+1];
			m_nodes(i)->position(2) = ptr[3*j+2];
			j += 1;
		}
	}

	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(i)->position() = m_nodes(j)->position() + m_nodes(i)->offset();
		}
	}
}

void NonEuclideanShell::setPositionVector(const gsl_vector *a_vec)
{
	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			m_nodes(i)->position(0) = gsl_vector_get(a_vec,3*j);
			m_nodes(i)->position(1) = gsl_vector_get(a_vec,3*j+1);
			m_nodes(i)->position(2) = gsl_vector_get(a_vec,3*j+2);
			j += 1;
		}
	}

	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(i)->position() = m_nodes(j)->position() + m_nodes(i)->offset();
		}
	}
}
/* ============================================================================== */
/* I/O of state and force. Needed for external optimization procedure */
void NonEuclideanShell::getPositionVector(gsl_vector *a_vec)
{
	int j = 0;
	double* ptr = gsl_vector_ptr(a_vec, 0);
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			ptr[3*j]   = m_nodes(i)->position(0);
			ptr[3*j+1] = m_nodes(i)->position(1);
			ptr[3*j+2] = m_nodes(i)->position(2);
			j += 1;
		}
	}
}
/* ============================================================================== */
/* return the energy given an array containing the position */
double NonEuclideanShell::getEnergy(const gsl_vector *a_vec)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergy()");

	const double* ptr = gsl_vector_const_ptr(a_vec, 0);
	setPositionVector(ptr);
	return energy();
}
/* ============================================================================== */
/* return the energy gradient given an array containing the position */
void NonEuclideanShell::getEnergyGradient(const gsl_vector *a_state, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyGradient()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	initializeForce();
	setForce();

	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(j)->force() += m_nodes(i)->force();
		}
	}

	j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			ret_ptr[3*j]   = m_nodes(i)->force(0);
			ret_ptr[3*j+1] = m_nodes(i)->force(1);
			ret_ptr[3*j+2] = m_nodes(i)->force(2);
			j += 1;
		}
	}
}
/* ============================================================================== */
/* return the energy and its gradient given an array containing the position */
void NonEuclideanShell::getEnergyAndEnergyGradient(const gsl_vector *a_state,
								  double *a_energy, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyAndEnergyGradient()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	initializeForce();
	setForce();

	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(j)->force() += m_nodes(i)->force();
		}
	}

	j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			ret_ptr[3*j]   = m_nodes(i)->force(0);
			ret_ptr[3*j+1] = m_nodes(i)->force(1);
			ret_ptr[3*j+2] = m_nodes(i)->force(2);
			j += 1;
		}
	}
	*a_energy = energy();
}
/* ============================================================================== */
/* return the energy gradient given an array containing the position, Full calculation */
void NonEuclideanShell::getEnergyGradientFull(const gsl_vector *a_state, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyGradientFull()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	double ep    = 1.e-6;
	double Einit = energy();
	int j = 0;

	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			for (int comp=0; comp<3; comp++)
			{
				m_nodes(i)->position(comp) += ep;
				double Eplus = energy();
				m_nodes(i)->position(comp) -= 2*ep;
				double Eminus = energy();
				ret_ptr[3*j+comp] = 0.5*(Eplus-Eminus)/ep;
				m_nodes(i)->position(comp) += ep;
			}
			j += 1;
		}
	}
}
/* ============================================================================== */
/* return the energy and its gradient given an array containing the position, Full calculation */
void NonEuclideanShell::getEnergyAndEnergyGradientFull(const gsl_vector *a_state,
								  double *a_energy, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyAndEnergyGradientFull()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	double ep    = 1.e-6;
	double Einit = energy();
	int j = 0;

	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			for (int comp=0; comp<3; comp++)
			{
				m_nodes(i)->position(comp) += ep;
				double Eplus = energy();
				m_nodes(i)->position(comp) -= 2*ep;
				double Eminus = energy();
				ret_ptr[3*j+comp] = 0.5*(Eplus-Eminus)/ep;
				m_nodes(i)->position(comp) += ep;
			}
			j += 1;
		}
	}
	*a_energy = energy();
}

/* ============================================================================== */
/* Test the gradient calculation                                                 */
void NonEuclideanShell::testGradient()
{
	int    noDOF = SizeOfOptimizationProblem();
	double grad[noDOF], gradFull[noDOF];

	/* calculate force */
	initializeForce();
	setForce();
	int j = 0;

	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			grad[3*j]   = m_nodes(i)->force(0);
			grad[3*j+1] = m_nodes(i)->force(1);
			grad[3*j+2] = m_nodes(i)->force(2);
			j += 1;
		}
	}

	j = 0;

	/* calculate force directly */
	double ep    = 1.e-6;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			for (int comp=0; comp<3; comp++)
			{
				m_nodes(i)->position(comp) += ep;
				double Eplus = energy();
				m_nodes(i)->position(comp) -= 2*ep;
				double Eminus = energy();
				gradFull[3*j+comp] = 0.5*(Eplus-Eminus)/ep;
				m_nodes(i)->position(comp) += ep;
			}
			j += 1;
		}
	}

	/* output */
	for (int i=0; i<noDOF; i++)
	{
		// std::cout << grad[i] << " " << gradFull[i] << " " <<  fabs(gradFull[i]-grad[i]) << std::endl;
		std::cout << grad[i] << " " << gradFull[i] << " " << fabs(gradFull[i]-grad[i])
          << " " << (fabs(gradFull[i]) > 1e-10 ? fabs(gradFull[i]-grad[i])/fabs(gradFull[i]) : 0.0)
          << std::endl;	

	}


	}
	const Vector<Face*>& NonEuclideanShell::getFaces() const {
		return m_faces;
	}
/* ============================================================================== */
void NonEuclideanShell::checkNodePositions() const {
    for (int i=0; i<m_nodes.length(); ++i) {
        for (int d=0; d<3; ++d) {
            if (!std::isfinite(m_nodes(i)->position(d))) {
                std::cerr << "Invalid node " << i << " pos: "
                          << m_nodes(i)->position(0) << ", "
                          << m_nodes(i)->position(1) << ", "
                          << m_nodes(i)->position(2) << std::endl;
                exit(1);
            }
        }
    }
}

