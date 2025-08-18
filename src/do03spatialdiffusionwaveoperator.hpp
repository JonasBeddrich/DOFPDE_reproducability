/** After spatial discretization, the wave model can be written as:
 *
 *     d^2u/dt^2 = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass,
 *  and K is the stiffness matrix.
 *
 *  Class WaveOperator represents the right-hand side of the above ODE.
 */


#include "mfem.hpp"

#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


class DO03SpatialDiffusionWaveOperator
{
protected:
   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.
   int height; 

   std::vector<Coefficient*> weights01;  
   std::vector<Coefficient*> poles01; 
   std::vector<Coefficient*> weights12;  
   std::vector<Coefficient*> poles12; 
   std::vector<Coefficient*> weights23;  
   std::vector<Coefficient*> poles23; 

   FunctionCoefficient rhs_fc; 

   ParBilinearForm *M;
   ParBilinearForm *K;

   HypreParMatrix Mmat, Kmat, SIPmat; 
   HypreParMatrix *T; // T = M + dt K
   std::vector<HypreParMatrix> GAMMAmat01, GAMMAmat12, GAMMAmat23; // 1 / (1 + lambda * h)
   std::vector<HypreParMatrix> BETAmat01, BETAmat12, BETAmat23; // wh / (1 + lambda * h)

   real_t t; 
   real_t dt;

   CGSolver M_solver; // Krylov solver for inverting the mass matrix M
   HypreBoomerAMG M_prec;  // Preconditioner for the mass matrix M

   CGSolver T_solver; // Implicit solver for T = M + fac0*K
   HypreBoomerAMG T_prec;  // Preconditioner for the implicit solver

   Coefficient *c2;
   mutable Vector zn; // auxiliary vector
   mutable Vector tmp; // auxiliary vector

public:
   DO03SpatialDiffusionWaveOperator(ParFiniteElementSpace &f, Array<int> &ess_bdr, real_t speed, 
      std::vector<Coefficient*> &w01, std::vector<Coefficient*> &p01, 
      std::vector<Coefficient*> &w12, std::vector<Coefficient*> &p12, 
      std::vector<Coefficient*> &w23, std::vector<Coefficient*> &p23, 
      FunctionCoefficient &rhs, real_t dt_
   );

   void StepIE(Vector &u, Vector &du_dt, Vector &d2u_dt2, Vector &d3u_dt3, 
      BlockVector &modes01, BlockVector &modes12, BlockVector &modes23, real_t &t, real_t &dt); 

   void ComputeRightHandSide(const Vector &u, const Vector &du_dt, const Vector &d2u_dt2, 
      const BlockVector &modes01, const BlockVector &modes12, const BlockVector &modes23, const real_t &t, const real_t &dt); 

   void SolveImplicitProblem(Vector &d3u_dt3, const real_t &dt); 

   void UpdateSolution(Vector &u, Vector &du_dt, Vector &d2u_dt2, const Vector &d3u_dt3, const real_t &dt); 
   
   void UpdateModes(const Vector &u, const Vector &du_dt, const Vector &d2u_dt2, const Vector &d3u_dt3, 
      BlockVector &modes01, BlockVector &modes12, BlockVector &modes23, real_t &dt); 

   void ComputeCoefficientMatrices();

   ~DO03SpatialDiffusionWaveOperator();
};