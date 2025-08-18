#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

class ReducedSystemOperator;

/** After spatial discretization, the hyperelastic model can be written as a
 *  system of ODEs:
 *     dv/dt = -M^{-1}*(H(x) + S*v)
 *     dx/dt = v,
 *  where x is the vector representing the deformation, v is the velocity field,
 *  M is the mass matrix, S is the viscosity matrix, and H(x) is the nonlinear
 *  hyperelastic operator.
 *
 *  Class DOHyperelasticOperator represents the right-hand side of the above
 *  system of ODEs. */

 class DOHyperelasticOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list;

   ParBilinearForm M, S;
   ParNonlinearForm H;
   real_t viscosity;
   HyperelasticModel *model;

   std::vector<real_t> &poles; 
   std::vector<real_t> &weights; 
   VectorFunctionCoefficient &rhs_cf; 
   VectorFunctionCoefficient &bd_cf; 

   HypreParMatrix *Mmat; // Mass matrix from ParallelAssemble()
   // CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   /** Nonlinear operator defining the reduced backward Euler equation for the
       velocity. Used in the implementation of method ImplicitSolve. */
   ReducedSystemOperator *reduced_oper;

   /// Newton solver for the reduced backward Euler equation
   NewtonSolver newton_solver;

   /// Solver for the Jacobian solve in the Newton method
   Solver *J_solver;
   /// Preconditioner for the Jacobian solve in the Newton method
   Solver *J_prec;

   mutable Vector z,k;  // auxiliary vector

public:

   CGSolver M_solver;  

   DOHyperelasticOperator(ParFiniteElementSpace &f, Array<int> &ess_bdr,
                        real_t visc, real_t mu, real_t K, 
                        std::vector<real_t> &weights, std::vector<real_t> &poles, 
                        VectorFunctionCoefficient &rhs_fc, 
                        VectorFunctionCoefficient &bd_fc);

   /// Compute the right-hand side of the ODE system.
   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   void ImplicitSolve(const real_t dt, const Vector &x, const BlockVector &modes, Vector &k);

   void IEStep(Vector &x, BlockVector &modes, real_t &t, real_t & dt); 

   // y = M^-1 H x  
   void TestH(Vector &x, Vector &y); 

   real_t ElasticEnergy(const ParGridFunction &x) const;
   real_t KineticEnergy(const ParGridFunction &v) const;
   void GetElasticEnergyDensity(const ParGridFunction &x,
                                ParGridFunction &w) const;

   ~DOHyperelasticOperator() override;
};

/** Nonlinear operator of the form:
    k --> (M + dt*S)*k + H(x + dt*v + dt^2*k) + S*v,
    where M and S are given BilinearForms, H is a given NonlinearForm, v and x
    are given vectors, and dt is a scalar. */
class ReducedSystemOperator : public Operator
{
private:
   ParBilinearForm *M, *S;
   ParNonlinearForm *H;
   mutable HypreParMatrix *Jacobian;
   real_t dt;
   const Vector *v, *x;
   const BlockVector *modes; 
   mutable Vector w, z;
   const Array<int> &ess_tdof_list;
   real_t cS; 

   std::vector<real_t> &poles; 
   std::vector<real_t> &weights; 

public:
   ReducedSystemOperator(ParBilinearForm *M_, ParBilinearForm *S_,
                         ParNonlinearForm *H_, const Array<int> &ess_tdof_list, 
                         std::vector<real_t> &weights, std::vector<real_t> &poles);

   /// Set current dt, v, x values - needed to compute action and Jacobian.
   void SetParameters(real_t dt_, const Vector *v_, const Vector *x_, const BlockVector *modes_);

   /// Compute y = H(x + dt (v + dt k)) + M k + S (v + dt k).
   void Mult(const Vector &k, Vector &y) const override;

   /// Compute J = M + dt S + dt^2 grad_H(x + dt (v + dt k)).
   Operator &GetGradient(const Vector &k) const override;

   ~ReducedSystemOperator() override;
};

/** Function representing the elastic energy density for the given hyperelastic
    model+deformation. Used in DOHyperelasticOperator::GetElasticEnergyDensity. */
class ElasticEnergyCoefficient : public Coefficient
{
private:
   HyperelasticModel     &model;
   const ParGridFunction &x;
   DenseMatrix            J;

public:
   ElasticEnergyCoefficient(HyperelasticModel &m, const ParGridFunction &x_)
      : model(m), x(x_) { }
   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override;
   ~ElasticEnergyCoefficient() override { }
};

class JCoefficient : public Coefficient 
{
   private:
      const ParGridFunction &x;
      DenseMatrix            J;
      int i,j; 
   
   public:
      JCoefficient(const ParGridFunction &x_, int i_, int j_)
         : x(x_), i(i_), j(j_) { }
      
      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override;
   
      ~JCoefficient() override { }
   };

class PCoefficient : public Coefficient 
{
   private:
      const ParGridFunction &x;
      HyperelasticModel &model; 
      DenseMatrix       J, P;
      int i,j; 
   
   public:
      PCoefficient(const ParGridFunction &x_, HyperelasticModel &model_, int i_, int j_)
         : x(x_), model(model_), i(i_), j(j_), J(3), P(3) { }
      
      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override;
   
      ~PCoefficient() override { }
   };