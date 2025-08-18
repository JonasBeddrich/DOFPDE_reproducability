/** After spatial discretization, the wave model can be written as:
 *
 *     d^2u/dt^2 = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass,
 *  and K is the stiffness matrix.
 *
 *  Class WaveOperator represents the right-hand side of the above ODE.
 */

#include <fstream>
#include <iostream>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

class DODiffusionWaveOperator {
 protected:
  ParFiniteElementSpace &fespace;
  Array<int> ess_tdof_list;  // this list remains empty for pure Neumann b.c.
  int height;

  std::vector<real_t> weights01;
  std::vector<real_t> poles01;
  std::vector<real_t> weights12;
  std::vector<real_t> poles12;

  FunctionCoefficient rhs_fc;

  ParBilinearForm *M;
  ParBilinearForm *K;

  HypreParMatrix Mmat, Kmat;
  HypreParMatrix *T;  // T = M + dt K
  real_t t;
  real_t current_dt;

  CGSolver M_solver;      // Krylov solver for inverting the mass matrix M
  HypreBoomerAMG M_prec;  // Preconditioner for the mass matrix M

  CGSolver T_solver;      // Implicit solver for T = M + fac0*K
  HypreBoomerAMG T_prec;  // Preconditioner for the implicit solver

  Coefficient *c2;
  mutable Vector zn;   // auxiliary vector
  mutable Vector tmp;  // auxiliary vector

 public:
  DODiffusionWaveOperator(ParFiniteElementSpace &f, Array<int> &ess_bdr,
                          real_t speed, std::vector<real_t> &w01,
                          std::vector<real_t> &p01, std::vector<real_t> &w12,
                          std::vector<real_t> &p12, FunctionCoefficient &rhs);

  void StepIE(Vector &u, Vector &du_dt, Vector &d2u_dt2, BlockVector &modes01,
              BlockVector &modes12, real_t &t, real_t &dt);

  // TODO
  // void StepRIIA2(Vector &u, Vector &du_dt, BlockVector &modes01, BlockVector
  // &modes12, real_t &t, real_t &dt);

  void ComputeRightHandSide(const Vector &u, const Vector &du_dt,
                            const BlockVector &modes01,
                            const BlockVector &modes12, const real_t &t,
                            const real_t &dt);

  void SolveImplicitProblem(Vector &d2u_dt2, const real_t &dt);

  void UpdateSolution(Vector &u, Vector &du_dt, Vector &d2u_dt2,
                      const real_t &dt);

  void UpdateModes(Vector &u, Vector &du_dt, Vector &d2u_dt2,
                   BlockVector &modes01, BlockVector &modes12, real_t &dt);

  void SetParameters(const Vector &u);

  ~DODiffusionWaveOperator();
};