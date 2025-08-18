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

class DODiffusionWaveOperatorRIIA2 {
 protected:
  ParFiniteElementSpace &fespace;
  Array<int> ess_tdof_list;  // this list remains empty for pure Neumann b.c.
  int height;

  int s;
  Array<int> blockTrueOffsets;
  DenseMatrix A, Ainv, A2, AB1Ainv, AB2Ainv, AB2AinvAinv, coefficients;
  Vector AB1Ainvone, AB2Ainvone, AB2AinvAinvone;
  Vector b;
  Vector c;
  Vector one;

  bool first;

  std::vector<real_t> weights01;
  std::vector<real_t> poles01;
  std::vector<real_t> weights12;
  std::vector<real_t> poles12;

  std::vector<Vector> AB1jone;
  std::vector<Vector> AB2jone;

  std::vector<Vector> AB1jAone;
  std::vector<Vector> AB2jAone;

  std::vector<DenseMatrix> B1j;
  std::vector<DenseMatrix> B2j;

  std::vector<DenseMatrix> AB1j;
  std::vector<DenseMatrix> AB2j;

  std::vector<DenseMatrix> AB1jA;
  std::vector<DenseMatrix> AB2jA;

  std::vector<DenseMatrix> AB1jAA;

  FunctionCoefficient rhs_fc;

  ParBilinearForm *M;
  ParBilinearForm *K;

  HypreParMatrix Mmat, Kmat;
  BlockOperator *T;  // T = M + dt K
  real_t t;
  real_t current_dt;

  CGSolver M_solver;      // Krylov solver for inverting the mass matrix M
  HypreBoomerAMG M_prec;  // Preconditioner for the mass matrix M

  std::vector<HypreBoomerAMG *> Block_prec;
  FGMRESSolver T_solver;  // Implicit solver for T = M + fac0*K
  BlockDiagonalPreconditioner *T_BLT_prec;  // Preconditioner for the implicit

  GMRESSolver T11_solver;
  GMRESSolver T22_solver;
  GMRESSolver T33_solver;

  HypreBoomerAMG *A11_prec;
  HypreBoomerAMG *A22_prec;
  HypreBoomerAMG *A33_prec;

  Operator *A21_block;
  // solver

  Coefficient *c2;
  // mutable Vector zn;            // auxiliary vector
  mutable Vector tmp;
  mutable BlockVector tmp_bv, zn_bv, unc_bv, vnc_bv, rnc_bv, k0_bv, k1_bv,
      k2_bv;  // auxiliary vector

 public:
  DODiffusionWaveOperatorRIIA2(ParFiniteElementSpace &f, Array<int> &ess_bdr,
                               real_t speed, std::vector<real_t> &w01,
                               std::vector<real_t> &p01,
                               std::vector<real_t> &w12,
                               std::vector<real_t> &p12,
                               FunctionCoefficient &rhs);

  void Step(Vector &u, Vector &du_dt, Vector &d2u_dt2, BlockVector &modes01,
            BlockVector &modes12, real_t &t, real_t &dt);

  void ComputeCoefficients(const real_t &dt);
  void ComputeRightHandSide(const Vector &u, const Vector &du_dt,
                            const Vector &d2u_dt2, const BlockVector &modes01,
                            const BlockVector &modes12, const real_t &t,
                            const real_t &dt);

  void SolveImplicitProblem(const real_t &dt);
  void ComputeVnc(const Vector &u, const real_t &dt);
  void ComputeRnc(const Vector &du_dt, const real_t &dt);
  void UpdateSolution(Vector &u, Vector &du_dt, Vector &d2u_dt2,
                      const real_t &dt);

  void UpdateModes(Vector &u, Vector &du_dt, Vector &d2u_dt2,
                   BlockVector &modes01, BlockVector &modes12, real_t &dt);

  void SetParameters(const Vector &u);

  ~DODiffusionWaveOperatorRIIA2();
};