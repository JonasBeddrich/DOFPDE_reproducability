#include "dodiffusionwaveoperator.hpp"

#include <fstream>
#include <iostream>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

DODiffusionWaveOperator::DODiffusionWaveOperator(
    ParFiniteElementSpace &f, Array<int> &ess_bdr, real_t speed,
    std::vector<real_t> &w01, std::vector<real_t> &p01,
    std::vector<real_t> &w12, std::vector<real_t> &p12,
    FunctionCoefficient &rhs)
    : fespace(f),
      height(fespace.GetTrueVSize()),
      weights01(w01),
      poles01(p01),
      weights12(w12),
      poles12(p12),
      rhs_fc(rhs),
      M(NULL),
      K(NULL),
      T(NULL),
      t(0.0),
      current_dt(0.0),
      M_solver(f.GetComm()),
      T_solver(f.GetComm()),
      zn(height),
      tmp(height) {
  const real_t rel_tol = 1e-12;
  c2 = new ConstantCoefficient(speed * speed);
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  // Assemble Mass matrix
  M = new ParBilinearForm(&fespace);
  M->AddDomainIntegrator(new MassIntegrator());
  M->Assemble(0);
  M->Finalize();
  M->FormSystemMatrix(ess_tdof_list, Mmat);

  // Assemble Laplace matrix
  K = new ParBilinearForm(&fespace);
  K->AddDomainIntegrator(new DiffusionIntegrator(*c2));
  K->Assemble(0);
  K->Finalize();
  K->FormSystemMatrix(ess_tdof_list, Kmat);

  // Configure preconditioner
  M_solver.iterative_mode = false;
  M_solver.SetRelTol(rel_tol);
  M_solver.SetAbsTol(0.0);
  M_solver.SetMaxIter(100);
  M_solver.SetPrintLevel(0);
  M_solver.SetPreconditioner(M_prec);
  M_solver.SetOperator(Mmat);

  // Configure solver
  T_solver.iterative_mode = false;
  T_solver.SetRelTol(rel_tol);
  T_solver.SetAbsTol(0.0);
  T_solver.SetMaxIter(100);
  T_solver.SetPrintLevel(0);
  // T_prec.SetPrintLevel(0);
  // T_solver.SetPreconditioner(T_prec);
}

void DODiffusionWaveOperator::StepIE(Vector &u, Vector &du_dt, Vector &d2u_dt2,
                                     BlockVector &modes01, BlockVector &modes12,
                                     real_t &t, real_t &dt) {
  // Compute all the stuff that is left
  ComputeRightHandSide(u, du_dt, modes01, modes12, t, dt);
  // Solve d2u_dt2
  SolveImplicitProblem(d2u_dt2, dt);
  // Compute u and dudt
  UpdateSolution(u, du_dt, d2u_dt2, dt);
  // Update the fractional modes
  UpdateModes(u, du_dt, d2u_dt2, modes01, modes12, dt);
  // Update time
  t += dt;
}

void DODiffusionWaveOperator::ComputeRightHandSide(
    const Vector &u, const Vector &du_dt, const BlockVector &modes01,
    const BlockVector &modes12, const real_t &t, const real_t &dt) {
  ParGridFunction rhs_gf(&fespace);
  rhs_fc.SetTime(t + dt);
  rhs_gf.ProjectCoefficient(rhs_fc);
  Vector rhs;
  rhs_gf.GetTrueDofs(rhs);

  // TODO check time t + dt for rhs

  tmp = rhs;

  real_t c{0.0};
  for (size_t i = 0; i < weights01.size(); ++i) {
    c += weights01[i] * dt / (1. + poles01[i] * dt);
  }
  tmp.Add(c, du_dt);

  for (size_t i = 0; i < poles01.size(); ++i) {
    tmp.Add(1. / (1. + poles01[i] * dt), modes01.GetBlock(i));
  }

  for (size_t i = 0; i < poles12.size(); ++i) {
    tmp.Add(1. / (1. + poles12[i] * dt), modes12.GetBlock(i));
  }

  tmp.Neg();
  Mmat.Mult(tmp, zn);

  // Sign ?
  Kmat.AddMult(u, zn, -1.0);
  Kmat.AddMult(du_dt, zn, -dt);
  zn.SetSubVector(ess_tdof_list, 0.0);
}

void DODiffusionWaveOperator::SolveImplicitProblem(Vector &d2u_dt2,
                                                   const real_t &dt) {
  if (!T) {
    real_t c{0.0};
    for (size_t i = 0; i < weights01.size(); ++i) {
      c += weights01[i] * dt * dt / (1. + poles01[i] * dt);
    }

    for (size_t i = 0; i < weights12.size(); ++i) {
      c += weights12[i] * dt / (1. + poles12[i] * dt);
    }
    T = Add(c, Mmat, dt * dt, Kmat);
    T_solver.SetOperator(*T);
  }
  T_solver.Mult(zn, d2u_dt2);
  d2u_dt2.SetSubVector(ess_tdof_list, 0.0);
}

void DODiffusionWaveOperator::UpdateSolution(Vector &u, Vector &du_dt,
                                             Vector &d2u_dt2,
                                             const real_t &dt) {
  du_dt.Add(dt, d2u_dt2);
  u.Add(dt, du_dt);
}

void DODiffusionWaveOperator::UpdateModes(Vector &u, Vector &du_dt,
                                          Vector &d2u_dt2, BlockVector &modes01,
                                          BlockVector &modes12, real_t &dt) {
  // Update the 01 modes
  for (size_t i = 0; i < weights01.size(); ++i) {
    modes01.GetBlock(i).Add(weights01[i] * dt, du_dt);
    modes01.GetBlock(i) /= (1 + poles01[i] * dt);
  }
  // Update the 12 modes
  for (size_t i = 0; i < weights12.size(); ++i) {
    modes12.GetBlock(i).Add(weights12[i] * dt, d2u_dt2);
    modes12.GetBlock(i) /= (1 + poles12[i] * dt);
  }
}

DODiffusionWaveOperator::~DODiffusionWaveOperator() {
  delete T;
  delete M;
  delete K;
  delete c2;
}