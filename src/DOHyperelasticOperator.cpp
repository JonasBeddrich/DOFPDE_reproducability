#include "DOHyperelasticOperator.hpp"

#include <fstream>
#include <iostream>
#include <memory>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

ReducedSystemOperator::ReducedSystemOperator(ParBilinearForm *M_,
                                             ParBilinearForm *S_,
                                             ParNonlinearForm *H_,
                                             const Array<int> &ess_tdof_list_,
                                             std::vector<real_t> &weights_,
                                             std::vector<real_t> &poles_)
    : Operator(M_->ParFESpace()->TrueVSize()),
      M(M_),
      S(S_),
      H(H_),
      Jacobian(NULL),
      dt(0.0),
      v(NULL),
      x(NULL),
      modes(NULL),
      w(height),
      z(height),
      ess_tdof_list(ess_tdof_list_),
      cS(0.0),
      weights(weights_),
      poles(poles_) {}

void ReducedSystemOperator::SetParameters(real_t dt_, const Vector *v_,
                                          const Vector *x_,
                                          const BlockVector *modes_) {
  dt = dt_;
  v = v_;
  x = x_;
  modes = modes_;
  cS = 0.;
  for (size_t i = 0; i < weights.size(); ++i) {
    cS += weights[i] * dt / (1. + poles[i] * dt);
  }
}

void ReducedSystemOperator::Mult(const Vector &k, Vector &y) const {
  add(*v, dt, k, w);  // w = v^n+1
  add(*x, dt, w, z);  // z = x^n+1

  H->Mult(z, y);
  M->TrueAddMult(k, y);
  z.Set(cS, w);
  for (size_t i = 0; i < weights.size(); ++i) {
    z.Add(1. / (1. + poles[i] * dt), modes->GetBlock(i));
  }
  S->TrueAddMult(z, y);
  y.SetSubVector(ess_tdof_list, 0.0);
}

Operator &ReducedSystemOperator::GetGradient(const Vector &k) const {
  delete Jacobian;
  SparseMatrix *localJ = Add(1.0, M->SpMat(), dt * cS, S->SpMat());
  add(*v, dt, k, w);
  add(*x, dt, w, z);
  localJ->Add(dt * dt, H->GetLocalGradient(z));
  Jacobian = M->ParallelAssemble(localJ);
  delete localJ;
  HypreParMatrix *Je = Jacobian->EliminateRowsCols(ess_tdof_list);
  delete Je;
  return *Jacobian;
}

ReducedSystemOperator::~ReducedSystemOperator() { delete Jacobian; }

DOHyperelasticOperator::DOHyperelasticOperator(
    ParFiniteElementSpace &f, Array<int> &ess_bdr, real_t visc, real_t mu,
    real_t K, std::vector<real_t> &weights_, std::vector<real_t> &poles_,
    VectorFunctionCoefficient &rhs_cf_, VectorFunctionCoefficient &bd_cf_)
    : TimeDependentOperator(2 * f.TrueVSize(), (real_t)0.0),
      fespace(f),
      M(&fespace),
      S(&fespace),
      H(&fespace),
      viscosity(visc),
      M_solver(f.GetComm()),
      newton_solver(f.GetComm()),
      z(height / 2),
      k(height),
      weights(weights_),
      poles(poles_),
      rhs_cf(rhs_cf_),
      bd_cf(bd_cf_) {
  const real_t rel_tol = 1e-8;
  const real_t newton_abs_tol = 1e-8;

  const int skip_zero_entries = 0;

  const real_t ref_density = 1.0;  // density in the reference configuration
  ConstantCoefficient rho0(ref_density);
  M.AddDomainIntegrator(new VectorMassIntegrator(rho0));
  M.Assemble(skip_zero_entries);
  M.Finalize(skip_zero_entries);
  Mmat = M.ParallelAssemble();
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  HypreParMatrix *Me = Mmat->EliminateRowsCols(ess_tdof_list);
  delete Me;

  M_solver.iterative_mode = false;
  M_solver.SetRelTol(rel_tol);
  M_solver.SetAbsTol(0.0);
  M_solver.SetMaxIter(30);
  M_solver.SetPrintLevel(0);
  M_prec.SetType(HypreSmoother::Jacobi);
  M_solver.SetPreconditioner(M_prec);
  M_solver.SetOperator(*Mmat);

  model = new NeoHookeanModel(mu, K);
  H.AddDomainIntegrator(new HyperelasticNLFIntegrator(model));
  H.SetEssentialTrueDofs(ess_tdof_list);

  ConstantCoefficient visc_coeff(viscosity);
  S.AddDomainIntegrator(new VectorDiffusionIntegrator(visc_coeff));
  S.Assemble(skip_zero_entries);
  S.Finalize(skip_zero_entries);

  reduced_oper =
      new ReducedSystemOperator(&M, &S, &H, ess_tdof_list, weights, poles);

  HypreSmoother *J_hypreSmoother = new HypreSmoother;
  J_hypreSmoother->SetType(HypreSmoother::l1Jacobi);
  J_hypreSmoother->SetPositiveDiagonal(true);
  J_prec = J_hypreSmoother;

  MINRESSolver *J_minres = new MINRESSolver(f.GetComm());
  J_minres->SetRelTol(rel_tol);
  J_minres->SetAbsTol(0.0);
  J_minres->SetMaxIter(300);
  J_minres->SetPrintLevel(-1);
  J_minres->SetPreconditioner(*J_prec);
  J_solver = J_minres;

  newton_solver.iterative_mode = false;
  newton_solver.SetSolver(*J_solver);
  newton_solver.SetOperator(*reduced_oper);
  newton_solver.SetPrintLevel(1);  // print Newton iterations
  newton_solver.SetRelTol(rel_tol);
  newton_solver.SetAbsTol(newton_abs_tol);
  newton_solver.SetAdaptiveLinRtol(2, 0.5, 0.9);
  newton_solver.SetMaxIter(10);
}

void DOHyperelasticOperator::IEStep(Vector &vx, BlockVector &modes, real_t &t,
                                    real_t &dt) {
  rhs_cf.SetTime(t + dt);
  bd_cf.SetTime(t + dt);
  ImplicitSolve(dt, vx, modes, k);  // solve for k: k = f(x + dt*k, t + dt)
  // update solution
  vx.Add(dt, k);
  // update modes
  int sc = height / 2;
  Vector v(vx.GetData() + 0, sc);
  for (size_t i = 0; i < weights.size(); ++i) {
    modes.GetBlock(i).Add(dt * weights[i], v);
    modes.GetBlock(i) *= 1. / (1. + dt * poles[i]);
  }
}

void DOHyperelasticOperator::TestH(Vector &vx, Vector &y) {
  int sc = height / 2;
  Vector v(vx.GetData() + 0, sc);
  Vector x(vx.GetData() + sc, sc);

  H.Mult(x, k);
  k *= -1;
  M_solver.Mult(k, y);
}

void DOHyperelasticOperator::ImplicitSolve(const real_t dt, const Vector &vx,
                                           const BlockVector &modes,
                                           Vector &dvx_dt) {
  // By eliminating kx from the coupled system:
  //    kv = -M^{-1}*[H(x + dt*kx) + S*(v + dt*kv)]
  //    kx = v + dt*kv
  // we reduce it to a nonlinear equation for kv, represented by the
  // reduced_oper. This equation is solved with the newton_solver
  // object (using J_solver and J_prec internally).

  int sc = height / 2;
  Vector v(vx.GetData() + 0, sc);
  Vector x(vx.GetData() + sc, sc);
  Vector dv_dt(dvx_dt.GetData() + 0, sc);
  Vector dx_dt(dvx_dt.GetData() + sc, sc);

  ParLinearForm rhs(&fespace);
  rhs.AddDomainIntegrator(new VectorDomainLFIntegrator(rhs_cf));
  rhs.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(bd_cf));
  rhs.Assemble();

  Vector rhs_vec(x.Size());
  rhs.ParallelAssemble(rhs_vec);

  //  HypreParMatrix tmp_sp;
  //  ParGridFunction tmp_gf(&fespace);
  //  tmp_gf = 0.;
  //  Vector TMP, RHS;
  //  M.FormLinearSystem(ess_tdof_list, tmp_gf, rhs, tmp_sp, TMP, RHS);

  reduced_oper->SetParameters(dt, &v, &x, &modes);
  newton_solver.Mult(rhs_vec, dv_dt);
  MFEM_VERIFY(newton_solver.GetConverged(), "Newton solver did not converge.");
  add(v, dt, dv_dt, dx_dt);
}

real_t DOHyperelasticOperator::ElasticEnergy(const ParGridFunction &x) const {
  return H.GetEnergy(x);
}

real_t DOHyperelasticOperator::KineticEnergy(const ParGridFunction &v) const {
  real_t energy = 0.5 * M.ParInnerProduct(v, v);
  return energy;
}

void DOHyperelasticOperator::GetElasticEnergyDensity(const ParGridFunction &x,
                                                     ParGridFunction &w) const {
  ElasticEnergyCoefficient w_coeff(*model, x);
  w.ProjectCoefficient(w_coeff);
}

DOHyperelasticOperator::~DOHyperelasticOperator() {
  delete J_solver;
  delete J_prec;
  delete reduced_oper;
  delete model;
  delete Mmat;
}

real_t ElasticEnergyCoefficient::Eval(ElementTransformation &T,
                                      const IntegrationPoint &ip) {
  model.SetTransformation(T);
  x.GetVectorGradient(T, J);
  // return model.EvalW(J);  // in reference configuration
  return model.EvalW(J) / J.Det();  // in deformed configuration
}

real_t JCoefficient::Eval(ElementTransformation &T,
                          const IntegrationPoint &ip) {
  x.GetVectorGradient(T, J);
  return J(i, j);
}

real_t PCoefficient::Eval(ElementTransformation &T,
                          const IntegrationPoint &ip) {
  model.SetTransformation(T);
  x.GetVectorGradient(T, J);
  model.EvalP(J, P);
  return P(i, j);
}