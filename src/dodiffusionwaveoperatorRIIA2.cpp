#include "dodiffusionwaveoperatorRIIA2.hpp"

#include <fstream>
#include <iostream>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

DODiffusionWaveOperatorRIIA2::DODiffusionWaveOperatorRIIA2(
    ParFiniteElementSpace &f, Array<int> &ess_bdr, real_t speed,
    std::vector<real_t> &w01, std::vector<real_t> &p01,
    std::vector<real_t> &w12, std::vector<real_t> &p12,
    FunctionCoefficient &rhs)
    : fespace(f),
      height(fespace.GetTrueVSize()),
      s(2),
      first(true),
      blockTrueOffsets(s + 1),
      weights01(w01),
      poles01(p01),
      weights12(w12),
      poles12(p12),
      AB1jone(weights01.size()),
      AB2jone(weights12.size()),
      B1j(weights01.size()),
      B2j(weights12.size()),
      AB1jA(weights01.size()),
      AB2j(weights12.size()),
      rhs_fc(rhs),
      M(NULL),
      K(NULL),
      T(NULL),
      t(0.0),
      current_dt(0.0),
      M_solver(f.GetComm()),
      T_solver(f.GetComm()),
      tmp(height) {
  c2 = new ConstantCoefficient(speed * speed);
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  coefficients.SetSize(s);
  A.SetSize(s);
  b.SetSize(s);
  c.SetSize(s);

  if (s == 1) {
    A(0, 0) = 1.;
    b(0) = 1.;
    c(0) = 1.;

  } else if (s == 2) {
    // real_t alpha = 1. / 3.;

    // A(0, 0) = alpha;
    // A(0, 1) = 0.;
    // A(1, 0) = 3. / 4.;
    // A(1, 1) = 1. / 4.;

    // b(0) = A(1, 0);
    // b(1) = A(1, 1);

    // c(0) = alpha;
    // c(1) = 1;

    A(0, 0) = 5. / 12.;
    A(0, 1) = -1. / 12.;
    A(1, 0) = 3. / 4.;
    A(1, 1) = 1. / 4.;

    b(0) = 3. / 4.;
    b(1) = 1. / 4.;

    c(0) = 1. / 3.;
    c(1) = 1;
  }

  A2.SetSize(s);
  Mult(A, A, A2);

  one.SetSize(s);
  one = 1.;

  for (int j = 0; j < weights01.size(); ++j) {
    AB1jone[j].SetSize(s);
    B1j[j].SetSize(s);
    AB1jA[j].SetSize(s);

    AB1jone[j] = 0.;
    B1j[j] = 0.;
    AB1jA[j] = 0.;
  }

  for (int j = 0; j < weights12.size(); ++j) {
    AB2jone[j].SetSize(s);
    B2j[j].SetSize(s);
    AB2j[j].SetSize(s);

    AB2jone[j] = 0.;
    B2j[j] = 0.;
    AB2j[j] = 0.;
  }

  blockTrueOffsets = fespace.TrueVSize();
  blockTrueOffsets.PartialSum();

  tmp_bv.Update(blockTrueOffsets);
  zn_bv.Update(blockTrueOffsets);
  k1_bv.Update(blockTrueOffsets);
  k0_bv.Update(blockTrueOffsets);

  tmp_bv = 0.;
  zn_bv = 0.;
  k1_bv = 0.;
  k0_bv = 0.;

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

  // Configure solver
  T_solver.iterative_mode = true;
  T_solver.SetRelTol(1e-12);
  T_solver.SetMaxIter(1000);
  T_solver.SetPrintLevel(0);

  // Configure solver
  T11_solver.iterative_mode = true;
  T11_solver.SetAbsTol(1e-12);
  T11_solver.SetMaxIter(1000);
  T11_solver.SetPrintLevel(0);

  // Configure solver
  T22_solver.iterative_mode = true;
  T22_solver.SetAbsTol(1e-12);
  T22_solver.SetMaxIter(1000);
  T22_solver.SetPrintLevel(0);
}

void DODiffusionWaveOperatorRIIA2::Step(Vector &u, Vector &du_dt,
                                        Vector &d2u_dt2, BlockVector &modes01,
                                        BlockVector &modes12, real_t &t,
                                        real_t &dt) {
  ComputeCoefficients(dt);
  // Compute all the stuff that is left
  ComputeRightHandSide(u, du_dt, d2u_dt2, modes01, modes12, t, dt);
  // Solve d2u_dt2
  SolveImplicitProblem(dt);
  // Compute u and dudt
  UpdateSolution(u, du_dt, d2u_dt2, dt);
  // Update the fractional modes
  UpdateModes(u, du_dt, d2u_dt2, modes01, modes12, dt);
  // Update time
  t += dt;
}

void DODiffusionWaveOperatorRIIA2::ComputeCoefficients(const real_t &dt) {
  Vector tmp_v(s);
  DenseMatrix tmp_m(s);
  tmp_v = 0.;
  tmp_m = 0.;

  for (size_t j = 0; j < poles01.size(); ++j) {
    B1j[j].Diag(1., s);
    B1j[j].Add(dt * poles01[j], A);
    B1j[j].Invert();
    B1j[j].Mult(one, tmp_v);

    A.Mult(tmp_v, AB1jone[j]);

    Mult(A, B1j[j], tmp_m);
    Mult(tmp_m, A, AB1jA[j]);
  }

  for (size_t j = 0; j < poles12.size(); ++j) {
    B2j[j].Diag(1., s);
    B2j[j].Add(dt * poles12[j], A);
    B2j[j].Invert();

    B2j[j].Mult(one, tmp_v);
    A.Mult(tmp_v, AB2jone[j]);

    Mult(A, B2j[j], AB2j[j]);
  }
}

void DODiffusionWaveOperatorRIIA2::ComputeRightHandSide(
    const Vector &u, const Vector &du_dt, const Vector &d2u_dt2,
    const BlockVector &modes01, const BlockVector &modes12, const real_t &t,
    const real_t &dt) {
  tmp_bv = 0.;

  for (size_t i = 0; i < s; ++i) {
    // sum h w_1j v^(1) AB_1j 1
    for (size_t j = 0; j < weights01.size(); ++j) {
      tmp_bv.GetBlock(i).Add(dt * weights01[j] * AB1jone[j](i), du_dt);
    }

    // - sum sum h lambda_ij v_ij^n A B_ij 1Mmat
    for (size_t j = 0; j < weights01.size(); ++j) {
      tmp_bv.GetBlock(i).Add(-dt * poles01[j] * AB1jone[j](i),
                             modes01.GetBlock(j));
    }
    for (size_t j = 0; j < weights12.size(); ++j) {
      tmp_bv.GetBlock(i).Add(-dt * poles12[j] * AB2jone[j](i),
                             modes12.GetBlock(j));
    }

    // sum sum v_ij^n 1
    for (size_t j = 0; j < weights01.size(); ++j) {
      tmp_bv.GetBlock(i).Add(1.0, modes01.GetBlock(j));
    }
    for (size_t j = 0; j < weights12.size(); ++j) {
      tmp_bv.GetBlock(i).Add(1.0, modes12.GetBlock(j));
    }

    tmp_bv.Neg();

    // rhs
    ParGridFunction rhs_gf(&fespace);
    rhs_fc.SetTime(t + c(i) * dt);
    rhs_gf.ProjectCoefficient(rhs_fc);
    Vector rhs;
    rhs_gf.GetTrueDofs(rhs);
    tmp_bv.GetBlock(i).Add(-1.0, rhs);

    // multiply with mass matrix
    Mmat.Mult(tmp_bv.GetBlock(i), zn_bv.GetBlock(i));
  }

  for (size_t i = 0; i < s; ++i) {
    Kmat.AddMult(u, zn_bv.GetBlock(i), -1.0);
    for (size_t l = 0; l < s; ++l) {
      Kmat.AddMult(du_dt, zn_bv.GetBlock(i), -dt * A(i, l));
    }
  }
}

void DODiffusionWaveOperatorRIIA2::SolveImplicitProblem(const real_t &dt) {
  T = new BlockOperator(blockTrueOffsets);

  for (size_t i = 0; i < s; ++i) {
    for (size_t l = 0; l < s; ++l) {
      real_t tmp_f{0.0};
      for (size_t j = 0; j < weights01.size(); ++j) {
        tmp_f += weights01[j] * dt * dt * AB1jA[j](i, l);
      }
      for (size_t j = 0; j < weights12.size(); ++j) {
        tmp_f += weights12[j] * dt * AB2j[j](i, l);
      }
      coefficients(i, l) = tmp_f;
    }
  }

  for (size_t i = 0; i < s; ++i) {
    for (size_t l = 0; l < s; ++l) {
      T->SetBlock(i, l,
                  Add(coefficients(i, l), Mmat, dt * dt * A2(i, l), Kmat));
    }
  }
  T_solver.SetOperator(*T);

  T_solver.Mult(zn_bv, k1_bv);

  for (size_t i = 0; i < s; ++i) {
    k1_bv.GetBlock(i).SetSubVector(ess_tdof_list, 0.0);
  }
  delete T;
}

void DODiffusionWaveOperatorRIIA2::UpdateSolution(Vector &u, Vector &du_dt,
                                                  Vector &d2u_dt2,
                                                  const real_t &dt) {
  for (size_t i = 0; i < s; ++i) {
    k0_bv.GetBlock(i).Set(1.0, du_dt);
    for (size_t l = 0; l < s; ++l) {
      k0_bv.GetBlock(i).Add(dt * A(i, l), k1_bv.GetBlock(l));
    }
  }

  for (size_t i = 0; i < s; ++i) {
    du_dt.Add(dt * b(i), k1_bv.GetBlock(i));
    u.Add(dt * b(i), k0_bv.GetBlock(i));
  }
}

void DODiffusionWaveOperatorRIIA2::UpdateModes(Vector &u, Vector &du_dt,
                                               Vector &d2u_dt2,
                                               BlockVector &modes01,
                                               BlockVector &modes12,
                                               real_t &dt) {
  for (size_t j = 0; j < weights01.size(); ++j) {
    tmp_bv = 0.;
    for (size_t i = 0; i < s; ++i) {
      for (size_t l = 0; l < s; ++l) {
        tmp_bv.GetBlock(i).Add(-poles01[j] * B1j[j](i, l), modes01.GetBlock(j));
        tmp_bv.GetBlock(i).Add(weights01[j] * B1j[j](i, l), k0_bv.GetBlock(l));
      }
    }
    for (size_t i = 0; i < s; ++i) {
      modes01.GetBlock(j).Add(dt * b(i), tmp_bv.GetBlock(i));
    }
  }

  for (size_t j = 0; j < weights12.size(); ++j) {
    tmp_bv = 0.;
    for (size_t i = 0; i < s; ++i) {
      for (size_t l = 0; l < s; ++l) {
        tmp_bv.GetBlock(i).Add(-poles12[j] * B2j[j](i, l), modes12.GetBlock(j));
        tmp_bv.GetBlock(i).Add(weights12[j] * B2j[j](i, l), k1_bv.GetBlock(l));
      }
    }

    for (size_t i = 0; i < s; ++i) {
      modes12.GetBlock(j).Add(dt * b(i), tmp_bv.GetBlock(i));
    }
  }
}

DODiffusionWaveOperatorRIIA2::~DODiffusionWaveOperatorRIIA2() {
  // delete T;
  delete M;
  delete K;
  delete c2;
}