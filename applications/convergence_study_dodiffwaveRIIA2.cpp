#include <fstream>
#include <iostream>

#include "dodiffusionwaveoperatorRIIA2.hpp"
#include "mfem.hpp"

using namespace std;
using namespace mfem;

real_t InitialSolution(const Vector &x) { return 0.0; }

real_t InitialRate(const Vector &x) { return 0.0; }

real_t rhs_function(const Vector &x, real_t t) {
  real_t sinsin = sin(4 * M_PI * x[0]) * sin(4 * M_PI * x[1]);
  real_t t5 = t * t * t * t * t;
  real_t t3 = t * t * t;
  real_t D_phi = 120 * (t5 - t3 * exp(-2)) / (1 + log(t));
  // return -sinsin * D_phi;
  return -sinsin * (D_phi + 32 * M_PI * M_PI * t5);
}

real_t solution(const Vector &x, real_t t) {
  return sin(4 * M_PI * x[0]) * sin(4 * M_PI * x[1]) * t * t * t * t * t;
}

real_t rate(const Vector &x, real_t t) {
  return sin(4 * M_PI * x[0]) * sin(4 * M_PI * x[1]) * 5. * t * t * t * t;
}

void grad_solution(const Vector &x, real_t t, Vector &y) {
  // y[0] = 4 * M_PI * cos(4 * M_PI * x[0]) * sin(4 * M_PI * x[1]) * t * t * t *
  //        t * t;
  // y[1] = 4 * M_PI * sin(4 * M_PI * x[0]) * cos(4 * M_PI * x[1]) * t * t * t *
  //        t * t;
}

void ReadWeightsPoles(std::vector<real_t> &weights, std::vector<real_t> &poles,
                      const char *filename) {
  std::ifstream file(filename);
  std::string line, item;

  if (std::getline(file, line)) {
    std::stringstream ss(line);
    while (std::getline(ss, item, ',')) {
      weights.push_back(std::stod(item));
    }
  }

  if (std::getline(file, line)) {
    std::stringstream ss(line);
    while (std::getline(ss, item, ',')) {
      poles.push_back(-1.0 * std::stod(item));
    }
  }
}

void runSimulation(real_t dt, int sr, int order);

int main(int argc, char *argv[]) {
  Mpi::Init(argc, argv);
  int num_procs = Mpi::WorldSize();
  int myid = Mpi::WorldRank();
  Hypre::Init();

  int order = 2;

  int spatial_refinement_levels = 5;
  int temporal_refinement_levels = 15;
  int refinement_levels = 10;
  real_t spatial_refinement_study_dt = std::pow(2, -20);
  real_t temporal_refinement_study_sr = 5;

  OptionsParser args(argc, argv);
  args.AddOption(&order, "-o", "--order", "");
  args.AddOption(&spatial_refinement_levels, "-srl",
                 "--spatial_refinement_levels", "");
  args.AddOption(&temporal_refinement_levels, "-trl",
                 "--temporal_refinement_levels", "");
  args.AddOption(&refinement_levels, "-rl", "--combined_refinement_levels", "");
  args.AddOption(&spatial_refinement_study_dt, "-srdt",
                 "--spatial_refinement_study_dt", "");
  args.AddOption(&temporal_refinement_study_sr, "-srsr",
                 "--temporal_refinement_study_sr", "");
  args.Parse();

  if (!args.Good()) {
    if (myid == 0) {
      args.PrintUsage(cout);
    }
    return 1;
  }
  if (myid == 0) {
    args.PrintOptions(cout);
  }

  if (myid == 0) {
    cout << endl;
  }

  //   //
  //   ************************************************************************************
  //   // Convergence study with respect to temporal refinement
  //   //
  //   ************************************************************************************

  //   if (myid == 0) {
  //     cout << "TEMPORAL CONVERGENCE STUDY" << endl;
  //   }

  //   for (int idx = 0; idx < temporal_refinement_levels; idx++) {
  //     runSimulation(std::pow(2, -idx), temporal_refinement_study_sr);
  //   }

  //   //
  //   ************************************************************************************
  //   // Convergence study with respect to spatial refinement
  //   //
  //   ************************************************************************************

  //   if (myid == 0) {
  //     cout << "SPATIAL CONVERGENCE STUDY" << endl;
  //   }

  //   for (int idx = 0; idx < spatial_refinement_levels; idx++) {
  //     runSimulation(spatial_refinement_study_dt, idx);
  //   }

  // ************************************************************************************
  // Combined convergence study
  // ************************************************************************************

  if (myid == 0) {
    cout << "COMBINED CONVERGENCE STUDY" << endl;
  }

  for (int idx = 0; idx < 7; idx++) {
    runSimulation(0.25 * std::pow(2, -idx), idx, order);
  }

  // runSimulation(0.1, 0, order);

  return 0;
}

// ************************************************************************************
//
// ************************************************************************************

void runSimulation(real_t dt, int sr, int order) {
  int myid = Mpi::WorldRank();

  const char *mesh_file = "../data/inline-quad.mesh";
  int par_ref_levels = 1;

  real_t t_final = 1;
  real_t speed = 1.;

  int precision = 8;
  cout.precision(precision);

  // 4. Read weights and poles
  // const char *filename01 = "../data/distExpGamma6_m_20_AAAtol_1.0e-10.csv";
  // const char *filename12 = "../data/distExpGamma5_m_20_AAAtol_1.0e-10.csv";

  const char *filename01 = "../data/distExpGamma6_m_45_AAAtol_1.0e-20.csv";
  const char *filename12 = "../data/distExpGamma5_m_44_AAAtol_1.0e-20.csv";

  // const char *filename01 = "../data/distExpGamma6_m_92_AAAtol_1.0e-40.csv";
  // const char *filename12 = "../data/distExpGamma5_m_91_AAAtol_1.0e-40.csv";

  std::vector<real_t> weights01, poles01, weights12, poles12;

  ReadWeightsPoles(weights01, poles01, filename01);
  ReadWeightsPoles(weights12, poles12, filename12);

  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();

  for (int lev = 0; lev < sr; lev++) {
    mesh->UniformRefinement();
  }

  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;
  for (int lev = 0; lev < par_ref_levels; lev++) {
    pmesh->UniformRefinement();
  }

  H1_FECollection fe_coll(order, dim);
  ParFiniteElementSpace fespace(pmesh, &fe_coll);

  HYPRE_BigInt fe_size = fespace.GlobalTrueVSize();
  if (myid == 0) {
    cout << endl;
    cout << "dofs: " << fe_size << endl;
    cout << "dt:   " << dt << endl;
  }

  ParGridFunction u_gf(&fespace);
  ParGridFunction dudt_gf(&fespace);
  ParGridFunction d2udt2_gf(&fespace);
  ParGridFunction error_gf(&fespace);

  // 6. Set the initial conditions for u. All boundaries are considered
  //    natural.
  FunctionCoefficient u_0(InitialSolution);
  u_gf.ProjectCoefficient(u_0);
  Vector u;
  u_gf.GetTrueDofs(u);

  FunctionCoefficient dudt_0(InitialRate);
  dudt_gf.ProjectCoefficient(dudt_0);
  Vector dudt;
  dudt_gf.GetTrueDofs(dudt);

  Vector d2udt2;
  d2udt2_gf.GetTrueDofs(d2udt2);

  Vector error;
  error_gf.GetTrueDofs(error);

  Array<int> block_trueOffsets01(weights01.size() + 1);
  block_trueOffsets01 = fespace.TrueVSize();
  block_trueOffsets01.PartialSum();
  BlockVector modes01(block_trueOffsets01);
  modes01 = 0.0;

  Array<int> block_trueOffsets12(weights12.size() + 1);
  block_trueOffsets12 = fespace.TrueVSize();
  block_trueOffsets12.PartialSum();
  BlockVector modes12(block_trueOffsets12);
  modes12 = 0.0;

  FunctionCoefficient rhs_fc(rhs_function);
  FunctionCoefficient sol_fc(solution);
  FunctionCoefficient rate_fc(rate);
  VectorFunctionCoefficient grad_sol_fc(dim, grad_solution);

  // 7. Initialize the wave operator and the visualization.
  Array<int> ess_bdr;
  if (pmesh->bdr_attributes.Size()) {
    ess_bdr.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr = 1;
  }

  DODiffusionWaveOperatorRIIA2 oper(fespace, ess_bdr, speed, weights01, poles01,
                                    weights12, poles12, rhs_fc);
  dudt_gf.SetFromTrueDofs(dudt);
  d2udt2_gf.SetFromTrueDofs(d2udt2);

  ParaViewDataCollection dataCollection(
      "conv_RIIA2" + std::to_string(sr) + "_" + std::to_string(dt), pmesh);
  dataCollection.SetPrefixPath("Paraview");
  dataCollection.SetLevelsOfDetail(1);
  dataCollection.SetDataFormat(VTKFormat::BINARY);
  dataCollection.SetHighOrderOutput(true);

  dataCollection.RegisterField("u", &u_gf);
  dataCollection.RegisterField("dudt", &dudt_gf);
  dataCollection.RegisterField("d2udt2", &d2udt2_gf);
  dataCollection.RegisterField("error", &error_gf);

  dataCollection.SetCycle(0);
  dataCollection.SetTime(0);
  dataCollection.Save();

  // Progress bar
  real_t pb = 0.0;
  bool last_pb = false;
  for (int ti = 1; !last_pb; ti++) {
    if (pb + dt >= t_final - dt / 2) {
      last_pb = true;
    }
    pb += dt;
    if (myid == 0) {
      cout << "I";
    }
  }
  if (myid == 0) {
    cout << endl;
  }

  // 8. Perform time-integration (looping over the time iterations, ti, with a
  //    time-step dt).

  real_t t = 0.0;
  bool last_step = false;

  real_t L2_inf{0.0};
  real_t L2_L1{0.0};
  real_t L2_L2{0.0};
  real_t L2error{0.0};
  real_t H1error{0.0};
  real_t rate_error{0.0};

  for (int ti = 1; !last_step; ti++) {
    if (t + dt >= t_final - dt / 2) {
      last_step = true;
    }
    if (myid == 0) {
      cout << "I" << flush;
    }
    oper.Step(u, dudt, d2udt2, modes01, modes12, t, dt);

    u_gf.SetFromTrueDofs(u);
    dudt_gf.SetFromTrueDofs(dudt);
    d2udt2_gf.SetFromTrueDofs(d2udt2);

    sol_fc.SetTime(t);
    rate_fc.SetTime(t);
    grad_sol_fc.SetTime(t);

    error_gf.ProjectCoefficient(sol_fc);
    error_gf -= u_gf;

    L2error = u_gf.ComputeL2Error(sol_fc);
    // cout << L2error << endl;
    L2_inf = max(L2_inf, L2error);
    L2_L1 += L2error * dt;
    L2_L2 += L2error * L2error * dt;
    rate_error = dudt_gf.ComputeL2Error(rate_fc);

    // dataCollection.SetCycle(ti);
    // dataCollection.SetTime(t);
    // dataCollection.Save();

    // if (myid == 0) {
    //   cout << "H1:   " << H1error << endl;
    // }
  }
  if (myid == 0) {
    cout << endl;
    cout << "L2_inf:   " << L2_inf << endl;
    cout << "L2_L1:   " << L2_L1 << endl;
    cout << "L2_L2:   " << sqrt(L2_L2) << endl;
    // cout << "L2_rate:   " << rate_error << endl;
  }

  delete pmesh;
}