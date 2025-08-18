#include <fstream>
#include <iostream>

#include "dodiffusionwaveoperator.hpp"
#include "mfem.hpp"

using namespace std;
using namespace mfem;

real_t InitialSolution(const Vector &x) { return 0.; }

real_t InitialRate(const Vector &x) { return 0.0; }

real_t InitialAcceleration(const Vector &x) { return 0.0; }

real_t rhs_function(const Vector &x, real_t t) {
  if (t < 0.1) {
    return 100 * sin(20 * M_PI * t) * exp(-1 / (10 * t * (1 - 10 * t)));
  }
  return 0.;
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

std::string GetDistribution(const char *filename) {
  std::string tmp(filename);
  int start = tmp.find("dist");
  int end = tmp.find("_m_");
  std::string result = tmp.substr(start + 4, end - start - 4);
  return result;
}

int main(int argc, char *argv[]) {
  Mpi::Init(argc, argv);
  int num_procs = Mpi::WorldSize();
  int myid = Mpi::WorldRank();
  Hypre::Init();

  const char *mesh_file = "../data/Starnbergersee/starnbergersee3D.msh";
  int scenario = 0;
  int ser_ref_levels = 1;
  int par_ref_levels = 1;

  int order = 2;
  real_t t_final = 1;
  real_t dt = 1.0e-2;
  real_t speed = 10;

  int vis_steps = 1;

  int precision = 8;
  cout.precision(precision);

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&scenario, "-sc", "--scenario", "Scenario of the simulation.");
  args.AddOption(&ser_ref_levels, "-sr", "--serial-refine",
                 "Number of times to refine the mesh uniformly.");
  args.AddOption(&par_ref_levels, "-pr", "--parallel-refine",
                 "Number of times to refine the mesh uniformly.");
  args.AddOption(&order, "-o", "--order",
                 "Order (degree) of the finite elements.");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&speed, "-c", "--speed", "Wave speed.");
  args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                 "Visualize every n-th timestep.");
  args.Parse();

  if (!args.Good()) {
    args.PrintUsage(cout);
    return 1;
  }
  if (myid == 0) {
    args.PrintOptions(cout);
  }

  // 4. Read weights and poles
  const char *filename01 = "../data/RA/distZero_m_.csv";
  const char *filename12;
  if (scenario == 0) {
    filename12 = "../data/RA/distr01Bump_center_1.0_m_20_AAAtol_5.0e-15.csv";
  } else if (scenario == 1) {
    filename12 = "../data/RA/distr05Bump_center_1.0_m_20_AAAtol_5.0e-15.csv";
  } else {
    filename12 = "../data/Starnbergersee/uniform_m_20_AAAtol_5.0e-15.csv";
  }

  std::vector<real_t> weights01, poles01, weights12, poles12, weights23,
      poles23;

  ReadWeightsPoles(weights01, poles01, filename01);
  ReadWeightsPoles(weights12, poles12, filename12);

  // normalize weights as the distribution is defined from (c-r, c+r)
  // but we use (c-r, c)
  if (scenario == 0 || scenario == 1) {
    for (real_t &weight : weights01) {
      weight *= 2.;
    }
    for (real_t &weight : weights12) {
      weight *= 2.;
    }
  }

  // 2. Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral and hexahedral meshes with the same code.
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();

  // 3. Refine the mesh to increase the resolution. In this example we do
  //    'ref_levels' of uniform refinement, where 'ref_levels' is a
  //    command-line parameter.
  for (int lev = 0; lev < ser_ref_levels; lev++) {
    mesh->UniformRefinement();
  }

  // parallel mesh
  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;
  for (int lev = 0; lev < par_ref_levels; lev++) {
    pmesh->UniformRefinement();
  }

  // 5. Define the vector finite element space representing the current and the
  //    initial temperature, u_ref.
  H1_FECollection fe_coll(order, dim);
  ParFiniteElementSpace fespace(pmesh, &fe_coll);

  HYPRE_BigInt fe_size = fespace.GlobalTrueVSize();
  if (myid == 0) {
    cout << "Number of spatial unknowns: " << fe_size << endl;
  }

  ParGridFunction u_gf(&fespace);
  ParGridFunction dudt_gf(&fespace);
  ParGridFunction d2udt2_gf(&fespace);
  ParGridFunction d3udt3_gf(&fespace);

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

  FunctionCoefficient d2udt2_0(InitialAcceleration);
  d2udt2_gf.ProjectCoefficient(d2udt2_0);
  Vector d2udt2;
  d2udt2_gf.GetTrueDofs(d2udt2);

  Vector d3udt3;
  d3udt3_gf.GetTrueDofs(d3udt3);

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

  // 7. Initialize the wave operator and the visualization.
  Array<int> ess_bdr;
  ess_bdr.SetSize(pmesh->bdr_attributes.Max());
  // cout << "bdr attributes max" << pmesh->bdr_attributes.Max() << endl;
  ess_bdr[0] =
      0;  // I have no idea, surface, sides, bottom and island are not it
  ess_bdr[1] = 0;  // surface
  ess_bdr[2] = 1;  // sides
  ess_bdr[3] = 1;  // bottom

  DODiffusionWaveOperator oper(fespace, ess_bdr, speed, weights01, poles01,
                               weights12, poles12, rhs_fc);

  u_gf.SetFromTrueDofs(u);
  dudt_gf.SetFromTrueDofs(dudt);
  d2udt2_gf.SetFromTrueDofs(d2udt2);

  std::string filename = "Starnbergersee";
  filename += "_c_" + std::to_string((int)speed);
  filename += "_sc_" + std::to_string(scenario);

  ParaViewDataCollection dataCollection(filename, pmesh);
  dataCollection.SetPrefixPath("Paraview");
  dataCollection.SetLevelsOfDetail(1);
  dataCollection.SetDataFormat(VTKFormat::BINARY);
  dataCollection.SetHighOrderOutput(true);

  dataCollection.RegisterField("u", &u_gf);
  dataCollection.RegisterField("dudt", &dudt_gf);
  dataCollection.RegisterField("d2udt2", &d2udt2_gf);

  dataCollection.SetCycle(0);
  dataCollection.SetTime(0);
  dataCollection.Save();

  // 8. Perform time-integration (looping over the time iterations, ti, with a
  //    time-step dt).

  real_t t = 0.0;
  bool last_step = false;

  for (int ti = 1; !last_step; ti++) {
    if (t + dt >= t_final - dt / 2) {
      last_step = true;
    }

    oper.StepIE(u, dudt, d2udt2, modes01, modes12, t, dt);

    u_gf.SetFromTrueDofs(u);
    dudt_gf.SetFromTrueDofs(dudt);
    d2udt2_gf.SetFromTrueDofs(d2udt2);

    if (myid == 0) {
      cout << "step " << ti << ", t = " << t << endl;
    }

    if (last_step || (ti % vis_steps) == 0) {
      dataCollection.SetCycle(ti);
      dataCollection.SetTime(t);
      dataCollection.Save();
    }
  }

  // Free the used memory.
  delete pmesh;

  return 0;
}
