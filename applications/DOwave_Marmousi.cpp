#include <fstream>
#include <iostream>

#include "do03spatialdiffusionwaveoperator.hpp"
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class CSVInterpolatedCoefficient : public Coefficient {
 private:
  std::vector<std::vector<double>> data;
  int N;
  double dx;

 public:
  CSVInterpolatedCoefficient(const std::string &filename) {
    // Read CSV into data array
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
      std::stringstream ss(line);
      std::string cell;
      std::vector<double> row;

      while (std::getline(ss, cell, ',')) {
        row.push_back(std::stod(cell));
      }

      if (!data.empty() && row.size() != data[0].size()) {
        mfem::err << "Inconsistent row length in CSV file." << endl;
        MFEM_ABORT("Invalid CSV format");
      }

      data.push_back(row);
    }

    N = data.size();
    if (N < 2 || data[0].size() != N) {
      mfem::err << "CSV must be a square grid with at least 2x2 values."
                << endl;
      MFEM_ABORT("Invalid CSV dimensions");
    }

    dx = 1.0 / (N - 1);
  }

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
    Vector trans_ip;
    T.Transform(ip, trans_ip);
    double x = trans_ip(0);
    double y = trans_ip(1);

    // Clamp to unit square
    x = std::max(0.0, std::min(1.0, x));
    y = std::max(0.0, std::min(1.0, y));

    // Determine the lower-left grid index
    int i = std::min(int(x / dx), N - 2);
    int j = std::min(int(y / dx), N - 2);

    double x0 = i * dx;
    double y0 = j * dx;

    // Compute weights
    double tx = (x - x0) / dx;
    double ty = (y - y0) / dx;

    // Bilinear interpolation
    double val = (1 - tx) * (1 - ty) * data[j][i] +
                 tx * (1 - ty) * data[j][i + 1] +
                 (1 - tx) * ty * data[j + 1][i] + tx * ty * data[j + 1][i + 1];

    return val;
  }
};

real_t InitialSolution_sc0(const Vector &x) {
  Vector tmp(2);
  real_t result = 0.;

  tmp[0] = x[0] - 0.63;
  tmp[1] = x[1] - 0.37;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.5;
  tmp[1] = x[1] - 0.167;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.65;
  tmp[1] = x[1] - 0.1;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.8;
  tmp[1] = x[1] - 0.55;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.2;
  tmp[1] = x[1] - 0.36;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  return result;
}

real_t InitialSolution_sc1(const Vector &x) {
  Vector tmp(2);
  real_t result = 0.;

  tmp[0] = x[0] - 0.47;
  tmp[1] = x[1] - 0.36;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.72;
  tmp[1] = x[1] - 0.63;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  return result;
}

real_t InitialSolution_sc2(const Vector &x) {
  Vector tmp(2);
  real_t result = 0.;

  tmp[0] = x[0] - 0.25;
  tmp[1] = x[1] - 0.25;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.25;
  tmp[1] = x[1] - 0.75;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.75;
  tmp[1] = x[1] - 0.25;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  tmp[0] = x[0] - 0.75;
  tmp[1] = x[1] - 0.75;
  result += exp(-1000 * tmp.Norml2() * tmp.Norml2());

  return result;
}

real_t InitialRate(const Vector &x) { return 0.0; }

real_t InitialAcceleration(const Vector &x) { return 0.0; }

real_t rhs_function(const Vector &x, real_t t) { return 0; }

int main(int argc, char *argv[]) {
  Mpi::Init(argc, argv);
  int num_procs = Mpi::WorldSize();
  int myid = Mpi::WorldRank();
  Hypre::Init();

  const char *mesh_file = "../data/inline-quad.mesh";
  int ser_ref_levels = 4;
  int par_ref_levels = 4;

  int scenario = 2;
  int m = 20;

  int order = 2;
  real_t t_final = 3;
  real_t dt = 1.0e-2;
  real_t speed = 0.2;

  bool dirichlet = false;
  int vis_steps = 1;

  int precision = 8;
  cout.precision(precision);

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&ser_ref_levels, "-sr", "--serial-refine",
                 "Number of times to refine the mesh uniformly.");
  args.AddOption(&par_ref_levels, "-pr", "--parallel-refine",
                 "Number of times to refine the mesh uniformly.");
  args.AddOption(
      &scenario, "-sc", "--scenario",
      "Choose the scenario: 0 - Marmousi, 1 - Random field, 2 - toy problem.");
  args.AddOption(&order, "-o", "--order",
                 "Order (degree) of the finite elements.");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&speed, "-c", "--speed", "Wave speed.");
  args.AddOption(&dirichlet, "-dir", "--dirichlet", "-neu", "--neumann",
                 "BC switch.");
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
  std::vector<Coefficient *> weights01_coeff, poles01_coeff, weights12_coeff,
      poles12_coeff, weights23_coeff, poles23_coeff;

  std::string path;
  if (scenario == 0) {
    path = "../data/Marmousi/weights_n_poles_sc0/";
  } else if (scenario == 1) {
    path = "../data/Marmousi/weights_n_poles_sc1/";
  } else {
    scenario = 2;
    path = "../data/Marmousi/weights_n_poles_sc2/";
  }

  for (int i = 0; i < m; i++) {
    weights01_coeff.push_back(new CSVInterpolatedCoefficient(
        path + "weights01_" + std::to_string(i) + ".csv"));
    weights12_coeff.push_back(new CSVInterpolatedCoefficient(
        path + "weights12_" + std::to_string(i) + ".csv"));
    weights23_coeff.push_back(new CSVInterpolatedCoefficient(
        path + "weights23_" + std::to_string(i) + ".csv"));
    poles01_coeff.push_back(new CSVInterpolatedCoefficient(
        path + "poles01_" + std::to_string(i) + ".csv"));
    poles12_coeff.push_back(new CSVInterpolatedCoefficient(
        path + "poles12_" + std::to_string(i) + ".csv"));
    poles23_coeff.push_back(new CSVInterpolatedCoefficient(
        path + "poles23_" + std::to_string(i) + ".csv"));
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

  ParGridFunction alpha_gf(&fespace);
  CSVInterpolatedCoefficient alpha_coeff("../data/Marmousi/alpha_sc" +
                                         std::to_string(scenario) + ".csv");
  alpha_gf.ProjectCoefficient(alpha_coeff);

  ParGridFunction u_gf(&fespace);
  ParGridFunction dudt_gf(&fespace);
  ParGridFunction d2udt2_gf(&fespace);
  ParGridFunction d3udt3_gf(&fespace);

  // 6. Set the initial conditions for u. All boundaries are considered
  //    natural.

  if (scenario == 0) {
    FunctionCoefficient u_0(InitialSolution_sc0);
    u_gf.ProjectCoefficient(u_0);
  } else if (scenario == 1) {
    FunctionCoefficient u_0(InitialSolution_sc1);
    u_gf.ProjectCoefficient(u_0);
  } else {
    FunctionCoefficient u_0(InitialSolution_sc2);
    u_gf.ProjectCoefficient(u_0);
  }
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

  Array<int> block_trueOffsets(m + 1);
  block_trueOffsets = fespace.TrueVSize();
  block_trueOffsets.PartialSum();

  BlockVector modes01(block_trueOffsets);
  BlockVector modes12(block_trueOffsets);
  BlockVector modes23(block_trueOffsets);

  modes01 = 0.0;
  modes12 = 0.0;
  modes23 = 0.0;
  FunctionCoefficient rhs_fc(rhs_function);

  // 7. Initialize the wave operator and the visualization.
  Array<int> ess_bdr;
  if (pmesh->bdr_attributes.Size()) {
    ess_bdr.SetSize(pmesh->bdr_attributes.Max());

    if (dirichlet) {
      ess_bdr = 1;
    } else {
      ess_bdr = 0;
    }
  }
  DO03SpatialDiffusionWaveOperator oper(
      fespace, ess_bdr, speed, weights01_coeff, poles01_coeff, weights12_coeff,
      poles12_coeff, weights23_coeff, poles23_coeff, rhs_fc, dt);

  u_gf.SetFromTrueDofs(u);
  dudt_gf.SetFromTrueDofs(dudt);
  d2udt2_gf.SetFromTrueDofs(d2udt2);
  d3udt3_gf.SetFromTrueDofs(d3udt3);

  std::string filename = "Marmousi";
  filename += std::to_string(scenario);
  filename += "_c_" + std::to_string(speed);

  ParaViewDataCollection dataCollection(filename, pmesh);
  dataCollection.SetPrefixPath("Paraview");
  dataCollection.SetLevelsOfDetail(1);
  dataCollection.SetDataFormat(VTKFormat::BINARY);
  dataCollection.SetHighOrderOutput(true);

  dataCollection.RegisterField("u", &u_gf);
  dataCollection.RegisterField("dudt", &dudt_gf);
  dataCollection.RegisterField("d2udt2", &d2udt2_gf);
  dataCollection.RegisterField("d3udt3", &d3udt3_gf);
  dataCollection.RegisterField("alpha", &alpha_gf);

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

    oper.StepIE(u, dudt, d2udt2, d3udt3, modes01, modes12, modes23, t, dt);

    u_gf.SetFromTrueDofs(u);
    dudt_gf.SetFromTrueDofs(dudt);
    d2udt2_gf.SetFromTrueDofs(d2udt2);
    d3udt3_gf.SetFromTrueDofs(d3udt3);

    if (last_step || (ti % vis_steps) == 0) {
      if (myid == 0) {
        cout << "step " << ti << ", t = " << t << endl;
      }
      dataCollection.SetCycle(ti);
      dataCollection.SetTime(t);
      dataCollection.Save();
    }
  }

  // Free the used memory.
  delete pmesh;
  return 0;
}