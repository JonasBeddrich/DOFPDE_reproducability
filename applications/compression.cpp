#include <fstream>
#include <iostream>
#include <memory>

#include "DOHyperelasticOperator.hpp"
#include "mfem.hpp"
#include "utils.cpp"

using namespace std;
using namespace mfem;

void InitialDeformation(const Vector &x, Vector &y);

void InitialVelocity(const Vector &x, Vector &v);

void rhs_function(const Vector &x, real_t t, Vector &y);

void circle_function(const Vector &x, real_t t, Vector &y);

void TUM_function(const Vector &x, real_t t, Vector &y);

void bunny_function(const Vector &x, real_t t, Vector &y);

void ammonite_function(const Vector &x, real_t t, Vector &y);

int main(int argc, char *argv[]) {
  // 1. Initialize MPI and HYPRE.
  Mpi::Init(argc, argv);
  int myid = Mpi::WorldRank();
  Hypre::Init();

  // 2. Parse command-line options.
  const char *mesh_file = "../data/ref-cube.mesh";
  int ser_ref_levels = 2;
  int par_ref_levels = 0;
  int order = 2;

  int scenario = 0;

  real_t t_final = 3.;
  real_t dt = 0.001;
  real_t visc = 1;
  real_t mu = 0.1;
  real_t K = 1;

  real_t X = 1.;
  real_t Y = 1.;
  real_t Z = 0.5;

  bool adaptive_lin_rtol = true;
  bool visualization = true;
  int vis_steps = 10;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                 "Number of times to refine the mesh uniformly in serial.");
  args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                 "Number of times to refine the mesh uniformly in parallel.");
  args.AddOption(&order, "-o", "--order",
                 "Order (degree) of the finite elements.");
  args.AddOption(&scenario, "-sc", "--scenario", "Choose scenario");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&visc, "-v", "--viscosity", "Viscosity coefficient.");
  args.AddOption(&mu, "-mu", "--shear-modulus",
                 "Shear modulus in the Neo-Hookean hyperelastic model.");
  args.AddOption(&K, "-K", "--bulk-modulus",
                 "Bulk modulus in the Neo-Hookean hyperelastic model.");
  args.AddOption(&adaptive_lin_rtol, "-alrtol", "--adaptive-lin-rtol",
                 "-no-alrtol", "--no-adaptive-lin-rtol",
                 "Enable or disable adaptive linear solver rtol.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                 "Visualize every n-th timestep.");
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

  const char *filename;
  if (scenario == 0) {
    filename = "../data/dist01TTPalpha_m_44_AAAtol_1.0e-20.csv";
  } else if (scenario == 1) {
    filename = "../data/dist05TTPalpha_m_44_AAAtol_1.0e-20.csv";
  } else if (scenario == 2) {
    filename = "../data/dist2TTPalpha_m_44_AAAtol_1.0e-20.csv";
  } else if (scenario == 3) {
    filename = "../data/RA/distr01Bump_center_1.0_m_20_AAAtol_5.0e-15.csv";
  }

  std::vector<real_t> weights, poles;
  ReadWeightsPoles(weights, poles, filename);

  // 3. Read the serial mesh from the given mesh file on all processors. We can
  //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
  //    with the same code.
  // Mesh *mesh = new Mesh(mesh_file, 1, 1);
  Mesh mesh = Mesh::MakeCartesian3D(100, 100, 50, Element::HEXAHEDRON, X, Y, Z);
  int dim = mesh.Dimension();

  // 5. Refine the mesh in serial to increase the resolution. In this example
  //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
  //    a command-line parameter.
  for (int lev = 0; lev < ser_ref_levels; lev++) {
    mesh.UniformRefinement();
  }

  // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
  //    this mesh further in parallel to increase the resolution. Once the
  //    parallel mesh is defined, the serial mesh can be deleted.
  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
  for (int lev = 0; lev < par_ref_levels; lev++) {
    pmesh->UniformRefinement();
  }

  // 7. Define the parallel vector finite element spaces representing the mesh
  //    deformation x_gf, the velocity v_gf, and the initial configuration,
  //    x_ref. Define also the elastic energy density, w_gf, which is in a
  //    discontinuous higher-order space. Since x and v are integrated in time
  //    as a system, we group them together in block vector vx, on the unique
  //    parallel degrees of freedom, with offsets given by array true_offset.
  H1_FECollection fe_coll(order, dim);
  ParFiniteElementSpace fespace(pmesh, &fe_coll, dim);

  HYPRE_BigInt glob_size = fespace.GlobalTrueVSize();
  if (myid == 0) {
    cout << "Number of velocity/deformation unknowns: " << glob_size << endl;
  }

  int true_size = fespace.TrueVSize();
  Array<int> true_offset(3);
  true_offset[0] = 0;
  true_offset[1] = true_size;
  true_offset[2] = 2 * true_size;

  BlockVector vx(true_offset);
  ParGridFunction v_gf, x_gf;
  v_gf.MakeTRef(&fespace, vx, true_offset[0]);
  x_gf.MakeTRef(&fespace, vx, true_offset[1]);

  Array<int> block_trueOffsets(weights.size() + 1);
  block_trueOffsets = fespace.TrueVSize();
  block_trueOffsets.PartialSum();
  BlockVector modes(block_trueOffsets);
  modes = 0.0;

  ParGridFunction x_ref(&fespace);
  pmesh->GetNodes(x_ref);

  L2_FECollection w_fec(order + 1, dim);
  ParFiniteElementSpace w_fespace(pmesh, &w_fec);
  ParGridFunction w_gf(&w_fespace);

  // 8. Set the initial conditions for v_gf, x_gf and vx, and define the
  //    boundary conditions on a beam-like mesh (see description above).
  VectorFunctionCoefficient velo(dim, InitialVelocity);
  v_gf.ProjectCoefficient(velo);
  v_gf.SetTrueVector();
  v_gf.SetFromTrueVector();

  VectorFunctionCoefficient deform(dim, InitialDeformation);
  x_gf.ProjectCoefficient(deform);
  x_gf.SetTrueVector();
  x_gf.SetFromTrueVector();

  Array<int> ess_bdr(fespace.GetMesh()->bdr_attributes.Max());
  ess_bdr = 0;
  ess_bdr[0] = 1;

  VectorFunctionCoefficient rhs_fc(dim, rhs_function);
  VectorFunctionCoefficient bd_fc(dim, ammonite_function);

  // 9. Initialize the hyperelastic operator, the GLVis visualization and print
  //    the initial energies.
  DOHyperelasticOperator oper(fespace, ess_bdr, visc, mu, K, weights, poles,
                              rhs_fc, bd_fc);

  real_t ee0 = oper.ElasticEnergy(x_gf);
  real_t ke0 = oper.KineticEnergy(v_gf);
  if (myid == 0) {
    cout << "initial elastic energy (EE) = " << ee0 << endl;
    cout << "initial kinetic energy (KE) = " << ke0 << endl;
    cout << "initial   total energy (TE) = " << (ee0 + ke0) << endl;
  }

  std::string paraviewname = "compression";
  paraviewname += "_sc_" + std::to_string(scenario);

  ParaViewDataCollection dataCollection(paraviewname, pmesh);
  dataCollection.SetPrefixPath("Paraview");
  dataCollection.SetLevelsOfDetail(order);
  dataCollection.SetDataFormat(VTKFormat::BINARY);
  dataCollection.SetHighOrderOutput(true);

  dataCollection.RegisterField("x", &x_gf);
  dataCollection.RegisterField("v", &v_gf);

  dataCollection.SetCycle(0);
  dataCollection.SetTime(0);
  dataCollection.Save();

  real_t t = 0.0;
  oper.SetTime(t);

  // 10. Perform time-integration
  //     (looping over the time iterations, ti, with a time-step dt).
  bool last_step = false;
  for (int ti = 1; !last_step; ti++) {
    real_t dt_real = min(dt, t_final - t);

    //   ode_solver->Step(vx, t, dt_real);

    oper.SetTime(t + dt);
    oper.IEStep(vx, modes, t, dt_real);
    t += dt;

    v_gf.SetFromTrueVector();
    x_gf.SetFromTrueVector();

    last_step = (t >= t_final - 1e-8 * dt);

    if (last_step || (ti % vis_steps) == 0) {
      real_t ee = oper.ElasticEnergy(x_gf);
      real_t ke = oper.KineticEnergy(v_gf);

      dataCollection.SetCycle(ti);
      dataCollection.SetTime(t);
      dataCollection.Save();

      if (myid == 0) {
        cout << "step " << ti << ", t = " << t << ", EE = " << ee
             << ", KE = " << ke << ", ΔTE = " << (ee + ke) - (ee0 + ke0)
             << endl;
      }
    }
  }

  // 11. Save the displaced mesh, the velocity and elastic energy.
  {
    v_gf.SetFromTrueVector();
    x_gf.SetFromTrueVector();
    GridFunction *nodes = &x_gf;
    int owns_nodes = 0;
    pmesh->SwapNodes(nodes, owns_nodes);

    ostringstream mesh_name, velo_name, ee_name;
    mesh_name << "deformed." << setfill('0') << setw(6) << myid;
    velo_name << "velocity." << setfill('0') << setw(6) << myid;
    ee_name << "elastic_energy." << setfill('0') << setw(6) << myid;

    ofstream mesh_ofs(mesh_name.str().c_str());
    mesh_ofs.precision(8);
    pmesh->Print(mesh_ofs);
    pmesh->SwapNodes(nodes, owns_nodes);
    ofstream velo_ofs(velo_name.str().c_str());
    velo_ofs.precision(8);
    v_gf.Save(velo_ofs);
    ofstream ee_ofs(ee_name.str().c_str());
    ee_ofs.precision(8);
    oper.GetElasticEnergyDensity(x_gf, w_gf);
    w_gf.Save(ee_ofs);
  }

  // 12. Free the used memory.
  delete pmesh;

  return 0;
}

void InitialDeformation(const Vector &x, Vector &y) {
  y = 0.;
  y += x;
}

void InitialVelocity(const Vector &x, Vector &v) { v = 0.; }

void rhs_function(const Vector &x, real_t t, Vector &y) { y = 0.; }

void circle_function(const Vector &x, real_t t, Vector &y) {
  if (x[2] > 0.2 - 1e-6 && t < 1. &&
      (x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) < 0.2) {
    y[2] = -1.;
  }
}

void TUM_function(const Vector &x, real_t t, Vector &y) {
  if (x[2] > 0.5 - 1e-8 && t < 1.) {
    if (x[0] > 0.25 && x[0] < 0.45 && x[1] > 0.4 && x[1] < 0.45) {
      y[2] = -1.;
    } else if (x[0] > 0.5 && x[0] < 0.75 && x[1] > 0.4 && x[1] < 0.45) {
      y[2] = -1.;
    } else if (x[0] > 0.3 && x[0] < 0.35 && x[1] > 0.45 && x[1] < 0.6) {
      y[2] = -1.;
    } else if (x[0] > 0.4 && x[0] < 0.45 && x[1] > 0.45 && x[1] < 0.6) {
      y[2] = -1.;
    } else if (x[0] > 0.5 && x[0] < 0.55 && x[1] > 0.45 && x[1] < 0.6) {
      y[2] = -1.;
    } else if (x[0] > 0.6 && x[0] < 0.65 && x[1] > 0.45 && x[1] < 0.6) {
      y[2] = -1.;
    } else if (x[0] > 0.7 && x[0] < 0.75 && x[1] > 0.45 && x[1] < 0.6) {
      y[2] = -1.;
    } else if (x[0] > 0.45 && x[0] < 0.5 && x[1] > 0.55 && x[1] < 0.6) {
      y[2] = -1.;
    }
  }
}
void bunny_function(const Vector &x, real_t t, Vector &y) {
  if (x[2] > 0.5 - 1e-8 && t < 1.) {
    if ((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.3) * (x[1] - 0.3) < 0.2 * 0.2) {
      y[2] = -1.;
    } else if ((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.6) * (x[1] - 0.6) <
               0.125 * 0.125) {
      y[2] = -1.;
    } else if (((x[0] - 0.6) * cos(1.25) + (x[1] - 0.8) * sin(1.25)) *
                       ((x[0] - 0.6) * cos(1.25) + (x[1] - 0.8) * sin(1.25)) /
                       0.02 +
                   ((x[1] - 0.8) * cos(1.25) - (x[0] - 0.6) * sin(1.25)) *
                       ((x[1] - 0.8) * cos(1.25) - (x[0] - 0.6) * sin(1.25)) /
                       0.002 <
               1) {
      y[2] = -1.;
    } else if (((x[0] - 0.4) * cos(1.88) + (x[1] - 0.8) * sin(1.88)) *
                       ((x[0] - 0.4) * cos(1.88) + (x[1] - 0.8) * sin(1.88)) /
                       0.02 +
                   ((x[1] - 0.8) * cos(1.88) - (x[0] - 0.4) * sin(1.88)) *
                       ((x[1] - 0.8) * cos(1.88) - (x[0] - 0.4) * sin(1.88)) /
                       0.002 <
               1) {
      y[2] = -1.;
    }
  }
}

void ammonite_function(const Vector &coord, real_t t, Vector &f) {
  if (coord[2] > 0.5 - 1e-8 && t < 2.) {
    real_t a = 0.025;
    real_t b = 0.1;

    real_t x_ = coord[0] + 0.019 - 0.5;
    real_t y_ = coord[1] + 0.053 - 0.5;
    real_t x = x_ * cos(M_PI / 4) - y_ * sin(M_PI / 4);
    real_t y = x_ * sin(M_PI / 4) + y_ * cos(M_PI / 4);

    real_t r = std::sqrt(x * x + y * y);
    real_t theta = std::atan2(y, x);

    int outer;
    bool flag_outer = false;
    for (int i = 0; i < 5; ++i) {
      if (!flag_outer && r < a * std::exp(b * (theta + i * 2. * M_PI))) {
        flag_outer = true;
        outer = i;
      }
    }

    real_t r_outer = a * std::exp(b * (theta + outer * 2. * M_PI));
    real_t r_inner = 0;
    if (outer > 0) {
      r_inner = a * std::exp(b * (theta + (outer - 1) * 2. * M_PI));
    }
    real_t r_center = (r_outer + r_inner) / 2.;

    real_t force = -std::sqrt((r_outer - r_inner) * (r_outer - r_inner) / 4. -
                              (r - r_center) * (r - r_center));
    force *= (10 + 3 * cos((11 * M_PI - theta - 2 * outer * M_PI) *
                           (11 * M_PI - theta - 2 * outer * M_PI)));
    if (std::isnan(force)) {
      force = 0.;
    }
    if (t < 1.) {
      f[2] = force * t;
    } else {
      f[2] = force;
    }
  }
}