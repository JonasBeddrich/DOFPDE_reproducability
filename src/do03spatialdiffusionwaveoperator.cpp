#include "mfem.hpp"
#include "do03spatialdiffusionwaveoperator.hpp"
#include "QuotientCoefficient.cpp"

#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

DO03SpatialDiffusionWaveOperator::DO03SpatialDiffusionWaveOperator(
   ParFiniteElementSpace &f, 
   Array<int> &ess_bdr, 
   real_t speed, 
   std::vector<Coefficient*> &w01, 
   std::vector<Coefficient*> &p01, 
   std::vector<Coefficient*> &w12, 
   std::vector<Coefficient*> &p12, 
   std::vector<Coefficient*> &w23, 
   std::vector<Coefficient*> &p23, 
   FunctionCoefficient &rhs, 
   real_t dt_
): 
   fespace(f), 
   height(fespace.GetTrueVSize()), 
   weights01(w01), 
   poles01(p01), 
   weights12(w12), 
   poles12(p12),
   weights23(w23), 
   poles23(p23),
   rhs_fc(rhs),  
   M(NULL), 
   K(NULL), 
   T(NULL), 
   t(0.0), 
   dt(dt_), 
   M_solver(f.GetComm()), 
   T_solver(f.GetComm()),
   zn(height), 
   tmp(height)
{

const real_t rel_tol = 1e-12;
c2 = new ConstantCoefficient(speed*speed);
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
T_prec.SetPrintLevel(0); 
T_solver.SetPreconditioner(T_prec); 

ComputeCoefficientMatrices(); 

}

void DO03SpatialDiffusionWaveOperator::StepIE(Vector &u, Vector &du_dt, Vector &d2u_dt2, Vector &d3u_dt3, BlockVector &modes01, BlockVector &modes12, BlockVector &modes23, real_t &t, real_t &dt){
   // Compute all the stuff that is left 
   ComputeRightHandSide(u, du_dt, d2u_dt2, modes01, modes12, modes23, t, dt); 
   // Solve d3u_dt3 
   SolveImplicitProblem(d3u_dt3, dt); 
   // Compute u and dudt 
   UpdateSolution(u, du_dt, d2u_dt2, d3u_dt3, dt); 
   // Update the fractional modes 
   UpdateModes(u, du_dt, d2u_dt2, d3u_dt3, modes01, modes12, modes23, dt); 
   // Update time 
   t += dt; 
}

void DO03SpatialDiffusionWaveOperator::ComputeRightHandSide(const Vector &u, const Vector &du_dt, const Vector &d2u_dt2, 
   const BlockVector &modes01, const BlockVector &modes12, const BlockVector &modes23, const real_t &t, const real_t &dt){

   ParGridFunction rhs_gf(&fespace); 
   rhs_fc.SetTime(t); 
   rhs_gf.ProjectCoefficient(rhs_fc); 
   Vector rhs;
   rhs_gf.GetTrueDofs(rhs); 
   Mmat.Mult(rhs, zn); 
   
   for (size_t i = 0; i < poles01.size(); ++i) {
      BETAmat01[i].AddMult(du_dt, zn);    
   }
   for (size_t i = 0; i < poles01.size(); ++i) {
      BETAmat01[i].AddMult(d2u_dt2, zn, dt);    
   }
   for (size_t i = 0; i < poles12.size(); ++i) {
      BETAmat12[i].AddMult(d2u_dt2, zn);    
   }

   for (size_t i = 0; i < poles01.size(); ++i) {
      GAMMAmat01[i].AddMult(modes01.GetBlock(i), zn);    
   }
   for (size_t i = 0; i < poles12.size(); ++i) {
      GAMMAmat12[i].AddMult(modes12.GetBlock(i), zn);    
   }
   for (size_t i = 0; i < poles23.size(); ++i) {
      GAMMAmat23[i].AddMult(modes23.GetBlock(i), zn);    
   }
   zn.Neg(); 
   
   Kmat.AddMult(u, zn, -1.0);
   Kmat.AddMult(du_dt, zn, -dt);
   Kmat.AddMult(d2u_dt2, zn, -dt*dt);
   zn.SetSubVector(ess_tdof_list, 0.0); 
}

void DO03SpatialDiffusionWaveOperator::SolveImplicitProblem(Vector &d3u_dt3, const real_t &dt){
   if (!T)
   {
   T = Add(1.0, SIPmat, dt * dt * dt, Kmat);
   T_solver.SetOperator(*T);
   }
   T_solver.Mult(zn, d3u_dt3);
   d3u_dt3.SetSubVector(ess_tdof_list, 0.0);
}

void DO03SpatialDiffusionWaveOperator::UpdateSolution(Vector &u, Vector &du_dt, Vector &d2u_dt2, const Vector &d3u_dt3, const real_t &dt){
   d2u_dt2.Add(dt, d3u_dt3); 
   du_dt.Add(dt, d2u_dt2); 
   u.Add(dt, du_dt); 
}

void DO03SpatialDiffusionWaveOperator::UpdateModes(const Vector &u, const Vector &du_dt, const Vector &d2u_dt2, const Vector &d3u_dt3, 
   BlockVector &modes01, BlockVector &modes12, BlockVector &modes23, real_t &dt){
   // Update the 01 modes 
   for (size_t i = 0; i < weights01.size(); ++i) {
      GAMMAmat01[i].Mult(modes01.GetBlock(i), tmp); 
      BETAmat01[i].AddMult(du_dt, tmp); 
      M_solver.Mult(tmp, modes01.GetBlock(i)); 
   } 
   // Update the 12 modes 
   for (size_t i = 0; i < weights12.size(); ++i) {
      GAMMAmat12[i].Mult(modes12.GetBlock(i), tmp); 
      BETAmat12[i].AddMult(d2u_dt2, tmp); 
      M_solver.Mult(tmp, modes12.GetBlock(i)); 
   } 
   // Update the 23 modes 
   for (size_t i = 0; i < weights23.size(); ++i) {
      GAMMAmat23[i].Mult(modes23.GetBlock(i), tmp); 
      BETAmat23[i].AddMult(d3u_dt3, tmp); 
      M_solver.Mult(tmp, modes23.GetBlock(i)); 
   } 
}

void DO03SpatialDiffusionWaveOperator::ComputeCoefficientMatrices(){
   HypreParMatrix TMPmat; 
   ConstantCoefficient con(1.0); 
   
   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator(con));
   M->Assemble(0);
   M->Finalize();    
   M->FormSystemMatrix(ess_tdof_list, SIPmat);
   SIPmat *= 0; 

   for (size_t i = 0; i < weights01.size(); ++i) {
      ProductCoefficient pro(dt, *weights01[i]); 
      SumCoefficient sum(1.0, *poles01[i], 1.0, dt); 
      QuotientCoefficient quo(pro, sum); 

      M = new ParBilinearForm(&fespace);
      M->AddDomainIntegrator(new MassIntegrator(quo));
      M->Assemble(0);
      M->Finalize(); 
      M->FormSystemMatrix(ess_tdof_list, TMPmat);

      BETAmat01.push_back(TMPmat); 
   }

   for (size_t i = 0; i < weights12.size(); ++i) {
      
      ProductCoefficient pro(dt, *weights12[i]); 
      SumCoefficient sum(1.0, *poles12[i], 1.0, dt); 
      QuotientCoefficient quo(pro, sum); 

      M = new ParBilinearForm(&fespace);
      M->AddDomainIntegrator(new MassIntegrator(quo));
      M->Assemble(0);
      M->Finalize(); 
      M->FormSystemMatrix(ess_tdof_list, TMPmat);   
      BETAmat12.push_back(TMPmat); 
   }

   for (size_t i = 0; i < weights23.size(); ++i) {      
      ProductCoefficient pro(dt, *weights23[i]); 
      SumCoefficient sum(1.0, *poles23[i], 1.0, dt); 
      QuotientCoefficient quo(pro, sum); 

      M = new ParBilinearForm(&fespace);
      M->AddDomainIntegrator(new MassIntegrator(quo));
      M->Assemble(0);
      M->Finalize(); 
      M->FormSystemMatrix(ess_tdof_list, TMPmat);   
      BETAmat23.push_back(TMPmat); 
   }

   for (size_t i = 0; i < poles01.size(); ++i) {
      SumCoefficient sum(1.0, *poles01[i], 1.0, dt); 
      QuotientCoefficient quo(con, sum); 

      M = new ParBilinearForm(&fespace);
      M->AddDomainIntegrator(new MassIntegrator(quo));
      M->Assemble(0);
      M->Finalize(); 
      M->FormSystemMatrix(ess_tdof_list, TMPmat);   
      GAMMAmat01.push_back(TMPmat); 
   }
   
   for (size_t i = 0; i < poles12.size(); ++i) {
      SumCoefficient sum(1.0, *poles12[i], 1.0, dt); 
      QuotientCoefficient quo(con, sum); 

      M = new ParBilinearForm(&fespace);
      M->AddDomainIntegrator(new MassIntegrator(quo));
      M->Assemble(0);
      M->Finalize(); 
      M->FormSystemMatrix(ess_tdof_list, TMPmat);   
      GAMMAmat12.push_back(TMPmat); 
   }

   for (size_t i = 0; i < poles23.size(); ++i) {
      SumCoefficient sum(1.0, *poles23[i], 1.0, dt); 
      QuotientCoefficient quo(con, sum); 

      M = new ParBilinearForm(&fespace);
      M->AddDomainIntegrator(new MassIntegrator(quo));
      M->Assemble(0);
      M->Finalize(); 
      M->FormSystemMatrix(ess_tdof_list, TMPmat);   
      GAMMAmat23.push_back(TMPmat); 
   }

   M = new ParBilinearForm(&fespace);
   std::vector<std::shared_ptr<ProductCoefficient>> pro01_coeffs, pro12_coeffs, pro23_coeffs; 
   std::vector<std::shared_ptr<SumCoefficient>> sum01_coeffs, sum12_coeffs, sum23_coeffs;
   std::vector<std::shared_ptr<QuotientCoefficient>> quo01_coeffs, quo12_coeffs, quo23_coeffs;
      
   for (size_t i = 0; i < weights01.size(); ++i) {
      auto pro = std::make_shared<ProductCoefficient>(dt * dt * dt, *weights01[i]);
      auto sum = std::make_shared<SumCoefficient>(1.0, *poles01[i], 1.0, dt);
      auto quo = std::make_shared<QuotientCoefficient>(*pro, *sum);
      pro01_coeffs.push_back(pro);
      sum01_coeffs.push_back(sum);
      quo01_coeffs.push_back(quo);
      M->AddDomainIntegrator(new MassIntegrator(*quo));
   }

   for (size_t i = 0; i < weights12.size(); ++i) {
      auto pro = std::make_shared<ProductCoefficient>(dt * dt, *weights12[i]);
      auto sum = std::make_shared<SumCoefficient>(1.0, *poles12[i], 1.0, dt);
      auto quo = std::make_shared<QuotientCoefficient>(*pro, *sum);
      pro01_coeffs.push_back(pro);
      sum01_coeffs.push_back(sum);
      quo01_coeffs.push_back(quo);
      M->AddDomainIntegrator(new MassIntegrator(*quo));
   }

   for (size_t i = 0; i < weights23.size(); ++i) {
      auto pro = std::make_shared<ProductCoefficient>(dt, *weights23[i]);
      auto sum = std::make_shared<SumCoefficient>(1.0, *poles23[i], 1.0, dt);
      auto quo = std::make_shared<QuotientCoefficient>(*pro, *sum);
      pro01_coeffs.push_back(pro);
      sum01_coeffs.push_back(sum);
      quo01_coeffs.push_back(quo);
      M->AddDomainIntegrator(new MassIntegrator(*quo));
   }

   M->Assemble(0);
   M->Finalize();
   M->FormSystemMatrix(ess_tdof_list, SIPmat);
}

DO03SpatialDiffusionWaveOperator::~DO03SpatialDiffusionWaveOperator()
{
delete T;
delete M;
delete K;
delete c2;
}