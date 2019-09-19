// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include "../general/forall.hpp"
#include "bilininteg.hpp"
#include "gridfunc.hpp"

using namespace std;
namespace mfem
{

static const IntegrationRule &DefaultGetRule(const FiniteElement &trial_fe,
                                             const FiniteElement &test_fe)
{
  int order;
  if (trial_fe.Space() == FunctionSpace::Pk)
    {
      order = trial_fe.GetOrder() + test_fe.GetOrder() - 2;
    }
  else
    {
      // order = 2*el.GetOrder() - 2;  // <-- this seems to work fine too
      order = trial_fe.GetOrder() + test_fe.GetOrder() + trial_fe.GetDim() - 1;
    }
  if (trial_fe.Space() == FunctionSpace::rQk)
    {
      return RefinedIntRules.Get(trial_fe.GetGeomType(), order);
    }
  return IntRules.Get(trial_fe.GetGeomType(), order);
}

void ConvectionIntegrator::AssemblePA(const FiniteElementSpace &fes)
{

  // Assuming the same element type
  Mesh *mesh = fes.GetMesh();
  if (mesh->GetNE() == 0) { return; }
  const FiniteElement &el = *fes.GetFE(0);
  ElementTransformation *T = mesh->GetElementTransformation(0);
  const IntegrationRule *ir = IntRule ? IntRule : &DefaultGetRule(el, el);
  dim = mesh->Dimension();
  ne = fes.GetMesh()->GetNE();
  nq = ir->GetNPoints();
  geom = mesh->GetGeometricFactors(*ir, GeometricFactors::COORDINATES |
                                   GeometricFactors::JACOBIANS);
  maps = &el.GetDofToQuad(*ir, DofToQuad::TENSOR);
  dofs1D = maps->ndof;
  quad1D = maps->nqpt;
  pa_data.SetSize(dim*ne*nq, Device::GetMemoryType());

  const int NE = ne;
  const int NQ = nq;
  auto w = ir->GetWeights().Read();
  auto J = Reshape(geom->J.Read(), NQ,2,2,NE);
  auto v = Reshape(pa_data.Write(), 2, NQ, NE);
  const double COEFF = -1.0;

  MFEM_FORALL(e, NE,
  {
    for(int q=0; q<NQ; ++q)
     {
        const double J11 = J(q,0,0,e);
        const double J21 = J(q,1,0,e);
        const double J12 = J(q,0,1,e);
        const double J22 = J(q,1,1,e);

        const double w_coeff = w[q] * COEFF;

        double cx = 1.0;
        double cy = 1.0;
        v(0,q,e) =    w_coeff*(cx * J22 - cy * J12);
        v(1,q,e) =  - w_coeff*(cx * J21 - cy * J11);
      }
  });

}

void ConvectionIntegrator::PAElementMatrix(const FiniteElementSpace &fes, 
                                           ElementTransformation &eTrans,
                                           DenseMatrix &mat,int i)
{

  Mesh *mesh = fes.GetMesh();
  const FiniteElement *el = fes.GetFE(0);
  const IntegrationRule *ir =  &DefaultGetRule(*el, *el);
  const DofToQuad *maps = &el->GetDofToQuad(*ir,DofToQuad::TENSOR);

  const int NE  = mesh->GetNE();
  const int DIM = mesh->Dimension();
  const int D1D = maps->ndof;
  const int Q1D = maps->nqpt;

  auto B1D = mfem::Reshape(maps->B.Read(), Q1D, D1D);
  auto Bt1D = mfem::Reshape(maps->Bt.Read(), D1D, Q1D);

  auto G1D = mfem::Reshape(maps->G.Read(),Q1D, D1D);
  auto Gt1D = mfem::Reshape(maps->Gt.Read(),D1D, Q1D);
    
  auto D = Reshape(pa_data.Read(), DIM, Q1D, Q1D, NE);

  const int max_dim = 10;
  auto Me = Reshape(mat.GetData(),D1D,D1D,D1D,D1D);

  /*
  for(int e=i; e<i+1; e++)
  {

    //Create B and G
    double B[max_dim][max_dim][max_dim];
    double BG[max_dim][max_dim][max_dim];
    double GB[max_dim][max_dim][max_dim];

    for(int j1=0; j1<D1D; ++j1){
      for(int i1=0; i1<D1D; ++i1){
        for(int k1=0; k1<Q1D; k1++){
          B[k1][i1][j1]  = B1D(k1,i1)*B1D(k1,j1);
          GB[k1][i1][j2] = G1D(k1,i1)*B1D(k1,j1);
          BG[k1][i1][j2] = B1D(k1,i1)*G1D(k1,j1);
        }
      }
    }

    
    double W[2][max_dim][max_dim][max_dim];    
    for(int j2=0; j2<D1D; ++j2) {
      for(int i2=0; i2<D1D; ++i2) {
        for(int k1=0;k1<Q1D;k1++) {      
          
          double dot(0.0);
          for(int k2=0; k2<Q2D; ++k2) {
            dot += B[k2][i2][j2]*D(0,k1,k2);            
          }
          W[0][k1][i2][j2];
        }
      }
    }

    for(int j2=0; j2<D1D; ++j2) {
      for(int j1=0; j1<D1D; ++j1) {
        
        for(int i2=0; i2<D1D; ++i2) {
          for(int i1=0; i1<D1D; ++i1) {

            double dot(0), dot2(0);
            for(int k1=0; k1<Q1D; ++k1) {        
              dot += BG[k1][i1][j1]*W[0][k1][i2][j2];
            }
            D(i1,i2,j1,j2) = dot;
          }
        }

      }
    }

  }//e loop
  */

}//assembly method


} //namespace mfem
