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

// PA Mass Integrator

// PA Mass Assemble kernel
void VectorMassIntegrator::Setup(const FiniteElementSpace &fes)
{
   // Assuming the same element type
   Mesh *mesh = fes.GetMesh();
   if (mesh->GetNE() == 0) { return; }
   const FiniteElement &el = *fes.GetFE(0);
   ElementTransformation *T = mesh->GetElementTransformation(0);
   const IntegrationRule *ir
      = IntRule ? IntRule : &MassIntegrator::GetRule(el, el, *T);
   dim = mesh->Dimension();
   ne = fes.GetMesh()->GetNE();
   nq = ir->GetNPoints();
   geom = mesh->GetGeometricFactors(*ir, GeometricFactors::COORDINATES |
                                    GeometricFactors::JACOBIANS);
   maps = &el.GetDofToQuad(*ir, DofToQuad::TENSOR);
   dofs1D = maps->ndof;
   quad1D = maps->nqpt;
   pa_data.SetSize(ne*nq, Device::GetMemoryType());
   double coeff = 1.0;
   if (Q)
   {
      ConstantCoefficient *cQ = dynamic_cast<ConstantCoefficient*>(Q);
      MFEM_VERIFY(cQ != NULL, "only ConstantCoefficient is supported!");
      coeff = cQ->constant;
   }
   if (dim==1) { MFEM_ABORT("Not supported yet... stay tuned!"); }
   if (dim==2)
   {
      const double constant = coeff;
      const int NE = ne;
      const int NQ = nq;
      auto w = ir->GetWeights().Read();
      auto J = Reshape(geom->J.Read(), NQ,2,2,NE);
      auto v = Reshape(pa_data.Write(), NQ, NE);
      MFEM_FORALL(e, NE,
      {
         for (int q = 0; q < NQ; ++q)
         {
            const double J11 = J(q,0,0,e);
            const double J12 = J(q,1,0,e);
            const double J21 = J(q,0,1,e);
            const double J22 = J(q,1,1,e);
            const double detJ = (J11*J22)-(J21*J12);
            v(q,e) =  w[q] * constant * detJ;
         }
      });
   }
   if (dim==3)
   {
      const double constant = coeff;
      const int NE = ne;
      const int NQ = nq;
      auto W = ir->GetWeights().Read();
      auto J = Reshape(geom->J.Read(), NQ,3,3,NE);
      auto v = Reshape(pa_data.Write(), NQ,NE);
      MFEM_FORALL(e, NE,
      {
         for (int q = 0; q < NQ; ++q)
         {
            const double J11 = J(q,0,0,e), J12 = J(q,0,1,e), J13 = J(q,0,2,e);
            const double J21 = J(q,1,0,e), J22 = J(q,1,1,e), J23 = J(q,1,2,e);
            const double J31 = J(q,2,0,e), J32 = J(q,2,1,e), J33 = J(q,2,2,e);
            const double detJ = J11 * (J22 * J33 - J32 * J23) -
            /* */               J21 * (J12 * J33 - J32 * J13) +
            /* */               J31 * (J12 * J23 - J22 * J13);
            v(q,e) = W[q] * constant * detJ;
         }
      });
   }
}

template<const int T_D1D = 0,
         const int T_Q1D = 0>
static void PAVectorMassApply2D(const int NE,
                                const Array<double> &_B,
                                const Array<double> &_Bt,
                                const Vector &_op,
                                const Vector &_x,
                                Vector &_y,
                                const int d1d = 0,
                                const int q1d = 0)
{
   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;
   constexpr int VDIM = 2;
   MFEM_VERIFY(D1D <= MAX_D1D, "");
   MFEM_VERIFY(Q1D <= MAX_Q1D, "");
   auto B = Reshape(_B.Read(), Q1D, D1D);
   auto Bt = Reshape(_Bt.Read(), D1D, Q1D);
   auto op = Reshape(_op.Read(), Q1D, Q1D, NE);
   auto x = Reshape(_x.Read(), D1D, D1D, VDIM, NE);
   auto y = Reshape(_y.ReadWrite(), D1D, D1D, VDIM, NE);
   MFEM_FORALL(e, NE,
   {
      const int D1D = T_D1D ? T_D1D : d1d; // nvcc workaround
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      // the following variables are evaluated at compile time
      constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
      constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
      double sol_xy[max_Q1D][max_Q1D];
      for (int c = 0; c < VDIM; ++c)
      {
         for (int qy = 0; qy < Q1D; ++qy)
         {
            for (int qx = 0; qx < Q1D; ++qx)
            {
               sol_xy[qy][qx] = 0.0;
            }
         }
         for (int dy = 0; dy < D1D; ++dy)
         {
            double sol_x[max_Q1D];
            for (int qy = 0; qy < Q1D; ++qy)
            {
               sol_x[qy] = 0.0;
            }
            for (int dx = 0; dx < D1D; ++dx)
            {
               const double s = x(dx,dy,c,e);
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  sol_x[qx] += B(qx,dx)* s;
               }
            }
            for (int qy = 0; qy < Q1D; ++qy)
            {
               const double d2q = B(qy,dy);
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  sol_xy[qy][qx] += d2q * sol_x[qx];
               }
            }
         }
         for (int qy = 0; qy < Q1D; ++qy)
         {
            for (int qx = 0; qx < Q1D; ++qx)
            {
               sol_xy[qy][qx] *= op(qx,qy,e);
            }
         }
         for (int qy = 0; qy < Q1D; ++qy)
         {
            double sol_x[max_D1D];
            for (int dx = 0; dx < D1D; ++dx)
            {
               sol_x[dx] = 0.0;
            }
            for (int qx = 0; qx < Q1D; ++qx)
            {
               const double s = sol_xy[qy][qx];
               for (int dx = 0; dx < D1D; ++dx)
               {
                  sol_x[dx] += Bt(dx,qx) * s;
               }
            }
            for (int dy = 0; dy < D1D; ++dy)
            {
               const double q2d = Bt(dy,qy);
               for (int dx = 0; dx < D1D; ++dx)
               {
                  y(dx,dy,c,e) += q2d * sol_x[dx];
               }
            }
         }
      }
   });
}

template<const int T_D1D = 0,
         const int T_Q1D = 0,
         const int T_NBZ = 0>
static void SmemPAVectorMassApply2D(const int NE,
                                    const Array<double> &b_,
                                    const Array<double> &bt_,
                                    const Vector &op_,
                                    const Vector &x_,
                                    Vector &y_,
                                    const int d1d = 0,
                                    const int q1d = 0)
{
   constexpr int VDIM = 2;
   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;
   constexpr int NBZ = T_NBZ ? T_NBZ : 1;
   constexpr int MQ1 = T_Q1D ? T_Q1D : MAX_Q1D;
   constexpr int MD1 = T_D1D ? T_D1D : MAX_D1D;
   MFEM_VERIFY(D1D <= MD1, "");
   MFEM_VERIFY(Q1D <= MQ1, "");
   auto b = Reshape(b_.Read(), Q1D, D1D);
   auto op = Reshape(op_.Read(), Q1D, Q1D, NE);
   auto x = Reshape(x_.Read(), D1D, D1D, VDIM, NE);
   auto y = Reshape(y_.ReadWrite(), D1D, D1D, VDIM, NE);
   MFEM_FORALL_2D(e, NE, Q1D, Q1D, NBZ,
   {
      const int tidz = MFEM_THREAD_ID(z);
      const int D1D = T_D1D ? T_D1D : d1d;
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      constexpr int NBZ = T_NBZ ? T_NBZ : 1;
      constexpr int MQ1 = T_Q1D ? T_Q1D : MAX_Q1D;
      constexpr int MD1 = T_D1D ? T_D1D : MAX_D1D;
      constexpr int MDQ = (MQ1 > MD1) ? MQ1 : MD1;
      MFEM_SHARED double BBt[MQ1*MD1];
      double (*B)[MD1] = (double (*)[MD1]) BBt;
      double (*Bt)[MQ1] = (double (*)[MQ1]) BBt;
      MFEM_SHARED double sm0[NBZ][MDQ*MDQ];
      MFEM_SHARED double sm1[NBZ][MDQ*MDQ];
      double (*X)[MD1] = (double (*)[MD1]) (sm0 + tidz);
      double (*DQ)[MQ1] = (double (*)[MQ1]) (sm1 + tidz);
      double (*QQ)[MQ1] = (double (*)[MQ1]) (sm0 + tidz);
      double (*QD)[MD1] = (double (*)[MD1]) (sm1 + tidz);
      for (int c = 0; c < VDIM; ++c)
      {
         MFEM_FOREACH_THREAD(dy,y,D1D)
         {
            MFEM_FOREACH_THREAD(dx,x,D1D)
            {
               X[dy][dx] = x(dx,dy,c,e);
            }
         }
         if (tidz == 0)
         {
            MFEM_FOREACH_THREAD(d,y,D1D)
            {
               MFEM_FOREACH_THREAD(q,x,Q1D)
               {
                  B[q][d] = b(q,d);
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(dy,y,D1D)
         {
            MFEM_FOREACH_THREAD(qx,x,Q1D)
            {
               double dq = 0.0;
               for (int dx = 0; dx < D1D; ++dx)
               {
                  dq += X[dy][dx] * B[qx][dx];
               }
               DQ[dy][qx] = dq;
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(qy,y,Q1D)
         {
            MFEM_FOREACH_THREAD(qx,x,Q1D)
            {
               double qq = 0.0;
               for (int dy = 0; dy < D1D; ++dy)
               {
                  qq += DQ[dy][qx] * B[qy][dy];
               }
               QQ[qy][qx] = qq * op(qx, qy, e);
            }
         }
         MFEM_SYNC_THREAD;
         if (tidz == 0)
         {
            MFEM_FOREACH_THREAD(d,y,D1D)
            {
               MFEM_FOREACH_THREAD(q,x,Q1D)
               {
                  Bt[d][q] = b(q,d);
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(qy,y,Q1D)
         {
            MFEM_FOREACH_THREAD(dx,x,D1D)
            {
               double dq = 0.0;
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  dq += QQ[qy][qx] * Bt[dx][qx];
               }
               QD[qy][dx] = dq;
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(dy,y,D1D)
         {
            MFEM_FOREACH_THREAD(dx,x,D1D)
            {
               double dd = 0.0;
               for (int qy = 0; qy < Q1D; ++qy)
               {
                  dd += (QD[qy][dx] * Bt[dy][qy]);
               }
               y(dx, dy, c, e) += dd;
            }
         }
      }
   });
}

template<const int T_D1D = 0,
         const int T_Q1D = 0>
static void PAVectorMassApply3D(const int NE,
                                const Array<double> &_B,
                                const Array<double> &_Bt,
                                const Vector &_op,
                                const Vector &_x,
                                Vector &_y,
                                const int d1d = 0,
                                const int q1d = 0)
{
   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;
   constexpr int VDIM = 3;
   MFEM_VERIFY(D1D <= MAX_D1D, "");
   MFEM_VERIFY(Q1D <= MAX_Q1D, "");
   auto B = Reshape(_B.Read(), Q1D, D1D);
   auto Bt = Reshape(_Bt.Read(), D1D, Q1D);
   auto op = Reshape(_op.Read(), Q1D, Q1D, Q1D, NE);
   auto x = Reshape(_x.Read(), D1D, D1D, D1D, VDIM, NE);
   auto y = Reshape(_y.ReadWrite(), D1D, D1D, D1D, VDIM, NE);
   MFEM_FORALL(e, NE,
   {
      const int D1D = T_D1D ? T_D1D : d1d;
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      constexpr int max_D1D = T_D1D ? T_D1D : MAX_D1D;
      constexpr int max_Q1D = T_Q1D ? T_Q1D : MAX_Q1D;
      double sol_xyz[max_Q1D][max_Q1D][max_Q1D];
      for (int c = 0; c < VDIM; ++ c)
      {
         for (int qz = 0; qz < Q1D; ++qz)
         {
            for (int qy = 0; qy < Q1D; ++qy)
            {
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  sol_xyz[qz][qy][qx] = 0.0;
               }
            }
         }
         for (int dz = 0; dz < D1D; ++dz)
         {
            double sol_xy[max_Q1D][max_Q1D];
            for (int qy = 0; qy < Q1D; ++qy)
            {
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  sol_xy[qy][qx] = 0.0;
               }
            }
            for (int dy = 0; dy < D1D; ++dy)
            {
               double sol_x[max_Q1D];
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  sol_x[qx] = 0;
               }
               for (int dx = 0; dx < D1D; ++dx)
               {
                  const double s = x(dx,dy,dz,c,e);
                  for (int qx = 0; qx < Q1D; ++qx)
                  {
                     sol_x[qx] += B(qx,dx) * s;
                  }
               }
               for (int qy = 0; qy < Q1D; ++qy)
               {
                  const double wy = B(qy,dy);
                  for (int qx = 0; qx < Q1D; ++qx)
                  {
                     sol_xy[qy][qx] += wy * sol_x[qx];
                  }
               }
            }
            for (int qz = 0; qz < Q1D; ++qz)
            {
               const double wz = B(qz,dz);
               for (int qy = 0; qy < Q1D; ++qy)
               {
                  for (int qx = 0; qx < Q1D; ++qx)
                  {
                     sol_xyz[qz][qy][qx] += wz * sol_xy[qy][qx];
                  }
               }
            }
         }
         for (int qz = 0; qz < Q1D; ++qz)
         {
            for (int qy = 0; qy < Q1D; ++qy)
            {
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  sol_xyz[qz][qy][qx] *= op(qx,qy,qz,e);
               }
            }
         }
         for (int qz = 0; qz < Q1D; ++qz)
         {
            double sol_xy[max_D1D][max_D1D];
            for (int dy = 0; dy < D1D; ++dy)
            {
               for (int dx = 0; dx < D1D; ++dx)
               {
                  sol_xy[dy][dx] = 0;
               }
            }
            for (int qy = 0; qy < Q1D; ++qy)
            {
               double sol_x[max_D1D];
               for (int dx = 0; dx < D1D; ++dx)
               {
                  sol_x[dx] = 0;
               }
               for (int qx = 0; qx < Q1D; ++qx)
               {
                  const double s = sol_xyz[qz][qy][qx];
                  for (int dx = 0; dx < D1D; ++dx)
                  {
                     sol_x[dx] += Bt(dx,qx) * s;
                  }
               }
               for (int dy = 0; dy < D1D; ++dy)
               {
                  const double wy = Bt(dy,qy);
                  for (int dx = 0; dx < D1D; ++dx)
                  {
                     sol_xy[dy][dx] += wy * sol_x[dx];
                  }
               }
            }
            for (int dz = 0; dz < D1D; ++dz)
            {
               const double wz = Bt(dz,qz);
               for (int dy = 0; dy < D1D; ++dy)
               {
                  for (int dx = 0; dx < D1D; ++dx)
                  {
                     y(dx,dy,dz,c,e) += wz * sol_xy[dy][dx];
                  }
               }
            }
         }
      }
   });
}

template<const int T_D1D = 0,
         const int T_Q1D = 0>
static void SmemPAVectorMassApply3D(const int NE,
                                    const Array<double> &b_,
                                    const Array<double> &bt_,
                                    const Vector &op_,
                                    const Vector &x_,
                                    Vector &y_,
                                    const int d1d = 0,
                                    const int q1d = 0)
{
   constexpr int VDIM = 3;
   const int D1D = T_D1D ? T_D1D : d1d;
   const int Q1D = T_Q1D ? T_Q1D : q1d;
   constexpr int M1Q = T_Q1D ? T_Q1D : MAX_Q1D;
   constexpr int M1D = T_D1D ? T_D1D : MAX_D1D;
   MFEM_VERIFY(D1D <= M1D, "");
   MFEM_VERIFY(Q1D <= M1Q, "");
   auto b = Reshape(b_.Read(), Q1D, D1D);
   auto op = Reshape(op_.Read(), Q1D, Q1D, Q1D, NE);
   auto x = Reshape(x_.Read(), D1D, D1D, D1D, VDIM, NE);
   auto y = Reshape(y_.ReadWrite(), D1D, D1D, D1D, VDIM, NE);
   MFEM_FORALL_3D(e, NE, Q1D, Q1D, Q1D,
   {
      const int tidz = MFEM_THREAD_ID(z);
      const int D1D = T_D1D ? T_D1D : d1d;
      const int Q1D = T_Q1D ? T_Q1D : q1d;
      constexpr int MQ1 = T_Q1D ? T_Q1D : MAX_Q1D;
      constexpr int MD1 = T_D1D ? T_D1D : MAX_D1D;
      constexpr int MDQ = (MQ1 > MD1) ? MQ1 : MD1;
      MFEM_SHARED double sDQ[MQ1*MD1];
      double (*B)[MD1] = (double (*)[MD1]) sDQ;
      double (*Bt)[MQ1] = (double (*)[MQ1]) sDQ;
      MFEM_SHARED double sm0[MDQ*MDQ*MDQ];
      MFEM_SHARED double sm1[MDQ*MDQ*MDQ];
      double (*X)[MD1][MD1]   = (double (*)[MD1][MD1]) sm0;
      double (*DDQ)[MD1][MQ1] = (double (*)[MD1][MQ1]) sm1;
      double (*DQQ)[MQ1][MQ1] = (double (*)[MQ1][MQ1]) sm0;
      double (*QQQ)[MQ1][MQ1] = (double (*)[MQ1][MQ1]) sm1;
      double (*QQD)[MQ1][MD1] = (double (*)[MQ1][MD1]) sm0;
      double (*QDD)[MD1][MD1] = (double (*)[MD1][MD1]) sm1;
      for (int c = 0; c < VDIM; ++ c)
      {
         MFEM_FOREACH_THREAD(dz,z,D1D)
         {
            MFEM_FOREACH_THREAD(dy,y,D1D)
            {
               MFEM_FOREACH_THREAD(dx,x,D1D)
               {
                  X[dz][dy][dx] = x(dx,dy,dz,c,e);
               }
            }
         }
         if (tidz == 0)
         {
            MFEM_FOREACH_THREAD(d,y,D1D)
            {
               MFEM_FOREACH_THREAD(q,x,Q1D)
               {
                  B[q][d] = b(q,d);
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(dz,z,D1D)
         {
            MFEM_FOREACH_THREAD(dy,y,D1D)
            {
               MFEM_FOREACH_THREAD(qx,x,Q1D)
               {
                  double u = 0.0;
                  for (int dx = 0; dx < D1D; ++dx)
                  {
                     u += X[dz][dy][dx] * B[qx][dx];
                  }
                  DDQ[dz][dy][qx] = u;
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(dz,z,D1D)
         {
            MFEM_FOREACH_THREAD(qy,y,Q1D)
            {
               MFEM_FOREACH_THREAD(qx,x,Q1D)
               {
                  double u = 0.0;
                  for (int dy = 0; dy < D1D; ++dy)
                  {
                     u += DDQ[dz][dy][qx] * B[qy][dy];
                  }
                  DQQ[dz][qy][qx] = u;
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(qz,z,Q1D)
         {
            MFEM_FOREACH_THREAD(qy,y,Q1D)
            {
               MFEM_FOREACH_THREAD(qx,x,Q1D)
               {
                  double u = 0.0;
                  for (int dz = 0; dz < D1D; ++dz)
                  {
                     u += DQQ[dz][qy][qx] * B[qz][dz];
                  }
                  QQQ[qz][qy][qx] = u * op(qx,qy,qz,e);
               }
            }
         }
         MFEM_SYNC_THREAD;
         if (tidz == 0)
         {
            MFEM_FOREACH_THREAD(d,y,D1D)
            {
               MFEM_FOREACH_THREAD(q,x,Q1D)
               {
                  Bt[d][q] = b(q,d);
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(qz,z,Q1D)
         {
            MFEM_FOREACH_THREAD(qy,y,Q1D)
            {
               MFEM_FOREACH_THREAD(dx,x,D1D)
               {
                  double u = 0.0;
                  for (int qx = 0; qx < Q1D; ++qx)
                  {
                     u += QQQ[qz][qy][qx] * Bt[dx][qx];
                  }
                  QQD[qz][qy][dx] = u;
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(qz,z,Q1D)
         {
            MFEM_FOREACH_THREAD(dy,y,D1D)
            {
               MFEM_FOREACH_THREAD(dx,x,D1D)
               {
                  double u = 0.0;
                  for (int qy = 0; qy < Q1D; ++qy)
                  {
                     u += QQD[qz][qy][dx] * Bt[dy][qy];
                  }
                  QDD[qz][dy][dx] = u;
               }
            }
         }
         MFEM_SYNC_THREAD;
         MFEM_FOREACH_THREAD(dz,z,D1D)
         {
            MFEM_FOREACH_THREAD(dy,y,D1D)
            {
               MFEM_FOREACH_THREAD(dx,x,D1D)
               {
                  double u = 0.0;
                  for (int qz = 0; qz < Q1D; ++qz)
                  {
                     u += QDD[qz][dy][dx] * Bt[dz][qz];
                  }
                  y(dx,dy,dz,c,e) += u;
               }
            }
         }
      }
   });
}

static void PAVectorMassApply(const int dim,
                              const int D1D,
                              const int Q1D,
                              const int NE,
                              const Array<double> &B,
                              const Array<double> &Bt,
                              const Vector &op,
                              const Vector &x,
                              Vector &y)
{
   if (dim == 2)
   {
      switch ((D1D << 4 ) | Q1D)
      {
         case 0x22: return SmemPAVectorMassApply2D<2,2,16>(NE, B, Bt, op, x, y);

         case 0x33: return SmemPAVectorMassApply2D<3,3,16>(NE, B, Bt, op, x, y);
         case 0x34: return SmemPAVectorMassApply2D<3,4,16>(NE, B, Bt, op, x, y);

         case 0x44: return SmemPAVectorMassApply2D<4,4,8>(NE, B, Bt, op, x, y);
         case 0x46: return SmemPAVectorMassApply2D<4,6,4>(NE, B, Bt, op, x, y);

         case 0x55: return SmemPAVectorMassApply2D<5,5,2>(NE, B, Bt, op, x, y);
         case 0x58: return SmemPAVectorMassApply2D<5,8,2>(NE, B, Bt, op, x, y);

         case 0x66: return SmemPAVectorMassApply2D<6,6,4>(NE, B, Bt, op, x, y);
         case 0x6A: return SmemPAVectorMassApply2D<6,10,2>(NE, B, Bt, op, x, y);

         case 0x77: return SmemPAVectorMassApply2D<7,7,4>(NE, B, Bt, op, x, y);
         case 0x7C: return SmemPAVectorMassApply2D<7,12,2>(NE, B, Bt, op, x, y);

         case 0x88: return SmemPAVectorMassApply2D<8,8,2>(NE, B, Bt, op, x, y);
         case 0x8E: return SmemPAVectorMassApply2D<8,14,1>(NE, B, Bt, op, x, y);

         case 0x99: return SmemPAVectorMassApply2D<9,9,2>(NE, B, Bt, op, x, y);
         case 0x90: return SmemPAVectorMassApply2D<9,16,1>(NE, B, Bt, op, x, y);
         default:   return PAVectorMassApply2D(NE, B, Bt, op, x, y, D1D, Q1D);
      }
   }
   if (dim == 3)
   {
      switch ((D1D << 4 ) | Q1D)
      {
         case 0x23: return SmemPAVectorMassApply3D<2,3>(NE, B, Bt, op, x, y);

         case 0x34: return SmemPAVectorMassApply3D<3,4>(NE, B, Bt, op, x, y);
         case 0x35: return SmemPAVectorMassApply3D<3,5>(NE, B, Bt, op, x, y);

         case 0x45: return SmemPAVectorMassApply3D<4,5>(NE, B, Bt, op, x, y);
         case 0x48: return SmemPAVectorMassApply3D<4,8>(NE, B, Bt, op, x, y);

         case 0x56: return SmemPAVectorMassApply3D<5,6>(NE, B, Bt, op, x, y);
         case 0x5A: return SmemPAVectorMassApply3D<5,10>(NE, B, Bt, op, x, y);

         case 0x67: return SmemPAVectorMassApply3D<6,7>(NE, B, Bt, op, x, y);
         // below, kernels use too much shared data @ compile or runtime
         //case 0x6D: return SmemPAVectorMassApply3D<6,13>(NE, B, Bt, op, x, y);

         case 0x78: return SmemPAVectorMassApply3D<7,8>(NE, B, Bt, op, x, y);
         //case 0x7F: return SmemPAVectorMassApply3D<7,15>(NE, B, Bt, op, x, y);

         case 0x89: return SmemPAVectorMassApply3D<8,9>(NE, B, Bt, op, x, y);
         // case 0x92: return SmemPAVectorMassApply3D<8,18>(NE, B, Bt, op, x, y);

         case 0x9A: return SmemPAVectorMassApply3D<9,10>(NE, B, Bt, op, x, y);
         //case 0x94: return SmemPAVectorMassApply3D<9,20>(NE, B, Bt, op, x, y);

         default:   return PAVectorMassApply3D(NE, B, Bt, op, x, y, D1D, Q1D);
      }
   }
   MFEM_ABORT("Unknown kernel.");
}

void VectorMassIntegrator::AddMultPA(const Vector &x, Vector &y) const
{
   PAVectorMassApply(dim, dofs1D, quad1D, ne, maps->B, maps->Bt, pa_data, x, y);
}

} // namespace mfem
