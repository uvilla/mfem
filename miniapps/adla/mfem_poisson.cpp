//                                MFEM Example 1
//
// Compile with: make adla_poisson.exe
//
// Sample runs:  ./adla_poisson.exe -m ../../data/inline-segment.mesh
//		 ./adla_poisson.exe -m ../../data/inline-tri.mesh
// 		 ./adla_poisson.exe -m ../../data/inline-quad.mesh
//		 ./adla_poisson.exe -m ../../data/inline-tet.mesh
//               ./adla_poisson.exe -m ../../data/inline-hex.mesh
//               ./adla_poisson.exe -m ../../data/square-disc.mesh
//               ./adla_poisson.exe -m ../../data/star.mesh
//               ./adla_poisson.exe -m ../../data/escher.mesh
//               ./adla_poisson.exe -m ../../data/fichera.mesh
//               ./adla_poisson.exe -m ../../data/square-disc-p2.vtk -o 2
//               ./adla_poisson.exe -m ../../data/square-disc-p3.mesh -o 3
//               ./adla_poisson.exe -m ../../data/amr-quad.mesh
//               ./adla_poisson.exe -m ../../data/amr-hex.mesh
//               ./adla_poisson.exe -m ../../data/fichera-amr.mesh
//
//
// Description:  This example computes
//               int_\Omega e^m *( grad(u), grad(p) ) dx
//               using MFEM objects

#include <mfem.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
#include <map>

double ufun(const mfem::Vector & x);
double pfun(const mfem::Vector & x);
double mfun(const mfem::Vector & x);

double mfem_functional(mfem::Mesh * mesh, const mfem::IntegrationRule * ir,
					   mfem::GridFunction & u, mfem::GridFunction & p, mfem::GridFunction & m)
{
	const int dim = mesh->Dimension();
	   mfem::IsoparametricTransformation eltrans;
	   mfem::Vector grad_u_local(dim), grad_p_local(dim);

	   double integral(0.);

	   for (int el = 0; el < mesh -> GetNE(); el++)
	   {
		   mesh->GetElementTransformation(el, &eltrans);

		   for(int iq(0); iq < ir->GetNPoints(); ++iq)
		   {
			   const mfem::IntegrationPoint & ipt = ir->IntPoint(iq);
			   eltrans.SetIntPoint(&ipt);
			   u.GetGradient(eltrans, grad_u_local);
			   p.GetGradient(eltrans, grad_p_local);
			   integral += exp(m.GetValue(el, ipt))*mfem::InnerProduct(grad_u_local, grad_p_local)*ipt.weight*eltrans.Weight();
		   }
	   }

	   return integral;
}


int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   std::map<std::string, std::string> info;
   const char *mesh_file = "../../data/inline-tri.mesh";
   int order_u = 1;
   int order_m = 1;
   bool visualization = true;

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order_u, "-ou", "--order_u",
                  "Finite element order (polynomial degree) for u, p");
   args.AddOption(&order_m, "-om", "--order_m",
                  "Finite element order (polynomial degree) for m");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(std::cout);
      return 1;
   }
   args.PrintOptions(std::cout);


   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   mfem::Mesh *mesh = new mfem::Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   {
	  int target_nelem;
	  switch(order_u+order_m)
	  {
	  case 2:
		  target_nelem = 10000;
		  break;
	  case 3:
		  target_nelem = 5000;
		  break;
	  case 4:
		  target_nelem = 2500;
		  break;
	  default:
		  target_nelem = 1000;
	  }

	  target_nelem *= (4 - dim);

      int ref_levels =
         (int)floor(log(double(target_nelem)/mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 5. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order.
   mfem::FiniteElementCollection *ufec = new mfem::H1_FECollection(order_u, dim);
   mfem::FiniteElementCollection *pfec = new mfem::H1_FECollection(order_u, dim);
   mfem::FiniteElementCollection *mfec = new mfem::H1_FECollection(order_m, dim);

   mfem::FiniteElementSpace *ufespace = new mfem::FiniteElementSpace(mesh, ufec);
   mfem::FiniteElementSpace *pfespace = new mfem::FiniteElementSpace(mesh, pfec);
   mfem::FiniteElementSpace *mfespace = new mfem::FiniteElementSpace(mesh, mfec);

   std::cout << "Dimension = " << dim << std::endl;
   std::cout << "Number of elements = " <<  mesh->GetNE() << std::endl;
   std::cout << "Number of u DOFS = " << ufespace->GetTrueVSize() << std::endl;
   std::cout << "Number of m DOFS = " << mfespace->GetTrueVSize() << std::endl;

   info.emplace("Order u", std::to_string(order_u));
   info.emplace("Order m", std::to_string(order_m));
   info.emplace("Dimension", std::to_string(dim));
   info.emplace("Elements", std::to_string(mesh->GetNE()));
   info.emplace("u DOFS", std::to_string(ufespace->GetTrueVSize()));
   info.emplace("m DOFS", std::to_string(mfespace->GetTrueVSize()));

   // 8. Define the finite element variables u, p, m
   mfem::GridFunction u(ufespace), p(pfespace), m(mfespace);
   mfem::FunctionCoefficient uCoeff(ufun), pCoeff(pfun), mCoeff(mfun);
   u.ProjectCoefficient(uCoeff);
   p.ProjectCoefficient(pCoeff);
   m.ProjectCoefficient(mCoeff);


   // Assemble functional using mfem
   const int geo_type = mesh->GetElementBaseGeometry(0);
   int quad_order = 2*(order_u-1) + 2*order_m;
   if(geo_type == mfem::Geometry::CUBE || geo_type == mfem::Geometry::SQUARE)
	   quad_order += 2;
   const mfem::IntegrationRule * ir = &(mfem::IntRules.Get(geo_type, quad_order) );

   const int elem_udof = ufec->FiniteElementForGeometry(static_cast<mfem::Geometry::Type>(geo_type))->GetDof();
   const int elem_mdof = mfec->FiniteElementForGeometry(static_cast<mfem::Geometry::Type>(geo_type))->GetDof();

   std::cout << "Geometry type = " << geo_type << std::endl;
   std::cout << "Quadrature order = " << quad_order << std::endl;
   std::cout << "Number of quadrature points = " << ir->GetNPoints() <<std::endl;
   std::cout << "Number of local u DOFS = " << elem_udof << std::endl;
   std::cout << "Number of local m DOFS = " << elem_mdof << std::endl;

   info.emplace("Geometry type", std::to_string(geo_type));
   info.emplace("Quadrature order", std::to_string(quad_order));
   info.emplace("Number of quadrature points", std::to_string(ir->GetNPoints()));
   info.emplace("local u DOFS", std::to_string(elem_udof));
   info.emplace("local m DOFS", std::to_string(elem_mdof));

   mfem::tic();
   double j_functional = mfem_functional(mesh, ir, u, p, m);
   const double time_mfem_functional = mfem::toc();

   mfem::Vector dJ_dm(mfespace->GetVSize());
   dJ_dm = 0.;
   mfem::tic();
   {
   	   mfem::IsoparametricTransformation eltrans;
   	   mfem::Vector grad_u_local(dim), grad_p_local(dim), local_m_tilde;
   	   mfem::Array<int> m_vdofs;

   	   for (int el = 0; el < mesh -> GetNE(); el++)
   	   {
   		   mesh->GetElementTransformation(el, &eltrans);

   		   mfespace->GetElementVDofs(el, m_vdofs);
   		   local_m_tilde.SetSize(m_vdofs.Size());

   		   const mfem::FiniteElement * fe_m = mfec->FiniteElementForGeometry(mesh->GetElementBaseGeometry(el));

   		   for(int iq(0); iq < ir->GetNPoints(); ++iq)
   		   {
   			   const mfem::IntegrationPoint & ipt = ir->IntPoint(iq);
   			   eltrans.SetIntPoint(&ipt);
   			   u.GetGradient(eltrans, grad_u_local);
   			   p.GetGradient(eltrans, grad_p_local);
   			   const double s = exp(m.GetValue(el, ipt))*mfem::InnerProduct(grad_u_local, grad_p_local)*ipt.weight*eltrans.Weight();
   			   fe_m->CalcShape(ipt, local_m_tilde);
   			   dJ_dm.AddElementVector(m_vdofs, s, local_m_tilde);
   		   }
   	   }
      }
   const double time_mfem_dJ_dm = mfem::toc();

   std::cout << "J (MFEM) = " << j_functional << std::endl;


   //FD Check
   mfem::Vector tilde_m(m.Size());
   tilde_m.Randomize(100);
   double eps = 1e-8;

   mfem::GridFunction mplus(mfespace);
   mplus = m;
   mplus.Add(eps, tilde_m);

   mfem::GridFunction mminus(mfespace);
   mminus = m;
   mminus.Add(-eps, tilde_m);

   const double j_plus = mfem_functional(mesh, ir, u, p, mplus);
   const double j_minus = mfem_functional(mesh, ir, u, p, mminus);

   double d_integral(j_plus-j_minus);
   double fd_g_tildem(d_integral/(2.*eps));
   double g_tildem(mfem::InnerProduct(dJ_dm, tilde_m) );
   std::cout << "FD gradient = " << fd_g_tildem << std::endl;
   std::cout << "MFEM gradient = " << g_tildem <<std::endl;
   std::cout << "FD error check= " << std::abs(fd_g_tildem-g_tildem) << std::endl;

   //Second derivatives
   mfem::SparseMatrix H(pfespace->GetVSize(), mfespace->GetVSize());
   mfem::tic();
   {
	   mfem::DenseMatrix local_H, grad_p_test;
   	   mfem::IsoparametricTransformation eltrans;
   	   mfem::Vector grad_u_local(dim), local_m_hat, grad_u_grad_p;
   	   mfem::Array<int> m_vdofs, p_vdofs;

	   for (int el = 0; el < mesh -> GetNE(); el++)
	   {
   		   mesh->GetElementTransformation(el, &eltrans);

   		   mfespace->GetElementVDofs(el, m_vdofs);
   		   pfespace->GetElementVDofs(el, p_vdofs);

   		   const int num_m_vdofs = m_vdofs.Size();
   		   const int num_p_vdofs = p_vdofs.Size();

   		   const mfem::FiniteElement * fe_m = mfec->FiniteElementForGeometry(mesh->GetElementBaseGeometry(el));
   		   const mfem::FiniteElement * fe_p = pfec->FiniteElementForGeometry(mesh->GetElementBaseGeometry(el));

   		   local_H.SetSize(num_p_vdofs, num_m_vdofs);
   		   local_H = 0.;

   		   local_m_hat.SetSize(num_m_vdofs);
   		   grad_p_test.SetSize(num_p_vdofs, dim);
   		   grad_u_grad_p.SetSize(num_p_vdofs);
   		   grad_u_grad_p = 0.;

   		   for(int iq(0); iq < ir->GetNPoints(); ++iq)
   		   {
   			   const mfem::IntegrationPoint & ipt = ir->IntPoint(iq);
   			   eltrans.SetIntPoint(&ipt);
   			   const double s = exp(m.GetValue(el, ipt))*ipt.weight*eltrans.Weight();
   			   u.GetGradient(eltrans, grad_u_local);
   			   fe_m->CalcShape(ipt, local_m_hat);
   			   fe_p->CalcPhysDShape(eltrans, grad_p_test);
   			   grad_p_test.Mult(grad_u_local, grad_u_grad_p);

   			   for(int pdof(0); pdof < num_p_vdofs; ++pdof)
   				   for(int mdof(0); mdof < num_m_vdofs; ++mdof)
   					   local_H(pdof, mdof) += s*grad_u_grad_p(pdof)*local_m_hat(mdof);
   		   }

		   H.AddSubMatrix(p_vdofs, m_vdofs, local_H, 1);

	   }
	   H.Finalize();
   }
   const double time_mfem_hessian = mfem::toc();
   {
	   std::ofstream fid("H_mfem.mtx");
	   H.PrintMatlab(fid);
   }

   std::cout << "MFEM functional time = " << time_mfem_functional << std::endl;
   std::cout << "MFEM dJ_dm time = " << time_mfem_dJ_dm << std::endl;
   std::cout << "MFEM ddJ_dpdm time = " << time_mfem_hessian << std::endl;

   info.emplace("MFEM J time", std::to_string(time_mfem_functional) );
   info.emplace("MFEM dJ_dm time", std::to_string(time_mfem_dJ_dm) );
   info.emplace("MFEM ddJ_dpdm time", std::to_string(time_mfem_hessian) );

   {
	   std::ofstream fid("mfem_poisson.json");
	   fid << "{";
	   for (const auto &p : info)
		   fid << '"' << p.first << '"'<<": " << p.second << ",";

	   fid << '"' << "Mesh" << '"'<<": " << '"' << mesh_file << '"';
	   fid << "}";
   }

   delete ufespace;
   delete pfespace;
   delete mfespace;
   delete ufec;
   delete pfec;
   delete mfec;
   delete mesh;

   return 0;
}

double ufun(const mfem::Vector & x)
{
	return sin(x[0])*cos(x[1]);
}

double pfun(const mfem::Vector & x)
{
	return cos(x[0])*sin(x[1]);
}

double mfun(const mfem::Vector & x)
{
	return -x.Norml2();
}
