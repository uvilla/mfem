#include <mfem.hpp>
#include <fstream>
#include <iostream>
#include "pde_filter.hpp"


double TFunc(const mfem::Vector& a){
    double rez=sin(2*sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]));
    if(std::fabs(rez)>0.5){ rez=1.0;} else {rez=0.0;}
	return rez;
}




int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
 
   
   // 2. Parse command-line options.
   const char *mesh_file = "../data/star.mesh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(std::cout);
      }
      MPI_Finalize();
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(std::cout);
   }
    
   //generate parallel mesh
   mfem::Mesh *mesh = new mfem::Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   {
      int ref_levels =
         (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }
   
   mfem::ParMesh *pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   /*
   {
      int par_ref_levels = 2;
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh->UniformRefinement();
      }
   }
   */
   std::map<int,double> bndr;
   mtop::PDEFilter* filt=new mtop::PDEFilterFEMP(pmesh, bndr);
   
   
   //define the input and the output
   mfem::FiniteElementCollection *ifec;
   ifec = new mfem::L2_FECollection(0,dim, mfem::BasisType::Positive, mfem::FiniteElement::VALUE);
   
   int vdim=1;
   mfem::ParFiniteElementSpace* ifesp = new mfem::ParFiniteElementSpace(pmesh, ifec, vdim);
   mfem::ParFiniteElementSpace* ofesp = new mfem::ParFiniteElementSpace(pmesh, filt->FEC(), vdim);
   
   
    mfem::FunctionCoefficient ifun(TFunc);
    
    mfem::ParGridFunction* igf = new mfem::ParGridFunction(ifesp);
    mfem::ParGridFunction* ogf = new mfem::ParGridFunction(ofesp);
    
    
    
    igf->ProjectCoefficient(ifun);
    
    ogf->ProjectGridFunction(*igf);
    
    //pde_filter 
    filt->FFilter(igf,ogf); 
    
    //filt->FFilter(&ifun, ogf);
    
    mfem::ParaViewDataCollection* pvdc=new mfem::ParaViewDataCollection("outdc", pmesh);
    pvdc->RegisterField("output",ogf);   
    pvdc->RegisterField("input",igf);
    pvdc->SetLevelsOfDetail(2);
    pvdc->SetTime(0.0);
    pvdc->SetTimeStep(0);
    pvdc->Save();
    /*
    filt->FFilter(&ifun, ogf);
    pvdc->SetTime(1.0);
    pvdc->SetTimeStep(1);
    pvdc->Save();
    */
    
    delete pvdc; 
    
    delete ogf;
    delete igf;
    
   
   delete ifesp;
   delete ifec;
   delete ofesp;
   delete filt;
   delete pmesh;
   MPI_Finalize();

   return 0;
}
