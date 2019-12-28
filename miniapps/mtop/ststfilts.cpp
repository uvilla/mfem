
#include <mfem.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
#include "pde_filter.hpp"

double TFunc(const mfem::Vector& a){
    double rez=sin(sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]));
    if(std::fabs(rez)>0.5){ rez=1.0;} else {rez=0.0;}
	return rez;
}



void VFunc4(const mfem::Vector& x, mfem::Vector& out){
  out[0]=std::sin(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
  out[1]=std::cos(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
  out[2]=std::sin(2*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
  out[3]=std::cos(2*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));

  if(fabs(out[0])>0.5){ out[0]=1.0;} else {out[0]=0.0;}
  if(fabs(out[1])>0.5){ out[1]=1.0;} else {out[1]=0.0;}
  if(fabs(out[2])>0.5){ out[2]=1.0;} else {out[2]=0.0;}
  if(fabs(out[3])>0.5){ out[3]=1.0;} else {out[3]=0.0;}
}



int main(int argc, char *argv[])
{
    // 1. Parse command-line options.
    const char *mesh_file = "../../data/beam-hex.mesh";
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
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);

    //read the mesh 
    mfem::Mesh *mesh = new mfem::Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();
    for(int ii = 0; ii < 1; ii++){
        mesh->UniformRefinement();}
    
    
    std::map<int,double> bndr;
    mtop::PDEFilter* filt=new mtop::PDEFilterFEMS(mesh, bndr);
    filt->SetCoeff(0.000*0.00);
    
    //define the input and the output
    mfem::FiniteElementCollection *ifec;
    ifec = new mfem::L2_FECollection(0,dim, mfem::BasisType::Positive, mfem::FiniteElement::VALUE);
    
    //ifec = new mfem::H1_FECollection(order, dim, mfem::BasisType::Positive);
    int vdim=1;
    mfem::FiniteElementSpace* ifesp = new mfem::FiniteElementSpace(mesh, ifec, vdim);
    mfem::FunctionCoefficient ifun(TFunc);
    
    mfem::FiniteElementCollection *ofec;
    ofec = new mfem::H1_FECollection(order, dim,
                                     mfem::BasisType::Positive);
    mfem::FiniteElementSpace* ofesp = new mfem::FiniteElementSpace(mesh, filt->FEC(), vdim);
    

    mfem::GridFunction* igf = new mfem::GridFunction(ifesp);
    mfem::GridFunction* ogf = new mfem::GridFunction(ofesp);
    //initialize the input GridFunction
    igf->ProjectCoefficient(ifun);
    
    ogf->ProjectGridFunction(*igf);
    
    //pde_filter 
    filt->FFilter(igf,ogf); 
    
    mfem::ParaViewDataCollection* pvdc=new mfem::ParaViewDataCollection("outdc",mesh);
    pvdc->RegisterField("output",ogf);   
    pvdc->RegisterField("input",igf);
    pvdc->SetLevelsOfDetail(2);
    pvdc->SetTime(0.0);
    pvdc->SetTimeStep(0);
    pvdc->Save();
    
    filt->FFilter(&ifun, ogf);
    pvdc->SetTime(1.0);
    pvdc->SetTimeStep(1);
    pvdc->Save();
    
    
    
    delete pvdc;    
    
    delete igf;
    delete ogf;
    
    
    delete ofesp;
    delete ifesp;
    
    delete ofec;
    delete ifec;
    
    delete filt;
    delete mesh;
    
    return 0;
}
