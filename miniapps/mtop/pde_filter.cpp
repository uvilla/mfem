#include "pde_filter.hpp"


namespace mtop{
    
void PDEFilterFEMS::Initialize(mfem::Mesh* mesh_, std::map<int, double>& bndr_)
{
    mfem_solver.device_config = "cpu";
    mfem_solver.order = 3;
    mfem_solver.pa = false;
    mfem_solver.static_cond = false; 

    fltr_conf.mesh = mesh_;
    fltr_conf.a = nullptr;
    
    fltr_conf.tt = nullptr;
    
    fltr_conf.vrhs = nullptr;
    fltr_conf.prec = nullptr;
    fltr_conf.fec = nullptr;
    fltr_conf.fespace = nullptr;
    fltr_conf.flag_assembl = false;
    
    //copy the boundary conditions
    bndr=bndr_; 
    
    InitCoeff();
}

void PDEFilterFEMS::Assemble()
{
    if(fltr_conf.flag_assembl){ //free the allocated mem
        fltr_conf.ess_tdof_list.clear();
        fltr_conf.all_ess_tdof_list.DeleteAll();
        
        delete fltr_conf.a;
        delete fltr_conf.tt;
        delete fltr_conf.vrhs;
        delete fltr_conf.fespace;
        delete fltr_conf.fec;
        
        if(fltr_conf.prec!=nullptr){ delete fltr_conf.prec;}
        fltr_conf.prec=nullptr;
    }
    
    //allocate mfem objects
    int dim = fltr_conf.mesh->Dimension();
    fltr_conf.fec = new mfem::H1_FECollection(mfem_solver.order,dim);
    //fltr_conf.fec = new mfem::L2_FECollection(mfem_solver.order,dim);
    fltr_conf.fespace = new mfem::FiniteElementSpace(fltr_conf.mesh, fltr_conf.fec);
    
    //allocate the bilinear form
    fltr_conf.a=new mfem::BilinearForm(fltr_conf.fespace);
    
    if (mfem_solver.pa) { 
        fltr_conf.a->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL); }
  
    if(fcoeff.mcoef_!=NULL){
        fltr_conf.a->AddDomainIntegrator(new mfem::DiffusionIntegrator(*fcoeff.mcoef_));}
    else if(fcoeff.scoef_!=NULL){
        fltr_conf.a->AddDomainIntegrator(new mfem::DiffusionIntegrator(*fcoeff.scoef_));}
    else{
        fltr_conf.a->AddDomainIntegrator(new mfem::DiffusionIntegrator(*fcoeff.one));}
    
    //add the mass matrix
    fltr_conf.a->AddDomainIntegrator(new mfem::MassIntegrator(*fcoeff.one));
    if (mfem_solver.static_cond) { fltr_conf.a->EnableStaticCondensation(); }
    
    //assemble the bilinear form
    fltr_conf.a->Assemble();
    
    //process the boundary conditions
    for(auto it = bndr.begin();it != bndr.end(); it++){
      if (fltr_conf.mesh->bdr_attributes.Size()){
           mfem::Array<int> ess_bdr(fltr_conf.mesh->bdr_attributes.Max());
           ess_bdr = false;
           ess_bdr[it->first] = true;
           fltr_conf.fespace->GetEssentialTrueDofs(ess_bdr, 
                                        fltr_conf.ess_tdof_list[it->first]);
           fltr_conf.all_ess_tdof_list.Append(fltr_conf.ess_tdof_list[it->first]);
      }
    }
    fltr_conf.all_ess_tdof_list.Sort();
    fltr_conf.all_ess_tdof_list.Unique();
    //end process the boundary conditions

    //allocate the RHS
    fltr_conf.vrhs = new mfem::GridFunction(fltr_conf.fespace); 
    
    *fltr_conf.vrhs=0;
    
    fltr_conf.flag_assembl=true;
}


PDEFilterFEMS::PDEFilterFEMS(mfem::Mesh* mesh_, std::map<int, double>& bnrd_)
{
    Initialize(mesh_,bnrd_);
    SetCoeff(1.0);
}

PDEFilterFEMS::~PDEFilterFEMS()
{
    fltr_conf.ess_tdof_list.clear();
    fltr_conf.all_ess_tdof_list.DeleteAll();
    
    if(fltr_conf.tt!=nullptr){
        delete fltr_conf.tt;
    }
        
        
    delete fltr_conf.a;
    delete fltr_conf.tt;
    delete fltr_conf.vrhs;
    delete fltr_conf.fespace;
    delete fltr_conf.fec;
        
    if(fltr_conf.prec!=nullptr){ delete fltr_conf.prec;}
    // free the coefficients 
    FreeCoeff();
}

//out is assumed to be in the FE space 
//obtained using FESpace() method
void PDEFilterFEMS::FFilter(mfem::FunctionCoefficient* in, mfem::GridFunction* out)
{
    if(fltr_conf.flag_assembl==false){ Assemble();}
    
    mfem::LinearForm *b = new mfem::LinearForm(fltr_conf.fespace);
    b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*in));
    b->Assemble();
    
    int copy_bc=1; //(1) coppy the BC otherwise (0) set them to zero
    
    //set the boundary conditions 
    for(auto it = bndr.begin();it != bndr.end(); it++){
        for(int ii = 0; ii < fltr_conf.ess_tdof_list[it->first].Size(); ii++){
            (*out)((fltr_conf.ess_tdof_list[it->first])[ii]) = it->second;
        }
    }
    
    fltr_conf.a->FormLinearSystem(fltr_conf.all_ess_tdof_list,
                                  *out,  *fltr_conf.vrhs,  
                                  fltr_conf.A, fltr_conf.X, fltr_conf.B, copy_bc);     
    
    if(fltr_conf.prec==nullptr){
        fltr_conf.prec = new mfem::GSSmoother((mfem::SparseMatrix&)(*fltr_conf.A));
    }
    
    mfem::PCG(*fltr_conf.A,*fltr_conf.prec, fltr_conf.B, fltr_conf.X, 1, 200, 1e-12, 0.0);
    
    fltr_conf.a->RecoverFEMSolution(fltr_conf.X, *fltr_conf.vrhs,*out);
    
    delete b;
}



void PDEFilterFEMS::FFilter(mfem::GridFunction* in, mfem::GridFunction* out)
{
    if(fltr_conf.flag_assembl==false){ Assemble();}
    
    if(fltr_conf.tt==nullptr){
        fltr_conf.tt=new mfem::MixedBilinearForm(in->FESpace(), fltr_conf.fespace); 
        fltr_conf.tt->AddDomainIntegrator(new mfem::MassIntegrator());
        fltr_conf.tt->Assemble();
    }
    
    fltr_conf.tt->Mult(*in, *fltr_conf.vrhs);
    
    int copy_bc=1; //(1) coppy the BC otherwise (0) set them to zero
    
    //set the boundary conditions 
    for(auto it = bndr.begin();it != bndr.end(); it++){
        for(int ii = 0; ii < fltr_conf.ess_tdof_list[it->first].Size(); ii++){
            (*out)((fltr_conf.ess_tdof_list[it->first])[ii]) = it->second;
        }
    }
    fltr_conf.a->FormLinearSystem(fltr_conf.all_ess_tdof_list,
                                  *out,  *fltr_conf.vrhs,  
                                  fltr_conf.A, fltr_conf.X, fltr_conf.B, copy_bc);     
    
    if(fltr_conf.prec==nullptr){
        fltr_conf.prec = new mfem::GSSmoother((mfem::SparseMatrix&)(*fltr_conf.A));
    }
    
    mfem::PCG(*fltr_conf.A,*fltr_conf.prec, fltr_conf.B, fltr_conf.X, 1, 200, 1e-12, 0.0);
    
    fltr_conf.a->RecoverFEMSolution(fltr_conf.X, *fltr_conf.vrhs,*out);
    
}

void PDEFilterFEMS::RFilter(mfem::GridFunction* in, mfem::GridFunction* out)
{
    if(!fltr_conf.flag_assembl){ Assemble();}
   
   
    //set the boundary conditions 
    for(auto it = bndr.begin();it != bndr.end(); it++){
        for(int ii = 0; ii < fltr_conf.ess_tdof_list[it->first].Size(); ii++){
            (*fltr_conf.vrhs)((fltr_conf.ess_tdof_list[it->first])[ii]) = 0.0;
        }
    }
    
    int copy_bc=1; //(1) coppy the BC otherwise (0) set them to zero
    fltr_conf.a->FormLinearSystem(fltr_conf.all_ess_tdof_list,
                                  *fltr_conf.vrhs, *in,
                                  fltr_conf.A, fltr_conf.X, fltr_conf.B, copy_bc);
    
    if(fltr_conf.prec==nullptr){
        fltr_conf.prec = new mfem::GSSmoother((mfem::SparseMatrix&)(*fltr_conf.A));
    }
    
    
    mfem::PCG(*fltr_conf.A,*fltr_conf.prec, fltr_conf.B, fltr_conf.X, 1, 200, 1e-12, 0.0);
    
    fltr_conf.a->RecoverFEMSolution(fltr_conf.X,*in, *fltr_conf.vrhs);
    //project the solution to the output GridFunction
    fltr_conf.tt->MultTranspose(*fltr_conf.vrhs,*out);
}

void PDEFilterFEMS::InitCoeff()
{
    fcoeff.flag_own = false;
    fcoeff.mcoef_ = nullptr;
    fcoeff.scoef_ = nullptr;
    fcoeff.one = new mfem::ConstantCoefficient(1.0);
}

void PDEFilterFEMS::FreeCoeff()
{
    if(fcoeff.flag_own){
      if(fcoeff.scoef_!=NULL){ delete fcoeff.scoef_;}
      if(fcoeff.mcoef_!=NULL){ delete fcoeff.mcoef_;}
    }
    delete fcoeff.one;
}

void PDEFilterFEMS::SetCoeff(double cf_)
{
     FreeCoeff();
     InitCoeff();
     fcoeff.scoef_=new mfem::ConstantCoefficient(cf_);
     fcoeff.flag_own=true;
     Assemble();
}

void PDEFilterFEMS::SetCoeff(mfem::DenseMatrix* mf_)
{
    FreeCoeff();
    InitCoeff();
    fcoeff.mcoef_=new  mfem::MatrixConstantCoefficient(*mf_);
    fcoeff.flag_own=true;
    Assemble();
}

void PDEFilterFEMS::SetCoeff(mfem::Coefficient* cf_)
{
    FreeCoeff();
    InitCoeff();
    fcoeff.scoef_=cf_;
    fcoeff.flag_own=false;
    Assemble();
}

void PDEFilterFEMS::SetCoeff(mfem::MatrixCoefficient* mf_)
{
    FreeCoeff();
    InitCoeff();
    fcoeff.mcoef_=mf_;
    fcoeff.flag_own=false;
    Assemble();
}


#ifdef MFEM_USE_MPI

void PDEFilterFEMP::Initialize(MPI_Comm comm_, mfem::Mesh* mesh_, std::map<int, double>& bndr_)
{
    mycomm=comm_;
    MPI_Comm_rank(mycomm,&myrank);
    MPI_Comm_size(mycomm,&nprocs);
    
    mfem_solver.order = 1;
    mfem_solver.static_cond = false;
    mfem_solver.device_config = "cpu";
    mfem_solver.pa = false;
    
    
    fltr_conf.aa = nullptr;
    fltr_conf.tt = nullptr;
    fltr_conf.vrhs = nullptr;
    fltr_conf.fec = nullptr;
    fltr_conf.fespace = nullptr;
    fltr_conf.prec = nullptr;
    fltr_conf.cg=nullptr;
    fltr_conf.flag_assembl = false;
    fltr_conf.mesh = dynamic_cast<mfem::ParMesh*>(mesh_);
    
    bndr=bndr_; //copy boundary conditions
    
    InitCoeff();
}

void PDEFilterFEMP::InitCoeff()
{
    fcoeff.flag_own = false;
    fcoeff.mcoef_ = nullptr;
    fcoeff.scoef_ = nullptr;
    fcoeff.one = new mfem::ConstantCoefficient(1.0);
}  

void PDEFilterFEMP::FreeCoeff()
{
    if(fcoeff.flag_own){
      if(fcoeff.scoef_!=NULL){ delete fcoeff.scoef_;}
      if(fcoeff.mcoef_!=NULL){ delete fcoeff.mcoef_;}
    }
    delete fcoeff.one;
}

void PDEFilterFEMP::SetCoeff(double c_)
{
     FreeCoeff();
     InitCoeff();
     fcoeff.scoef_ = new mfem::ConstantCoefficient(c_);
     fcoeff.flag_own = true;
     Assemble();
}

void PDEFilterFEMP::SetCoeff(mfem::DenseMatrix* cf_)
{
    FreeCoeff();
    InitCoeff();
    fcoeff.mcoef_ = new mfem::MatrixConstantCoefficient(*cf_);
    fcoeff.flag_own = true;
    Assemble();
}

void PDEFilterFEMP::SetCoeff(mfem::MatrixCoefficient* cd_)
{
    FreeCoeff();
    InitCoeff();
    fcoeff.mcoef_ = cd_;
    fcoeff.flag_own = false;
    Assemble();
}

void PDEFilterFEMP::SetCoeff(mfem::Coefficient* cf_)
{
    FreeCoeff();
    InitCoeff();
    fcoeff.scoef_ = cf_;
    fcoeff.flag_own = false;
    Assemble();
}

void PDEFilterFEMP::Assemble()
{
    
    if(fltr_conf.flag_assembl){ //free the allocated mem
        fltr_conf.ess_tdof_list.clear();
        fltr_conf.all_ess_tdof_list.DeleteAll();
        
        delete fltr_conf.tt;
        delete fltr_conf.vrhs;
        delete fltr_conf.fespace;
        delete fltr_conf.fec;
        
        fltr_conf.tt = nullptr;
        
        if(fltr_conf.prec!=nullptr){ delete fltr_conf.prec;}
        fltr_conf.prec = nullptr;
        
        if(fltr_conf.cg!=nullptr){ delete fltr_conf.cg;}
        fltr_conf.cg = nullptr;
        
        delete fltr_conf.aa;
            
    }
    
    //allocate mfem objects
    int dim = fltr_conf.mesh->Dimension();
    fltr_conf.fec = new mfem::H1_FECollection(mfem_solver.order,dim);
    fltr_conf.fespace = new mfem::ParFiniteElementSpace(fltr_conf.mesh, fltr_conf.fec);
    
    //allocate the bilinear form
    fltr_conf.aa=new mfem::ParBilinearForm(fltr_conf.fespace);

     if (mfem_solver.pa) { 
        fltr_conf.aa->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL); }
  
    if(fcoeff.mcoef_!=NULL){
        fltr_conf.aa->AddDomainIntegrator(new mfem::DiffusionIntegrator(*fcoeff.mcoef_));}
    else if(fcoeff.scoef_!=NULL){
        fltr_conf.aa->AddDomainIntegrator(new mfem::DiffusionIntegrator(*fcoeff.scoef_));}
    else{
        fltr_conf.aa->AddDomainIntegrator(new mfem::DiffusionIntegrator(*fcoeff.one));}
    
    //add the mass matrix
    fltr_conf.aa->AddDomainIntegrator(new mfem::MassIntegrator(*fcoeff.one));
    if (mfem_solver.static_cond) { fltr_conf.aa->EnableStaticCondensation(); }
    
    //assemble the bilinear form
    fltr_conf.aa->Assemble();

    //process the boundary conditions
    for(auto it = bndr.begin();it != bndr.end(); it++){
      if (fltr_conf.mesh->bdr_attributes.Size()){
           mfem::Array<int> ess_bdr(fltr_conf.mesh->bdr_attributes.Max());
           ess_bdr = false;
           ess_bdr[it->first] = true;
           fltr_conf.fespace->GetEssentialTrueDofs(ess_bdr, 
                                        fltr_conf.ess_tdof_list[it->first]);
           fltr_conf.all_ess_tdof_list.Append(fltr_conf.ess_tdof_list[it->first]);
      }
    }
    fltr_conf.all_ess_tdof_list.Sort();
    fltr_conf.all_ess_tdof_list.Unique();
    //end process the boundary conditions

    //allocate the RHS
    fltr_conf.vrhs = new mfem::ParGridFunction(fltr_conf.fespace); 
    
    *fltr_conf.vrhs=0;
    
    fltr_conf.flag_assembl=true;
}

PDEFilterFEMP::~PDEFilterFEMP()
{
    
    fltr_conf.ess_tdof_list.clear();
    fltr_conf.all_ess_tdof_list.DeleteAll();
    
    if(fltr_conf.tt!=nullptr){ delete fltr_conf.tt;}
    if(fltr_conf.cg!=nullptr){ delete fltr_conf.cg;}
    if(fltr_conf.prec!=nullptr){ delete fltr_conf.prec;}
    
    delete fltr_conf.aa;
    delete fltr_conf.vrhs;
    delete fltr_conf.fespace;
    delete fltr_conf.fec;
    
    // free the coefficients 
    FreeCoeff();
}

void PDEFilterFEMP::FFilter(mfem::FunctionCoefficient* in, mfem::GridFunction* out_)
{
    mfem::ParGridFunction* out = dynamic_cast<mfem::ParGridFunction*>(out_);
        
    mfem::ParLinearForm* b = new mfem::ParLinearForm(fltr_conf.fespace);
    b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*in));
    b->Assemble();
    
    int copy_bc=1; //(1) coppy the BC otherwise (0) set them to zero
        
    //set the boundary conditions 
    for(auto it = bndr.begin();it != bndr.end(); it++){
        for(int ii = 0; ii < fltr_conf.ess_tdof_list[it->first].Size(); ii++){
            (*out)((fltr_conf.ess_tdof_list[it->first])[ii]) = it->second;
        }
    }

    fltr_conf.aa->FormLinearSystem(fltr_conf.all_ess_tdof_list,
                                  *out,  *b,  
                                  fltr_conf.A, fltr_conf.X, fltr_conf.B, copy_bc); 

   if (!mfem_solver.pa) { fltr_conf.prec = new mfem::HypreBoomerAMG(); }
   if(fltr_conf.cg==nullptr)
   { 
       fltr_conf.cg=new mfem::CGSolver(mycomm);
       fltr_conf.cg->SetRelTol(1e-12);
       fltr_conf.cg->SetMaxIter(200);
       fltr_conf.cg->SetPrintLevel(1);
       if (fltr_conf.prec) { fltr_conf.cg->SetPreconditioner(*fltr_conf.prec); }
       fltr_conf.cg->SetOperator(*fltr_conf.A);
   }
   fltr_conf.cg->Mult(fltr_conf.B, fltr_conf.X);
   fltr_conf.aa->RecoverFEMSolution(fltr_conf.X, *b,*out);
   delete b;
}


void PDEFilterFEMP::FFilter(mfem::GridFunction* in_, mfem::GridFunction* out_)
{
    mfem::ParGridFunction* in = dynamic_cast<mfem::ParGridFunction*>(in_);
    mfem::ParGridFunction* out = dynamic_cast<mfem::ParGridFunction*>(out_);
    
    if(fltr_conf.flag_assembl==false){ Assemble();}
    
    if(fltr_conf.tt==nullptr){
        fltr_conf.tt = new mfem::ParMixedBilinearForm(in->ParFESpace(),fltr_conf.fespace); 
        fltr_conf.tt->AddDomainIntegrator(new mfem::MassIntegrator());
        fltr_conf.tt->Assemble();
    }
    
    fltr_conf.tt->Mult(*in,*fltr_conf.vrhs);
    
    int copy_bc=1; //(1) coppy the BC otherwise (0) set them to zero
    
    //set the boundary conditions 
    for(auto it = bndr.begin();it != bndr.end(); it++){
        for(int ii = 0; ii < fltr_conf.ess_tdof_list[it->first].Size(); ii++){
            (*out)((fltr_conf.ess_tdof_list[it->first])[ii]) = it->second;
        }
    }
    
    fltr_conf.aa->FormLinearSystem(fltr_conf.all_ess_tdof_list,
                                  *out,  *fltr_conf.vrhs,  
                                  fltr_conf.A, fltr_conf.X, fltr_conf.B, copy_bc);     
    
    if(fltr_conf.prec==nullptr){
        fltr_conf.prec = new mfem::HypreBoomerAMG(); }
        
    if(fltr_conf.cg==nullptr)
    { 
       fltr_conf.cg=new mfem::CGSolver(mycomm);
       fltr_conf.cg->SetRelTol(1e-12);
       fltr_conf.cg->SetMaxIter(200);
       fltr_conf.cg->SetPrintLevel(1);
       if (fltr_conf.prec) { fltr_conf.cg->SetPreconditioner(*fltr_conf.prec); }
       fltr_conf.cg->SetOperator(*fltr_conf.A);
    }
    fltr_conf.cg->Mult(fltr_conf.B, fltr_conf.X);
    fltr_conf.aa->RecoverFEMSolution(fltr_conf.X, *fltr_conf.vrhs ,*out);
}

void PDEFilterFEMP::RFilter(mfem::GridFunction* in_, mfem::GridFunction* out_)
{
    mfem::ParGridFunction* in = dynamic_cast<mfem::ParGridFunction*>(in_);
    mfem::ParGridFunction* out = dynamic_cast<mfem::ParGridFunction*>(out_);
    
    if(fltr_conf.flag_assembl==false){ Assemble();}
    
    int copy_bc = 1; //(1) coppy the BC otherwise (0) set them to zero
    *fltr_conf.vrhs = 0.0;
    
    fltr_conf.aa->FormLinearSystem(fltr_conf.all_ess_tdof_list,
                                  *fltr_conf.vrhs, *in, 
                                  fltr_conf.A, fltr_conf.X, fltr_conf.B, copy_bc);
    if(fltr_conf.prec==nullptr){
        fltr_conf.prec = new mfem::HypreBoomerAMG(); }
        
    if(fltr_conf.cg==nullptr)
    { 
       fltr_conf.cg=new mfem::CGSolver(mycomm);
       fltr_conf.cg->SetRelTol(1e-12);
       fltr_conf.cg->SetMaxIter(200);
       fltr_conf.cg->SetPrintLevel(1);
       if (fltr_conf.prec) { fltr_conf.cg->SetPreconditioner(*fltr_conf.prec); }
       fltr_conf.cg->SetOperator(*fltr_conf.A);
    }
    fltr_conf.cg->Mult(fltr_conf.B, fltr_conf.X);
    fltr_conf.aa->RecoverFEMSolution(fltr_conf.X, *in, *fltr_conf.vrhs);
        
    if(fltr_conf.tt==nullptr){
        fltr_conf.tt = new mfem::ParMixedBilinearForm(in->ParFESpace(),fltr_conf.fespace); 
        fltr_conf.tt->AddDomainIntegrator(new mfem::MassIntegrator());
        fltr_conf.tt->Assemble();
    }
    
    fltr_conf.tt->MultTranspose(*fltr_conf.vrhs,*out);
    
}


PDEFilterFEMP::PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int, double>& bndr_)
{
    Initialize(mesh->GetComm(), mesh, bndr_);
    SetCoeff(0.0);
}

PDEFilterFEMP::PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int, double>& bndr_, double cf_)
{
    Initialize(mesh->GetComm(), mesh, bndr_);
    SetCoeff(cf_);
}

PDEFilterFEMP::PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int, double>& bndr_, mfem::Coefficient* cf_)
{
    Initialize(mesh->GetComm(), mesh, bndr_);
    SetCoeff(cf_);
}

PDEFilterFEMP::PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int, double>& bndr_, mfem::DenseMatrix* cf_)
{
    Initialize(mesh->GetComm(), mesh, bndr_);
    SetCoeff(cf_);
}

PDEFilterFEMP::PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int, double>& bndr_, mfem::MatrixCoefficient* cd_)
{
    Initialize(mesh->GetComm(), mesh, bndr_);
    SetCoeff(cd_);
}


#endif //end MFEM_USE_MPI

}

