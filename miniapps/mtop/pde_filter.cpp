#include "pde_filter.hpp"


namespace mtop{
    
void PDEFilterFEMS::Initialize(mfem::Mesh* mesh_, std::map<int, double>& bndr_)
{
    mfem_solver.device_config = "cpu";
    mfem_solver.order = 1;
    mfem_solver.pa = false;
    mfem_solver.static_cond = false; 
    mfem_solver.visualization = true;

    
    fltr_conf.mesh = nullptr;
    fltr_conf.a = nullptr;
    fltr_conf.b = nullptr;
    
    fltr_conf.prec = nullptr;
    
    
}

PDEFilterFEMS::PDEFilterFEMS(mfem::Mesh* mesh_, std::map<int, double>& bnrd_)
{

    
}

PDEFilterFEMS::~PDEFilterFEMS()
{
    // free the coefficients 
    FreeCoeff();
}

void PDEFilterFEMS::FFilter(mfem::GridFunction& in, mfem::GridFunction& out)
{
    
}

void PDEFilterFEMS::RFilter(mfem::GridFunction& in, mfem::GridFunction& out)
{

    
}

void PDEFilterFEMS::InitCoeff()
{
    fcoeff.flag_own=false;
    fcoeff.mcoef_=nullptr;
    fcoeff.scoef_=nullptr;
}

void PDEFilterFEMS::FreeCoeff()
{
    if(fcoeff.flag_own){
      if(fcoeff.scoef_!=NULL){ delete fcoeff.scoef_;}
      if(fcoeff.mcoef_!=NULL){ delete fcoeff.mcoef_;}
    }
    InitCoeff();
}

void PDEFilterFEMS::SetCoeff(double cf_)
{
     FreeCoeff(); 
     fcoeff.scoef_=new mfem::ConstantCoefficient(cf_);
     fcoeff.flag_own=true;
}

void PDEFilterFEMS::SetCoeff(mfem::DenseMatrix& mf_)
{
    FreeCoeff();
    fcoeff.mcoef_=new  mfem::MatrixConstantCoefficient(mf_);
    fcoeff.flag_own=true;
}

void PDEFilterFEMS::SetCoeff(mfem::Coefficient& cf_)
{
    FreeCoeff();
    fcoeff.scoef_=&cf_;
    fcoeff.flag_own=false;
}

void PDEFilterFEMS::SetCoeff(mfem::MatrixCoefficient& mf_)
{
    FreeCoeff();
    fcoeff.mcoef_=&mf_;
    fcoeff.flag_own=false;
}

    
} // end namespace mtop
