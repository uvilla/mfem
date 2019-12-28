#ifndef PDE_FILTER_H
#define PDE_FILTER_H


#include <mfem.hpp>
#include <map>

// Standrard PDE filter implementation 

namespace mtop{

// base class
class Filter{

private:

public:
   Filter(){}
   virtual ~Filter(){}
   
   virtual void FFilter(mfem::GridFunction* in, mfem::GridFunction* out)=0;
   virtual void RFilter(mfem::GridFunction* in, mfem::GridFunction* out)=0;
   virtual void FFilter(mfem::FunctionCoefficient* in,  mfem::GridFunction* out)=0;
   
};


class PDEFilter: public Filter{
    
public:
    PDEFilter(){};
    virtual ~PDEFilter(){};
    
    virtual const mfem::FiniteElementCollection* FEC()=0;
    virtual const mfem::FiniteElementSpace* FESpace()=0;
    
    virtual void SetCoeff(double cf_)=0;
    virtual void SetCoeff(mfem::DenseMatrix* cf_)=0;
    virtual void SetCoeff(mfem::Coefficient* cf_)=0;
    virtual void SetCoeff(mfem::MatrixCoefficient* cd_)=0;
    
};



// Standard finite element discretization
// serial implementation
class PDEFilterFEMS: public PDEFilter{

    
public:
    PDEFilterFEMS(mfem::Mesh* mesh_, std::map<int,double>& bndr_);
 
    virtual const mfem::FiniteElementCollection* FEC() override {return fltr_conf.fec;}
    virtual const mfem::FiniteElementSpace* FESpace() override {return fltr_conf.fespace;} 

    virtual void SetCoeff(double cf_) override;  // set contant coefficient
    virtual void SetCoeff(mfem::DenseMatrix* mf_) override; //set constant tensor
    virtual void SetCoeff(mfem::Coefficient* cf_) override;
    virtual void SetCoeff(mfem::MatrixCoefficient* mf_) override;
    
    
    virtual ~PDEFilterFEMS() override;
    
    // forward filter
    // in - original density parametrizaiton
    // out - filtered filed H1 funtion
    virtual void FFilter(mfem::GridFunction* in, mfem::GridFunction* out) override;
    virtual void FFilter(mfem::FunctionCoefficient* in,  mfem::GridFunction* out) override;
    // reverse gradients 
    // in - gradients with respect to the filtered field
    // out - gradient with respect to the orifinal density parametrizaiton
    virtual void RFilter(mfem::GridFunction* in, mfem::GridFunction* out) override;

private:
    
    void Initialize(mfem::Mesh* mesh_, std::map<int,double>& bndr_);
    void Assemble(); 
    
    //boundary conditions
    std::map<int,double> bndr;
    

    //mfem objects
    struct{
     int order = 1;
     bool static_cond = false;
     bool pa = false;
     const char *device_config = "cpu";
    } mfem_solver;

    //configuration values
    struct{
        mfem::Mesh *mesh;// mesh for the optimization
        mfem::FiniteElementCollection *fec; // H1 finite element finite element collection
        mfem::FiniteElementSpace *fespace;
        
        mfem::BilinearForm *a;
        mfem::MixedBilinearForm *tt;
        
        
        std::map<int, mfem::Array<int>> ess_tdof_list;
        mfem::Array<int> all_ess_tdof_list;
    
        mfem::OperatorPtr A;
        mfem::Vector B;
        mfem::Vector X;
        
        mfem::GridFunction *vrhs; //RHS
        
        mfem::GSSmoother *prec; // GS smoother
        
        bool flag_assembl; // true of the bilinear form is assembled
                           // false otherwise
    } fltr_conf;
    
    struct{
        bool flag_own; //false if the Coefficients are owned by the class
        mfem::Coefficient *scoef_; //scalar Coefficient
        mfem::MatrixCoefficient *mcoef_; //matrix Coefficient 
        mfem::ConstantCoefficient *one; //constant coefficient one
    } fcoeff; 
    
    void InitCoeff();
    void FreeCoeff();
    
};


#ifdef MFEM_USE_MPI

class PDEFilterFEMP: public PDEFilter
{
public:
    
    PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int,double>& bndr_);
    PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int,double>& bndr_, double cf_);
    PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int,double>& bndr_, mfem::Coefficient* cf_);
    PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int,double>& bndr_, mfem::DenseMatrix* cf_);
    PDEFilterFEMP(mfem::ParMesh* mesh, std::map<int,double>& bndr_, mfem::MatrixCoefficient* cd_);
    
    virtual ~PDEFilterFEMP() override;
    
    virtual void FFilter(mfem::FunctionCoefficient* in, mfem::GridFunction* out) override;
    virtual void FFilter(mfem::GridFunction* in, mfem::GridFunction* out) override;
    virtual void RFilter(mfem::GridFunction* in, mfem::GridFunction* out) override;
    
    virtual void SetCoeff(double c_) override;
    virtual void SetCoeff(mfem::Coefficient* cf_) override;
    virtual void SetCoeff(mfem::DenseMatrix* cf_) override;
    virtual void SetCoeff(mfem::MatrixCoefficient* cd_) override;
    
    virtual const mfem::FiniteElementCollection* FEC() override {return fltr_conf.fec;}
    virtual const mfem::FiniteElementSpace* FESpace() override {return fltr_conf.fespace;}
    
private:
    
    void Initialize(MPI_Comm comm_, mfem::Mesh* mesh_, std::map<int,double>& bndr_);
    void Assemble(); 
    
    //boundary conditions
    std::map<int,double> bndr;
    

    //mfem objects
    struct{
     int order = 1;
     bool static_cond = false;
     bool pa = false;
     const char *device_config = "cpu";
    } mfem_solver;

    int myrank;
    int nprocs;
    MPI_Comm mycomm;
    
    struct{
      
        mfem::ParMesh *mesh;
        mfem::FiniteElementCollection* fec;
        mfem::ParFiniteElementSpace* fespace;
        
        mfem::ParBilinearForm* aa;
        mfem::ParMixedBilinearForm* tt;
        
        std::map<int, mfem::Array<int>> ess_tdof_list;
        mfem::Array<int> all_ess_tdof_list;
    
        mfem::OperatorPtr A;
        mfem::Vector B;
        mfem::Vector X;
        
        mfem::ParGridFunction *vrhs; //RHS
        mfem::Solver *prec;
        mfem::CGSolver *cg;
    
        bool flag_assembl; // true of the bilinear form is assembled
                           // false otherwise
    } fltr_conf;
    
    
    struct{
        bool flag_own; //false if the Coefficients are owned by the class
        mfem::Coefficient *scoef_; //scalar Coefficient
        mfem::MatrixCoefficient *mcoef_; //matrix Coefficient 
        mfem::ConstantCoefficient *one; //constant coefficient one
    } fcoeff;
    
    void InitCoeff();
    void FreeCoeff();
    
    
};
#endif
    
    
}

#endif
