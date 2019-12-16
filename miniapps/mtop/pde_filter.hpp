#ifndef PDE_FILTER_H
#define PDE_FILTER_H


#include <mfem.hpp>
#include <map>

// Standrard PDE filter implementation 

namespace mtop{

// base class
class PDEFilter{

private:

public:
   PDEFilter(){}
   virtual ~PDEFilter(){}
   
   virtual void FFilter(mfem::GridFunction& in, mfem::GridFunction& out)=0;
   virtual void RFilter(mfem::GridFunction& in, mfem::GridFunction& out)=0;
   
};


// Standard finite element discretization
// serial implementation
class PDEFilterFEMS: public PDEFilter{

    
public:
    PDEFilterFEMS(mfem::Mesh* mesh_, std::map<int,double>& bnrd_);

    virtual ~PDEFilterFEMS() override;
    
    // forward filter
    // in - original density parametrizaiton
    // out - filtered filed H1 funtion
    virtual void FFilter(mfem::GridFunction& in, mfem::GridFunction& out) override;
    // reverse gradients 
    // in - gradients with respect to the filtered field
    // out - gradient with respect to the orifinal density parametrizaiton
    virtual void RFilter(mfem::GridFunction& in, mfem::GridFunction& out) override;

private:
    
    void Initialize(mfem::Mesh* mesh_, std::map<int,double>& bndr_);
    
    std::map<int,double> bndr;

    //mfem objects
    struct{
     int order = 1;
     bool static_cond = false;
     bool pa = false;
     const char *device_config = "cpu";
     bool visualization = true;  
    } mfem_solver;

    //configuration values
    struct{
        mfem::Mesh *mesh;// mesh for the optimization
        mfem::FiniteElementCollection *fec; // H1 finite element finite element collection
        mfem::FiniteElementSpace *fespace;
        
        mfem::LinearForm *b;
        mfem::BilinearForm *a;
        
        mfem::OperatorPtr A;
        mfem::Vector B;
        mfem::Vector X;
        
        mfem::GSSmoother *prec; // GS smoother
        
        bool flag_assembl; // true of the bilinear form is assembled
                           // false otherwise
    } fltr_conf;
    
    struct{
        bool flag_own; //false if the Coefficients are owned by the class
        mfem::Coefficient *scoef_; //scalar Coefficient
        mfem::MatrixCoefficient *mcoef_; //matrix Coefficient 
    } fcoeff; 
    
    void InitCoeff();
    void FreeCoeff();
    void SetCoeff(double cf_);  // set contant coefficient
    void SetCoeff(mfem::DenseMatrix& mf_); //set constant tensor
    void SetCoeff(mfem::Coefficient& cf_);
    void SetCoeff(mfem::MatrixCoefficient& mf_);
    
};


    
    
}

#endif
