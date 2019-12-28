#ifndef ELAST_SOLVER_H
#define ELAST_SOLVER_H


#include <mfem.hpp>
#include <map>

namespace mtop{
     
    
//Voigt notation for stress/strain
//[11 22 33 23 13 12] -> [0 ,1 ,2 ,3 ,4 ,5]    
class ElastStressStrain{
public:
    ElastStressStrain(){};
    virtual ~ElastStressStrain(){};
    virtual void CalcStress(double* strain_, double* stress_)=0;
    virtual void CalcStrain(double* stress_, double* strain_)=0;
    virtual void CalcComplMatrix(double* stress_, double* cmat_)=0;
    virtual void CalcStiffMatrix(double* strain_, double* smat_)=0;
    virtual double CalcEnergyL(double* stress_)=0;
    virtual double CalcEnergyU(double* strain_)=0;
    virtual double CalcEnergy(double* stress_, double* strain_)=0;
        
};

class IsotropicLinElast: public ElastStressStrain{
public:
    IsotropicLinElast();
    IsotropicLinElast(double EE_, double nu_);
    virtual ~IsotropicLinElast();
    virtual void CalcStress(double* strain, double* stress);
    virtual void CalcStrain(double* stress, double* strain);
    virtual void CalcComplMatrix(double* stress_, double* cmat_);
    virtual void CalcStiffMatrix(double* strain_, double* dmat_);
    virtual double CalcEnergyL(double* stress_);
    virtual double CalcEnergyU(double* strain_);
    virtual double CalcEnergy(double* strain_, double* stress_);
private:
    double EE;
    double nu;
};

class AnisotropicLinElast: public ElastStressStrain{
public:
    AnisotropicLinElast(double* CC_);
    AnisotropicLinElast(double* CC_, double* DD_);
    ~AnisotropicLinElast();
    
    virtual void CalcStress(double* strain_, double* stress_);
    virtual void CalcStrain(double* stress_, double* strain_);
    virtual void CalcComplMatrix(double* stress_, double* cmat);
    virtual void CalcStiffMatrix(double* strian_, double* dmat);
    virtual double CalcEnergyL(double* stress_);
    virtual double CalcEnergyU(double* strain_);
    virtual double CalcEnergy(double* strain_, double* stress_);
    
private:
    double DD[21]; //point stiffness matrix
    double CC[21]; //point compliance matrix
};

enum ForceLoadType{
  PNT=0, // point load
  VOL=1, // volumetric distributed
  SRF=2, // surface distributed
  SNO=3, // surface nornal - pressure
  ACC=4, // acceleration
  FNL=5  //last number 
};


class ForceLoad{
public:
    ForceLoad(){};
    ~ForceLoad(){};
    virtual ForceLoadType GetType()=0;
    virtual mfem::GridFunction* GetLoad()=0;
private:
    
};


class MTopElasticityIntegrator : public mfem::BilinearFormIntegrator
{
    
public:
    MTopElasticityIntegrator(std::map<int, ElastStressStrain*> *matmap_)
                                {matmap=matmap_;}
    virtual void AssembleElementMatrix(const mfem::FiniteElement &,
                                        mfem::ElementTransformation &,
                                        mfem::DenseMatrix &);
    
private:
    std::map<int, ElastStressStrain*> *matmap;
};




class ElastSolverFEMS
{
public:
    ElastSolverFEMS();
    virtual ~ElastSolverFEMS();
    
    void AddBC(int mark_, int dof_); //set BC to zero
    void AddBC(int mark_, int dof_, int val_);
    
    void AddMat(int mark_, ElastStressStrain* mat_);
    
    void AddLoad(int mark_, ForceLoad* ld_);
    
   
    
private:
    std::map<int, ElastStressStrain*> matmap;
    
    
    
};//end ElastSolverFEMS    
    
    
class  ElastSolverFEMP: public ElastSolverFEMS
{
public:
    
private:
       
};//end class ElastSolverFEMP

}

#endif
