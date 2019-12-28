#include "elast_solver.hpp"

namespace mtop{
IsotropicLinElast::IsotropicLinElast()
{
    EE=1.0;
    nu=0.20;
}
    
IsotropicLinElast::IsotropicLinElast(double EE_, double nu_)
{
    EE=EE_;
    nu=nu_;
}

IsotropicLinElast::~IsotropicLinElast()
{
    
}

void IsotropicLinElast::CalcStress(double* strain_, double* stress_)
{
    double ll=EE*nu/((1.0+nu)*(1.0-2*nu));
    double mm=EE/(2.0*(1.0+nu));
    stress_[0]=(ll+2.0*mm)*strain_[0]+ll*(strain_[1]+strain_[2]);
    stress_[1]=(ll+2.0*mm)*strain_[1]+ll*(strain_[0]+strain_[2]);
    stress_[2]=(ll+2.0*mm)*strain_[2]+ll*(strain_[0]+strain_[1]);
    stress_[3]=2.0*mm*strain_[3];
    stress_[4]=2.0*mm*strain_[4];
    stress_[5]=2.0*mm*strain_[5];
}

void IsotropicLinElast::CalcStrain(double* stress, double* strain)
{
    //GG2=2.0*GG
    double GG2=EE/((1.0+nu));
    strain[0]=stress[0]/EE-(stress[1]+stress[2])*nu/EE;
    strain[1]=stress[1]/EE-(stress[0]+stress[2])*nu/EE;
    strain[2]=stress[2]/EE-(stress[0]+stress[1])*nu/EE;
    strain[3]=stress[3]/GG2;
    strain[4]=stress[4]/GG2;
    strain[5]=stress[5]/GG2;
}

void IsotropicLinElast::CalcStiffMatrix(double* strain, double* dmat)
{
    for(int ii=0;ii<36;ii++){ dmat[ii] = 0.0;}
    
    double ll=EE*nu/((1.0+nu)*(1.0-2*nu));
    double mm2=EE/(1.0+nu); //mm2 = mu/2.0
    
    dmat[0] = (ll+mm2); dmat[1] = ll; dmat[2] = ll; 
    dmat[6] = ll; dmat[7] = (ll+mm2); dmat[8] = ll;
    dmat[12] = ll; dmat[13] = ll; dmat[14] = (ll+mm2);
    
    dmat[21] = mm2;
    dmat[28] = mm2;
    dmat[35] = mm2;
}


void IsotropicLinElast::CalcComplMatrix(double* stress, double* cmat)
{
    for(int ii=0;ii<35;ii++){ cmat[ii] = 0.0;}
    
    double iE = 1.0/EE;
    double iM = (1.0+nu)/EE;
      
    cmat[0] = iE;       cmat[1] = -nu*iE;   cmat[2] = -nu*iE;
    cmat[6] = -nu*iE;   cmat[7] = iE;       cmat[8] = -nu*iE;
    cmat[12] = -nu*iE;  cmat[13] = -nu*iE;  cmat[14] = iE;
    
    cmat[21] = iM;
    cmat[28] = iM;
    cmat[35] = iM;
}

double IsotropicLinElast::CalcEnergy(double* strain_, double* stress_)
{
    double rez = 0.0;
    for(int i=0;i<6;i++){ rez=rez+ strain_[i]*stress_[i];}
    return rez;
}

double IsotropicLinElast::CalcEnergyL(double* stress_)
{
    double rez = 0.0;
    double strain[6];
    CalcStrain(stress_,strain);
    for(int i=0;i<6;i++){ rez=rez+ strain[i]*stress_[i];}
    return rez;
}

double IsotropicLinElast::CalcEnergyU(double* strain_)
{
    double rez = 0.0;
    double stress[3];
    CalcStress(strain_,stress);
    for(int i=0;i<6;i++){ rez=rez+ strain_[i]*stress[i];}
    return rez;
}



void MTopElasticityIntegrator::AssembleElementMatrix(const mfem::FiniteElement& el,
                                                    mfem::ElementTransformation& trans,
                                                    mfem::DenseMatrix& elmat)
{
    
    int dom_num = trans.Attribute; //get the domain number
    ElastStressStrain* lmat = (*matmap)[dom_num-1]; //local material
    const mfem::IntegrationRule *ir = NULL;
    int order = 2 * trans.OrderGrad(&el) - 1; // correct order?
    ir = &mfem::IntRules.Get(el.GetGeomType(), order);
    
    int dof = el.GetDof();
    int dim = el.GetDim();
    elmat.SetSize(dof * dim);
    elmat = 0.0;

    double w;
    
    mfem::DenseMatrix dshape(dof, dim), gshape(dof, dim);
    
    int ix[3]={0,5,4};
    int iy[3]={5,1,3};
    int iz[3]={4,3,2};
    
    mfem::DenseMatrix crn(dim, dim);
    mfem::DenseMatrix DD(6, 6);
    mfem::DenseMatrix tmpm(dof, dim);
    mfem::DenseMatrix tmpk(dof, dof);
    mfem::Vector ev(6); 
    ev = 0.0;
    
    for (int i = 0; i < ir -> GetNPoints(); i++)
    {
        const mfem::IntegrationPoint &ip = ir->IntPoint(i);
        el.CalcDShape(ip, dshape);
        trans.SetIntPoint(&ip);
        w = ip.weight * trans.Weight();
        mfem::Mult(dshape, trans.InverseJacobian(), gshape);
        
        lmat->CalcStiffMatrix(ev.GetData(),DD.GetData());
        //add the ip contribution
        //xx - block
        for(int ii = 0; ii < dim; ii++){
        for(int jj = 0; jj < dim; jj++){
            crn(ii,jj) = DD(ix[ii],ix[jj]);
        }}
        
        mfem::Mult(crn,gshape,tmpm);
        mfem::MultAtB(gshape,tmpm, tmpk);
        elmat.AddMatrix(w,tmpk,0,0);
        //xy - block
        for(int ii = 0; ii < 3; ii++){
        for(int jj = 0; jj < 3; jj++){
            crn(ii,jj) = DD(ix[ii],iy[jj]);
        }}
        mfem::Mult(crn,gshape,tmpm);
        mfem::MultAtB(gshape,tmpm, tmpk);
        elmat.AddMatrix(w,tmpk,0,dof);
        tmpk.Transpose();
        elmat.AddMatrix (w,tmpk,dof,0);
        //xz - block
        for(int ii = 0; ii < 3; ii++){
        for(int jj = 0; jj < 3; jj++){
            crn(ii,jj) = DD(ix[ii],iz[jj]);
        }}
        mfem::Mult(crn,gshape,tmpm);
        mfem::MultAtB(gshape,tmpm, tmpk);
        elmat.AddMatrix(w,tmpk,0,2*dof);
        tmpk.Transpose();
        elmat.AddMatrix (w,tmpk,2*dof,0);
        //yy - block
        for(int ii = 0; ii < 3; ii++){
        for(int jj = 0; jj < 3; jj++){
            crn(ii,jj) = DD(iy[ii],iy[jj]);
        }}
        mfem::Mult(crn,gshape,tmpm);
        mfem::MultAtB(gshape,tmpm, tmpk);
        elmat.AddMatrix(w,tmpk,dof,dof);
        //yz - block
        for(int ii = 0; ii < 3; ii++){
        for(int jj = 0; jj < 3; jj++){
            crn(ii,jj) = DD(iy[ii],iz[jj]);
        }}
        mfem::Mult(crn,gshape,tmpm);
        mfem::MultAtB(gshape,tmpm, tmpk);
        elmat.AddMatrix(w,tmpk,dof,2*dof);
        tmpk.Transpose();
        elmat.AddMatrix (w,tmpk,2*dof,dof);
        //zz - block
        for(int ii = 0; ii < 3; ii++){
        for(int jj = 0; jj < 3; jj++){
            crn(ii,jj) = DD(iz[ii],iz[jj]);
        }}
        mfem::Mult(crn,gshape,tmpm);
        mfem::MultAtB(gshape,tmpm, tmpk);
        elmat.AddMatrix(w,tmpk,2*dof,2*dof);
    }
    
    
}







    
    
    
    
    
    
}
