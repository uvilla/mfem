#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <mfem.hpp>
#include <map>
#include <vector>

namespace mtop{

class HeatSolverFEMS{
public:
    HeatSolverFEMS();
    virtual ~HeatSolverFEMS();
private:    
    
};


class HeatSolverFEMP: public HeatSolverFEMS{
public:
    HeatSolverFEMP();
    virtual ~HeatSolverFEMP();
    
    
};



}
#endif
