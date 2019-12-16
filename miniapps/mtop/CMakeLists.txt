# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the MFEM library. For more information and source code
# availability see http://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

if (MFEM_USE_MPI)
  add_mfem_miniapp(pmtop
    MAIN pmtop.cpp
    EXTRA_SOURCES elast_solver.cpp pde_filter.cpp
    EXTRA_HEADERS elast_solver.hpp pde_filter.hpp
    LIBRARIES mfem)

endif()


add_mfem_miniapp(smtop
	MAIN smtop.cpp
 	EXTRA_SOURCES elast_solver.cpp pde_filter.cpp
    	EXTRA_HEADERS elast_solver.hpp pde_filter.hpp
    	LIBRARIES mfem)

