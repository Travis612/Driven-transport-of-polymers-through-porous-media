/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(test/bc2,FixTestBc2)

#else

#ifndef LMP_FIX_TEST_BC2_H
#define LMP_FIX_TEST_BC2_H


#include "fix.h"

namespace LAMMPS_NS {

class FixTestBc2 : public Fix {
 public:
  FixTestBc2(class LAMMPS *, int, char **);
  virtual ~FixTestBc2() ;
  int setmask();
  virtual void init();
  void init_list(int, class NeighList *);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  virtual void reset_dt();
  double memory_usage();

 protected:
  double dtv,dtf,rhow;
  int mass_require;
  int maxatom;
  double **nw;
  class NeighList *list;

  


  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
