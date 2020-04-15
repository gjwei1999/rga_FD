/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DT_H_GUARD
#define DT_H_GUARD
#include <cmath>
#include <iostream>
#include "branches.hpp"

class Delta_T {
 private:
  std::shared_ptr<Branches12> _data;
  float _sc_t_v = NAN;
  float _sc_r_v = NAN;
  float _vertex = NAN;

  float _ctof_t_v = NAN;
  float _ctof_r_v = NAN;
  float _ctof_vertex = NAN;

  float _sc_t = NAN;
  float _sc_r = NAN;

  float _ctof_t = NAN;
  float _ctof_r = NAN;

  float _beta = NAN;
  float _momentum = NAN;
  bool _ctof = false;

  float _vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta);
  float _deltat(int num);
  float _ctof_deltat(int pid);

 public:
  Delta_T(std::shared_ptr<Branches12> data);
  ~Delta_T();

  void dt_calc(int i);
  void dt_calc_ctof(int i);

  float dt_E(int i);
  float dt_P(int i);
  float dt_Pi(int i);
  float dt_K(int i);
  float dt_E();
  float dt_P();
  float dt_Pi();
  float dt_K();
  float dt(int pid);

  float dt_ctof_E(int i);
  float dt_ctof_P(int i);
  float dt_ctof_Pi(int i);
  float dt_ctof_K(int i);
  float dt_ctof_E();
  float dt_ctof_P();
  float dt_ctof_Pi();
  float dt_ctof_K();
  float dt_ctof(int pid);

  float momentum();
  bool ctof();
};

#endif
