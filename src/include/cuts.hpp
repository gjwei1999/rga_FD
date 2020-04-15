/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "branches.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "physics.hpp"

class Cuts {
 private:
  std::shared_ptr<Branches12> _data = nullptr;
  std::shared_ptr<Delta_T> _dt = nullptr;

 public:
  Cuts(const std::shared_ptr<Branches12>& data);
  Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt);
  ~Cuts();

  bool ElectronCuts();
  bool FiducialCuts();
  bool IsPip(int i);
  bool IsProton(int i);
  bool IsPim(int i);
  bool IsPositive(int i);
  bool IsFDPositive(int i);
  bool IsCDPositive(int i);
};

#endif
