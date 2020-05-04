/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = std::make_shared<TCanvas>("def");
  if (getenv("BEAM_E") != NULL) {
    if (atof(getenv("BEAM_E")) < 3) {
      q2_max = 1.0;
      w_max = 3.5;
      p_max = 3.0;
    } else if (atof(getenv("BEAM_E")) < 8) {
      q2_max = 4.0;
      w_max = 4.0;
      p_max = 4.0;
    } else if (atof(getenv("BEAM_E")) < 9) {
      q2_max = 7.0;
      w_max = 7.0;
      p_max = 7.0;
    }
  }
  const Int_t ndims = 7;

  Int_t bins_[ndims] = {25, bins, bins, bins, bins, bins, bins};
  Double_t xmin[ndims] = {p_min, p_min, 1.1, 0.2, 0., 0, 0};
  Double_t xmax[ndims] = {w_max, q2_max, 3.0, 2.2, 180, 360, 360};

  // Kinematic
  sevenDHist = new THnSparseF("7D_hist", "7D_hist", ndims, bins_, xmin, xmax);

  momentum = std::make_shared<TH1D>("mom", "mom", bins, p_min, p_max);
  W_hist = std::make_shared<TH1D>("W", "W", bins, zero, w_max);

  Theta_vs_mom_x_mu = std::make_shared<TH2D>("theta_vs_mom_x_mu_all_sec", "theta_vs_mom_x_mu_all_sec", bins, zero,
                                             p_max, bins, zero, 120);

  W_P2pi_hist = std::make_shared<TH1D>("W_P2pi", "W_P2pi", bins, zero, w_max);

  Q2_hist = std::make_shared<TH1D>("Q2", "Q2", bins, zero, q2_max);
  W_vs_q2 = std::make_shared<TH2D>("W_vs_q2", "W_vs_q2", bins, zero, w_max, bins, zero, q2_max);

  W_thrown = std::make_shared<TH1D>("W_thrown", "W_thrown", bins, zero, w_max);
  W_vs_Q2_thrown = std::make_shared<TH2D>("W_vs_q2_thrown", "W_vs_q2_thrown", bins, zero, w_max, bins, zero, q2_max);

  theta_prot = std::make_shared<TH1D>("theta_prot", "theta_prot", bins, 0.5, 180);
  theta_pip = std::make_shared<TH1D>("theta_pip", "theta_pip", bins, 0.5, 180);
  theta_pim = std::make_shared<TH1D>("theta_pim", "theta_pim", bins, 0.5, 180);

  Phi_gamma = std::make_shared<TH1D>("Phi_gamma", "Phi_gamma", bins, 0, 360);
  Phi_prot = std::make_shared<TH1D>("Phi_prot", "Phi_prot", bins, 0, 360);
  Phi_pip = std::make_shared<TH1D>("Phi_pip", "Phi_pip", bins, 0, 360);
  Phi_pim = std::make_shared<TH1D>("Phi_pim", "Phi_pim", bins, 0, 360);

  alpha_prot = std::make_shared<TH1D>("alpha_prot", "alpha_prot", bins, 0.5, 360);
  alpha_pip = std::make_shared<TH1D>("alpha_pip", "alpha_pip", bins, 0.5, 360);
  alpha_pim = std::make_shared<TH1D>("alpha_pim", "alpha_pim", bins, 0.5, 360);

  theta_prot_mc = std::make_shared<TH1D>("theta_prot_mc", "theta_prot_mc", bins, 0, 180);
  theta_pip_mc = std::make_shared<TH1D>("theta_pip_mc", "theta_pip_mc", bins, 0, 180);
  theta_pim_mc = std::make_shared<TH1D>("theta_pim_mc", "theta_pim_mc", bins, 0, 180);

  Phi_gamma_mc = std::make_shared<TH1D>("Phi_gamma_mc", "Phi_gamma_mc", bins, 0, 360);
  Phi_prot_mc = std::make_shared<TH1D>("Phi_prot_mc", "Phi_prot_mc", bins, 0, 360);
  Phi_pip_mc = std::make_shared<TH1D>("Phi_pip_mc", "Phi_pip_mc", bins, 0, 360);
  Phi_pim_mc = std::make_shared<TH1D>("Phi_pim_mc", "Phi_pim_mc", bins, 0, 360);

  alpha_prot_mc = std::make_shared<TH1D>("alpha_prot_mc", "alpha_prot_mc", bins, 0, 360);
  alpha_pip_mc = std::make_shared<TH1D>("alpha_pip_mc", "alpha_pip_mc", bins, 0, 360);
  alpha_pim_mc = std::make_shared<TH1D>("alpha_pim_mc", "alpha_pim_mc", bins, 0, 360);

  theta_prot_thrown = std::make_shared<TH1D>("theta_prot_thrown", "theta_prot_thrown", bins, 0, 180);
  theta_pip_thrown = std::make_shared<TH1D>("theta_pip_thrown", "theta_pip_thrown", bins, 0, 180);
  theta_pim_thrown = std::make_shared<TH1D>("theta_pim_thrown", "theta_pim_thrown", bins, 0, 180);

  Phi_gamma_thrown = std::make_shared<TH1D>("Phi_gamma_thrown", "Phi_gamma_thrown", bins, 0, 360);
  Phi_prot_thrown = std::make_shared<TH1D>("Phi_prot_thrown", "Phi_prot_thrown", bins, 0, 360);
  Phi_pip_thrown = std::make_shared<TH1D>("Phi_pip_thrown", "Phi_pip_thrown", bins, 0, 360);
  Phi_pim_thrown = std::make_shared<TH1D>("Phi_pim_thrown", "Phi_pim_thrown", bins, 0, 360);

  alpha_prot_thrown = std::make_shared<TH1D>("alpha_prot_thrown", "alpha_prot_thrown", bins, 0, 360);
  alpha_pip_thrown = std::make_shared<TH1D>("alpha_pip_thrown", "alpha_pip_thrown", bins, 0, 360);
  alpha_pim_thrown = std::make_shared<TH1D>("alpha_pim_thrown", "alpha_pim_thrown", bins, 0, 360);

  MM_neutron = std::make_shared<TH1D>("missMass", "missMass", bins, zero, 4.0);
  MM_twoPi = std::make_shared<TH1D>("missMass_twoPi", "missMass_twoPi", bins, -0.4, 0.4);
  MM2_twoPi = std::make_shared<TH1D>("missMass_SQ_twoPi", "missMass_SQ_twoPi", bins, -0.2, 0.2);
  W_hist_twoPi = std::make_shared<TH1D>("W_twoPi", "W_twoPi", bins, zero, w_max);
  Q2_hist_twoPi = std::make_shared<TH1D>("Q2_twoPi", "Q2_twoPi", bins, zero, q2_max);
  W_vs_q2_twoPi = std::make_shared<TH2D>("W_vs_q2_twoPi", "W_vs_q2_twoPi", bins, zero, w_max, bins, zero, q2_max);
  
  W_q2_15_twoPi = std::make_shared<TH1D>("W_q2_15_twoPi", "W_q2_15_twoPi", bins, zero, w_max);
  W_q2_2_twoPi = std::make_shared<TH1D>("W_q2_2_twoPi", "W_q2_2_twoPi", bins, zero, w_max);
  W_q2_23_twoPi = std::make_shared<TH1D>("W_q2_23_twoPi", "W_q2_23_twoPi", bins, zero, w_max);
  W_q2_34_twoPi = std::make_shared<TH1D>("W_q2_34_twoPi", "W_q2_34_twoPi", bins, zero, w_max);
  
  W_q2_15_photon_twoPi = std::make_shared<TH1D>("W_q2_15_photon_twoPi", "W_q2_15_photon_twoPi", bins, zero, w_max);
  W_q2_2_photon_twoPi = std::make_shared<TH1D>("W_q2_2_photon_twoPi", "W_q2_2_photon_twoPi", bins, zero, w_max);
  W_q2_23_photon_twoPi = std::make_shared<TH1D>("W_q2_23_photon_twoPi", "W_q2_23_photon_twoPi", bins, zero, w_max);
  W_q2_34_photon_twoPi = std::make_shared<TH1D>("W_q2_34_photon_twoPi", "W_q2_34_photon_twoPi", bins, zero, w_max);
  /////
  /////
  /////modified by Jiawei
  invmass_Ppip_twoPi = std::make_shared<TH1D>("invmass_Ppip_twoPi", "invmass_Ppip_twoPi", bins, 1.0, 3.5);
  invmass_Ppim_twoPi = std::make_shared<TH1D>("invmass_Ppim_twoPi", "invmass_Ppim_twoPi", bins, 1.0, 3.5);
  invmass_pippim_twoPi = std::make_shared<TH1D>("invmass_pippim_twoPi", "invmass_pippim_twoPi", bins, 0.0, 2.5);
  
  theta_Ppip_vs_invmass_twoPi = std::make_shared<TH2D>("theta_Ppip_vs_invmass_twoPi", "theta_Ppip_vs_invmass_twoPi", bins, 1.0, 3.5, bins, 0, 180);
  theta_Ppim_vs_invmass_twoPi = std::make_shared<TH2D>("theta_Ppim_vs_invmass_twoPi", "theta_Ppim_vs_invmass_twoPi", bins, 1.0, 3.5, bins, 0, 180);
  theta_pippim_vs_invmass_twoPi = std::make_shared<TH2D>("theta_pippim_vs_invmass_twoPi", "theta_pippim_vs_invmass_twoPi", bins, 0.0, 2.5, bins, 0, 180);
  
  p_vs_phi_Prot_twoPi = std::make_shared<TH2D>("p_vs_phi_Prot_twoPi", "p_vs_phi_Prot_twoPi", bins, -180, 180, bins, 0, 8.0);
  p_vs_phi_pip_twoPi = std::make_shared<TH2D>("p_vs_phi_pip_twoPi", "p_vs_phi_pip_twoPi", bins, -180, 180, bins, 0, 8.0);
  p_vs_phi_pim_twoPi = std::make_shared<TH2D>("p_vs_phi_pim_twoPi", "p_vs_phi_pim_twoPi", bins, -180, 180, bins, 0, 8.0);

  p_vs_theta_Prot_twoPi = std::make_shared<TH2D>("p_vs_theta_Prot_twoPi", "p_vs_theta_Prot_twoPi", bins, 0, 180, bins, 0, 8.0);
  p_vs_theta_pip_twoPi = std::make_shared<TH2D>("p_vs_theta_pip_twoPi", "p_vs_theta_pip_twoPi", bins, 0, 180, bins, 0, 8.0);
  p_vs_theta_pim_twoPi = std::make_shared<TH2D>("p_vs_theta_pim_twoPi", "p_vs_theta_pim_twoPi", bins, 0, 180, bins, 0, 8.0);

  theta_vs_phi_Prot_twoPi = std::make_shared<TH2D>("theta_vs_phi_Prot_twoPi", "theta_vs_phi_Prot_twoPi", bins, -180, 180, bins, 0, 180);
  theta_vs_phi_pip_twoPi = std::make_shared<TH2D>("theta_vs_phi_pip_twoPi", "theta_vs_phi_pip_twoPi", bins, -180, 180, bins, 0, 180);
  theta_vs_phi_pim_twoPi = std::make_shared<TH2D>("theta_vs_phi_pim_twoPi", "theta_vs_phi_pim_twoPi", bins, -180, 180, bins, 0, 180);
/////
  /////
  /////

  W_hist_singlePi = std::make_shared<TH1D>("W_singlePi", "W_singlePi", bins, zero, w_max);
  Q2_hist_singlePi = std::make_shared<TH1D>("Q2_singlePi", "Q2_singlePi", bins, zero, q2_max);
  W_vs_q2_singlePi =
      std::make_shared<TH2D>("W_vs_q2_singlePi", "W_vs_q2_singlePi", bins, zero, w_max, bins, zero, q2_max);

  makeHists_sector();
  makeHists_deltat();
  makeHists_MomVsBeta();
  makeHists_electron_cuts();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;
  
 //  Write_EC();
  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  //TDirectory* WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  //WvsQ2_folder->cd();
  Write_WvsQ2();


  // THnSparse 7D
  std::cerr << BOLDBLUE << " Hists_7D()" << DEF << std::endl;
  TDirectory* THnSparse_7D_folder = RootOutputFile->mkdir("THnSparse_7D");
  THnSparse_7D_folder->cd();
  writeHists7D();

//  //  Write_EC();
//  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
//  TDirectory* WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
//  WvsQ2_folder->cd();
//  Write_WvsQ2();

  std::cerr << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
  TDirectory* Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
  Write_MomVsBeta_folder->cd();
  Write_MomVsBeta();

  std::cerr << BOLDBLUE << "Write_Electron_cuts()" << DEF << std::endl;
  TDirectory* Electron_Cuts = RootOutputFile->mkdir("Electron_Cuts");
  Electron_Cuts->cd();
  Write_Electron_cuts();

  std::cerr << BOLDBLUE << "Write_deltat()" << DEF << std::endl;
  // TDirectory* Write_deltat_folder = RootOutputFile->mkdir("Delta_t");
  // Write_deltat_folder->cd();
  Write_deltat();

  std::cerr << BOLDBLUE << "SF_CUTS()" << DEF << std::endl;
  TDirectory* SF_CUTS = RootOutputFile->mkdir("SF_CUTS");
  SF_CUTS->cd();
  Write_SF();

  std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::Fill_histSevenD(const std::shared_ptr<Reaction>& _e) {
  // fill it
  const Int_t ndims = 7;
  Double_t x[ndims];
  x[0] = _e->W();
  x[1] = _e->Q2();
  x[2] = _e->inv_Ppip();
  x[3] = _e->inv_pip_pim();
  x[4] = _e->pim_theta();
  x[5] = _e->pim_Phi();
  x[6] = _e->alpha_ppip_pipim();

  TThread::Lock();
  sevenDHist->Fill(x);
  TThread::UnLock();
  sevenDHist->GetNbins();
}
void Histogram::writeHists7D() { sevenDHist->Write(); }

void Histogram::Fill_x_mu(const std::shared_ptr<Reaction>& _e) {
  E_x_mu->Fill(_e->E_x_mu());
  P_x_mu->Fill(_e->P_x_mu());
  diff_E_P_x_mu->Fill(-_e->E_x_mu() + _e->P_x_mu());
  mom_vs_E_x_mu->Fill(_e->E_x_mu(), _e->P_x_mu());
  diff_theta_in_x_mu->Fill(-_e->theta_beam() + _e->theta_x_mu());
  Theta_vs_mom_x_mu->Fill(_e->P_x_mu(), _e->theta_x_mu());
}
void Histogram::Fill_SF(const std::shared_ptr<Branches12>& _d) {
  sf_hist->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}
void Histogram::Write_SF() {
  // set stat

  sf_hist->SetOption("COLZ");
  gStyle->SetOptFit(1111);
  sf_hist->Write();

  E_x_mu->SetXTitle("Energy comp (GeV)");
  E_x_mu->Write();

  P_x_mu->SetXTitle("3 momentum comp (GeV)");
  P_x_mu->Write();

  diff_E_P_x_mu->SetXTitle("diif (E-P_mag) (GeV)");
  diff_E_P_x_mu->Write();

  diff_theta_in_x_mu->SetXTitle("theta_x_mu - theta_beam");
  diff_theta_in_x_mu->Write();

  mom_vs_E_x_mu->SetXTitle("E comp (GeV)");
  mom_vs_E_x_mu->SetYTitle("P comp (GeV)");
  mom_vs_E_x_mu->SetOption("COLZ");
  mom_vs_E_x_mu->Write();

  Dthtea_vs_Dphi->SetXTitle("#Delata#Theta (deg)");
  Dthtea_vs_Dphi->SetYTitle("#Delata#Phi (deg)");
  Dthtea_vs_Dphi->SetOption("COLZ");
  Dthtea_vs_Dphi->Write();
}

void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction>& _e) {
  W_P2pi_hist->Fill(_e->w_P2pi_rec(), _e->weight());
  Q2_hist->Fill(_e->Q2(), _e->weight());

  // short sec = _e->sec();
  // if (sec > 0 && sec <= 6) {
  //         W_vs_q2_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
  //         W_sec[sec - 1]->Fill(_e->W()), _e->weight();
  // }

  short det = _e->det();
  if (det == 1 && _e->W() <= 3.5) {
    W_det[0]->Fill(_e->W(), _e->weight());
    WQ2_det[0]->Fill(_e->W(), _e->Q2(), _e->weight());
  } else if (det == 2) {
    W_det[1]->Fill(_e->W(), _e->weight());
    WQ2_det[1]->Fill(_e->W(), _e->Q2(), _e->weight());
  } else {
    W_det[2]->Fill(_e->W(), _e->weight());
    WQ2_det[2]->Fill(_e->W(), _e->Q2(), _e->weight());
  }
}
void Histogram::Fill_WvsQ2(const std::shared_ptr<MCReaction>& _e) {
  W_vs_q2->Fill(_e->W(), _e->Q2(), _e->weight());
  W_hist->Fill(_e->W(), _e->weight());
  Q2_hist->Fill(_e->Q2(), _e->weight());

  short sec = _e->sec();
  if (sec > 0 && sec <= 6) {
    W_vs_q2_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
    W_sec[sec - 1]->Fill(_e->W(), _e->weight());
  }

  short det = _e->det();
  if (det == 1 && _e->W() <= 3.5) {
    W_det[0]->Fill(_e->W(), _e->weight());
    WQ2_det[0]->Fill(_e->W(), _e->Q2(), _e->weight());
  } else if (det == 2) {
    W_det[1]->Fill(_e->W(), _e->weight());
    WQ2_det[1]->Fill(_e->W(), _e->Q2(), _e->weight());
  } else {
    W_det[2]->Fill(_e->W(), _e->weight());
    WQ2_det[2]->Fill(_e->W(), _e->Q2(), _e->weight());
  }
}

void Histogram::Fill_WvsQ2_singlePi(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();

  W_vs_q2_singlePi->Fill(_e->W(), _e->Q2(), _e->weight());
  W_hist_singlePi->Fill(_e->W(), _e->weight());
  Q2_hist_singlePi->Fill(_e->Q2(), _e->weight());
  MM_neutron->Fill(_e->MM(), _e->weight());
  if (sec > 0 && sec <= 6) {
    W_vs_MM_singlePi[sec - 1]->Fill(_e->W(), _e->MM(), _e->weight());
    W_vs_q2_singlePi_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
    W_singlePi_sec[sec - 1]->Fill(_e->W(), _e->weight());
    MM_neutron_sec[sec - 1]->Fill(_e->MM(), _e->weight());
  }
}

void Histogram::Fill_WvsQ2_twoPi(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();

  theta_prot_mc->Fill(_e->prot_theta(), _e->weight());
  theta_pip_mc->Fill(_e->pip_theta(), _e->weight());
  theta_pim_mc->Fill(_e->pim_theta(), _e->weight());

  Phi_gamma_mc->Fill(_e->gamma_Phi(), _e->weight());
  Phi_prot_mc->Fill(_e->prot_Phi(), _e->weight());
  Phi_pip_mc->Fill(_e->pip_Phi(), _e->weight());
  Phi_pim_mc->Fill(_e->pim_Phi(), _e->weight());

  alpha_prot_mc->Fill(_e->alpha_pippim_pipf(), _e->weight());
  alpha_pip_mc->Fill(_e->alpha_ppim_pipip(), _e->weight());
  alpha_pim_mc->Fill(_e->alpha_ppip_pipim(), _e->weight());

  W_vs_q2_twoPi->Fill(_e->W(), _e->Q2(), _e->weight());
  W_hist_twoPi->Fill(_e->W(), _e->weight());
  Q2_hist_twoPi->Fill(_e->Q2(), _e->weight());
  MM_twoPi->Fill(_e->MM(), _e->weight());
  MM2_twoPi->Fill(_e->MM2(), _e->weight());
  
  if(_e->Q2()<1.5){
  W_q2_15_twoPi->Fill(_e->W(), _e->weight());
  W_q2_15_photon_twoPi->Fill(_e->W(), _e->weight()/_e->flux());
  }
  if(_e->Q2()<2.0){
  W_q2_2_twoPi->Fill(_e->W(), _e->weight());
  W_q2_2_photon_twoPi->Fill(_e->W(), _e->weight()/_e->flux());
  }
  else if(_e->Q2()>2.0 && _e->Q2()<3.0){
  W_q2_23_twoPi->Fill(_e->W(), _e->weight());
  W_q2_23_photon_twoPi->Fill(_e->W(), _e->weight()/_e->flux());
  }
  else if(_e->Q2()>3.0 && _e->Q2()<4.0){
  W_q2_34_twoPi->Fill(_e->W(), _e->weight());
  W_q2_34_photon_twoPi->Fill(_e->W(), _e->weight()/_e->flux());
  }
  //
  //
  //modified by Jiawei
  invmass_Ppip_twoPi->Fill(_e->inv_Ppip(), _e->weight());
  invmass_Ppim_twoPi->Fill( _e->inv_Ppim(), _e->weight());
  invmass_pippim_twoPi->Fill(_e->inv_pip_pim(), _e->weight());
  
  theta_Ppip_vs_invmass_twoPi->Fill(_e->inv_Ppip(), _e->theta_Ppip_lab(), _e->weight());
  theta_Ppim_vs_invmass_twoPi->Fill( _e->inv_Ppim(), _e->theta_Ppim_lab(), _e->weight());
  theta_pippim_vs_invmass_twoPi->Fill( _e->inv_pip_pim(), _e->theta_pippim_lab(), _e->weight());
  
  p_vs_phi_Prot_twoPi->Fill(_e->prot_phi_lab(), _e->prot_p(), _e->weight());
  p_vs_phi_pip_twoPi->Fill( _e->pip_phi_lab(), _e->pip_p(), _e->weight());
  p_vs_phi_pim_twoPi->Fill( _e->pim_phi_lab(), _e->pim_p(), _e->weight());

  p_vs_theta_Prot_twoPi->Fill(_e->prot_theta_lab(), _e->prot_p(), _e->weight());
  p_vs_theta_pip_twoPi->Fill( _e->pip_theta_lab(), _e->pip_p(), _e->weight());
  p_vs_theta_pim_twoPi->Fill( _e->pim_theta_lab(), _e->pim_p(), _e->weight());
  
  theta_vs_phi_Prot_twoPi->Fill(_e->prot_phi_lab(), _e->prot_theta_lab(), _e->weight());
  theta_vs_phi_pip_twoPi->Fill( _e->pip_phi_lab(), _e->pip_theta_lab(), _e->weight());
  theta_vs_phi_pim_twoPi->Fill( _e->pim_phi_lab(), _e->pim_theta_lab(), _e->weight());
  //
  //
  //

  
  //theta_prot->Fill(_e->prot_theta(), _e->weight());
  //theta_pip->Fill(_e->pip_theta(), _e->weight());
  //theta_pim->Fill(_e->pim_theta(), _e->weight());
  
  //
  //
  //modified by Jiawei
  theta_prot->Fill(_e->prot_theta_lab(), _e->weight());
  theta_pip->Fill(_e->pip_theta_lab(), _e->weight());
  theta_pim->Fill(_e->pim_theta_lab(), _e->weight());
  //
  //
  //


  if (_e->gamma_Phi() != 0) Phi_gamma->Fill(_e->gamma_Phi(), _e->weight());
  if (_e->prot_Phi() != 0) Phi_prot->Fill(_e->prot_Phi(), _e->weight());
  if (_e->pip_Phi() != 0) Phi_pip->Fill(_e->pip_Phi(), _e->weight());
  if (_e->pim_Phi() != 0) Phi_pim->Fill(_e->pim_Phi(), _e->weight());

  alpha_prot->Fill(_e->alpha_pippim_pipf(), _e->weight());
  alpha_pip->Fill(_e->alpha_ppim_pipip(), _e->weight());
  alpha_pim->Fill(_e->alpha_ppip_pipim(), _e->weight());

  if (sec > 0 && sec <= 6) {
    W_vs_MM_twoPi[sec - 1]->Fill(_e->W(), _e->MM(), _e->weight());
    W_vs_q2_twoPi_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
    W_twoPi_sec[sec - 1]->Fill(_e->W(), _e->weight());
    MM_twoPi_sec[sec - 1]->Fill(_e->MM(), _e->weight());
    MM2_twoPi_sec[sec - 1]->Fill(_e->MM2(), _e->weight());
  }
  if (_e->W() < 2.5) {
    if (_e->Q2() < 4.5) {
      inv_mass_P_pip[0][0]->Fill(_e->inv_Ppip(), _e->weight());
      inv_mass_P_pim[0][0]->Fill(_e->inv_Ppim(), _e->weight());
      inv_mass_pip_pim[0][0]->Fill(_e->inv_pip_pim(), _e->weight());
      // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI <<  '\n';
      theta_P_vs_mass_pip_pim[0][0]->Fill(_e->inv_pip_pim(), _e->prot_theta(), _e->weight());
      theta_pip_vs_mass_Ppim[0][0]->Fill(_e->inv_Ppim(), _e->pip_theta(), _e->weight());
      theta_pim_vs_mass_Ppip[0][0]->Fill(_e->inv_Ppip(), _e->pim_theta(), _e->weight());

      theta_P_lab_vs_mass_pip_pim[0][0]->Fill(_e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim[0][0]->Fill(_e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip[0][0]->Fill(_e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
    } else if (_e->Q2() > 4.5) {
      inv_mass_P_pip[0][1]->Fill(_e->inv_Ppip(), _e->weight());
      inv_mass_P_pim[0][1]->Fill(_e->inv_Ppim(), _e->weight());
      inv_mass_pip_pim[0][1]->Fill(_e->inv_pip_pim(), _e->weight());
      // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI <<  '\n';
      theta_P_vs_mass_pip_pim[0][1]->Fill(_e->inv_pip_pim(), _e->prot_theta(), _e->weight());
      theta_pip_vs_mass_Ppim[0][1]->Fill(_e->inv_Ppim(), _e->pip_theta(), _e->weight());
      theta_pim_vs_mass_Ppip[0][1]->Fill(_e->inv_Ppip(), _e->pim_theta(), _e->weight());

      theta_P_lab_vs_mass_pip_pim[0][1]->Fill(_e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim[0][1]->Fill(_e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip[0][1]->Fill(_e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
    }
  } else if (_e->W() > 2.5) {
    if (_e->Q2() < 4.5) {
      inv_mass_P_pip[1][0]->Fill(_e->inv_Ppip(), _e->weight());
      inv_mass_P_pim[1][0]->Fill(_e->inv_Ppim(), _e->weight());
      inv_mass_pip_pim[1][0]->Fill(_e->inv_pip_pim(), _e->weight());
      theta_P_vs_mass_pip_pim[1][0]->Fill(_e->inv_pip_pim(), _e->prot_theta(), _e->weight());
      theta_pip_vs_mass_Ppim[1][0]->Fill(_e->inv_Ppim(), _e->pip_theta(), _e->weight());
      theta_pim_vs_mass_Ppip[1][0]->Fill(_e->inv_Ppip(), _e->pim_theta(), _e->weight());

      theta_P_lab_vs_mass_pip_pim[1][0]->Fill(_e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim[1][0]->Fill(_e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip[1][0]->Fill(_e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
    } else if (_e->Q2() > 4.5) {
      inv_mass_P_pip[1][1]->Fill(_e->inv_Ppip(), _e->weight());
      inv_mass_P_pim[1][1]->Fill(_e->inv_Ppim(), _e->weight());
      inv_mass_pip_pim[1][1]->Fill(_e->inv_pip_pim(), _e->weight());
      // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI <<  '\n';
      theta_P_vs_mass_pip_pim[1][1]->Fill(_e->inv_pip_pim(), _e->prot_theta(), _e->weight());
      theta_pip_vs_mass_Ppim[1][1]->Fill(_e->inv_Ppim(), _e->pip_theta(), _e->weight());
      theta_pim_vs_mass_Ppip[1][1]->Fill(_e->inv_Ppip(), _e->pim_theta(), _e->weight());

      theta_P_lab_vs_mass_pip_pim[1][1]->Fill(_e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim[1][1]->Fill(_e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip[1][1]->Fill(_e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
    }
  }
}

void Histogram::Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<MCReaction>& _e) {
  short sec = _e->sec();

  // W_vs_q2_twoPi_thrown->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
  // W_hist_twoPi_thrown->Fill(_e->W_mc(), _e->weight());
  // Q2_hist_twoPi_thrown->Fill(_e->Q2_mc(), _e->weight());
  // MM_twoPi_thrown->Fill(_e->MM_mc(), _e->weight());
  // MM2_twoPi_thrown->Fill(_e->MM2_mc(), _e->weight());
  W_vs_Q2_thrown->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
  W_thrown->Fill(_e->W_mc(), _e->weight());

  theta_prot_thrown->Fill(_e->MCprot_theta_thrown(), _e->weight());
  theta_pip_thrown->Fill(_e->MCpip_theta_thrown(), _e->weight());
  theta_pim_thrown->Fill(_e->MCpim_theta_thrown(), _e->weight());

  Phi_gamma_thrown->Fill(_e->MCgamma_Phi_thrown(), _e->weight());
  Phi_prot_thrown->Fill(_e->MCprot_Phi_thrown(), _e->weight());
  Phi_pip_thrown->Fill(_e->MCpip_Phi_thrown(), _e->weight());
  Phi_pim_thrown->Fill(_e->MCpim_Phi_thrown(), _e->weight());

  alpha_prot_thrown->Fill(_e->MCalpha_pippim_pipf_thrown(), _e->weight());
  alpha_pip_thrown->Fill(_e->MCalpha_ppim_pipip_thrown(), _e->weight());
  alpha_pim_thrown->Fill(_e->MCalpha_ppip_pipim_thrown(), _e->weight());

  if (sec > 0 && sec <= 6) {
    W_vs_MM_twoPi_thrown[sec - 1]->Fill(_e->W_mc(), _e->MM_mc(), _e->weight());
    W_vs_q2_twoPi_sec_thrown[sec - 1]->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
    W_twoPi_sec_thrown[sec - 1]->Fill(_e->W_mc(), _e->weight());
    MM_twoPi_sec_thrown[sec - 1]->Fill(_e->MM_mc(), _e->weight());
    MM2_twoPi_sec_thrown[sec - 1]->Fill(_e->MM2_mc(), _e->weight());
  }
  if (_e->W_mc() < 2.5) {
    if (_e->Q2_mc() < 4.5) {
      inv_mass_P_pip_thrown[0][0]->Fill(_e->MCinv_Ppip(), _e->weight());
      inv_mass_P_pim_thrown[0][0]->Fill(_e->MCinv_Ppim(), _e->weight());
      inv_mass_pip_pim_thrown[0][0]->Fill(_e->MCinv_pip_pim(), _e->weight());
      theta_P_vs_mass_pip_pim_thrown[0][0]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta_thrown(), _e->weight());
      theta_pip_vs_mass_Ppim_thrown[0][0]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta_thrown(), _e->weight());
      theta_pim_vs_mass_Ppip_thrown[0][0]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta_thrown(), _e->weight());
      theta_P_lab_vs_mass_pip_pim_thrown[0][0]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim_thrown[0][0]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip_thrown[0][0]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta_lab(), _e->weight());
    } else if (_e->Q2_mc() > 4.5) {
      inv_mass_P_pip_thrown[0][1]->Fill(_e->MCinv_Ppip(), _e->weight());
      inv_mass_P_pim_thrown[0][1]->Fill(_e->MCinv_Ppim(), _e->weight());
      inv_mass_pip_pim_thrown[0][1]->Fill(_e->MCinv_pip_pim(), _e->weight());
      theta_P_vs_mass_pip_pim_thrown[0][1]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta(), _e->weight());
      theta_pip_vs_mass_Ppim_thrown[0][1]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta(), _e->weight());
      theta_pim_vs_mass_Ppip_thrown[0][1]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta(), _e->weight());
      theta_P_lab_vs_mass_pip_pim_thrown[0][1]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim_thrown[0][1]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip_thrown[0][1]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta_lab(), _e->weight());
    }

  } else if (_e->W_mc() > 2.5) {
    if (_e->Q2_mc() < 4.5) {
      inv_mass_P_pip_thrown[1][0]->Fill(_e->MCinv_Ppip(), _e->weight());
      inv_mass_P_pim_thrown[1][0]->Fill(_e->MCinv_Ppim(), _e->weight());
      inv_mass_pip_pim_thrown[1][0]->Fill(_e->MCinv_pip_pim(), _e->weight());
      // std::cout << "inv_macc_pip_pim_mc " << _e->MCinv_pip_pim() << '\n';
      theta_P_vs_mass_pip_pim_thrown[1][0]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta(), _e->weight());
      theta_pip_vs_mass_Ppim_thrown[1][0]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta(), _e->weight());
      theta_pim_vs_mass_Ppip_thrown[1][0]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta(), _e->weight());
      theta_P_lab_vs_mass_pip_pim_thrown[1][0]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim_thrown[1][0]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip_thrown[1][0]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta_lab(), _e->weight());
    } else if (_e->Q2_mc() > 4.5) {
      inv_mass_P_pip_thrown[1][1]->Fill(_e->MCinv_Ppip(), _e->weight());
      inv_mass_P_pim_thrown[1][1]->Fill(_e->MCinv_Ppim(), _e->weight());
      inv_mass_pip_pim_thrown[1][1]->Fill(_e->MCinv_pip_pim(), _e->weight());
      theta_P_vs_mass_pip_pim_thrown[1][1]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta(), _e->weight());
      theta_pip_vs_mass_Ppim_thrown[1][1]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta(), _e->weight());
      theta_pim_vs_mass_Ppip_thrown[1][1]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta(), _e->weight());
      theta_P_lab_vs_mass_pip_pim_thrown[1][1]->Fill(_e->MCinv_pip_pim(), _e->MCprot_theta_lab(), _e->weight());
      theta_pip_lab_vs_mass_Ppim_thrown[1][1]->Fill(_e->MCinv_Ppim(), _e->MCpip_theta_lab(), _e->weight());
      theta_pim_lab_vs_mass_Ppip_thrown[1][1]->Fill(_e->MCinv_Ppip(), _e->MCpim_theta_lab(), _e->weight());
    }
  }
}

// W and Q^2
void Histogram::Fill_WvsQ2_Npip(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();
  if (sec > 0 && sec <= 6) {
    W_vs_q2_Npip_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
    W_Npip_sec[sec - 1]->Fill(_e->W(), _e->weight());
    MM_Npip_sec[sec - 1]->Fill(_e->MM(), _e->weight());
  }
}

void Histogram::Write_WvsQ2() {
  for (short i = 0; i < 3; i++) {
    WQ2_det[i]->SetXTitle("W (GeV)");
    WQ2_det[i]->SetYTitle("Q^{2} (GeV^2)");
    WQ2_det[i]->SetOption("COLZ1");
    if (WQ2_det[i]->GetEntries()) WQ2_det[i]->Write();
    W_det[i]->SetXTitle("W (GeV)");
    if (W_det[i]->GetEntries()) W_det[i]->Write();
  }
  auto WvsQ2_can = std::make_unique<TCanvas>("WvsQ2_can", "W vs Q2 sectors", 1920, 1080);
  WvsQ2_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_sec[i]->SetOption("COLZ1");
    WvsQ2_can->cd(i + 1);
    W_vs_q2_sec[i]->Draw("same");
  }
  WvsQ2_can->Write();

  auto W_can = std::make_unique<TCanvas>("W_can", "W sectors", 1920, 1080);
  W_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
    W_sec[i]->SetXTitle("W (GeV)");
    W_can->cd(i + 1);

    W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.85, 1.05);
    gStyle->SetOptFit(01);
    W_sec[i]->Draw("same");
  }
  W_can->Write();

  W_vs_q2->SetXTitle("W (GeV)");
  W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2->SetOption("COLZ1");
  if (W_vs_q2->GetEntries()) W_vs_q2->Write();

  W_hist->SetXTitle("W (GeV)");
  if (W_hist->GetEntries()) W_hist->Write();

  W_P2pi_hist->SetXTitle("W_P2pi (GeV)");
  if (W_P2pi_hist->GetEntries()) W_P2pi_hist->Write();

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  if (Q2_hist->GetEntries()) Q2_hist->Write();

  W_vs_Q2_thrown->SetXTitle("W (GeV)");
  W_vs_Q2_thrown->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_Q2_thrown->SetOption("COLZ1");
  if (W_vs_Q2_thrown->GetEntries()) W_vs_Q2_thrown->Write();

  W_thrown->SetXTitle("W_thrown (GeV)");
  if (W_thrown->GetEntries()) W_thrown->Write();

  theta_prot->SetXTitle("theta_prot (Deg)");
  if (theta_prot->GetEntries()) theta_prot->Write();
  theta_pip->SetXTitle("theta_pip (Deg)");
  if (theta_pip->GetEntries()) theta_pip->Write();
  theta_pim->SetXTitle("theta_pim (Deg)");
  if (theta_pim->GetEntries()) theta_pim->Write();

  Phi_gamma->SetXTitle("Phi_gamma (Deg)");
  if (Phi_gamma->GetEntries()) Phi_gamma->Write();
  Phi_prot->SetXTitle("Phi_prot (Deg)");
  if (Phi_prot->GetEntries()) Phi_prot->Write();
  Phi_pip->SetXTitle("Phi_pip (Deg)");
  if (Phi_pip->GetEntries()) Phi_pip->Write();
  Phi_pim->SetXTitle("Phi_pim (Deg)");
  if (Phi_pim->GetEntries()) Phi_pim->Write();

  alpha_prot->SetXTitle("alpha_prot (Deg)");
  if (alpha_prot->GetEntries()) alpha_prot->Write();
  alpha_pip->SetXTitle("alpha_pip (Deg)");
  if (alpha_pip->GetEntries()) alpha_pip->Write();
  alpha_pim->SetXTitle("alpha_pim (Deg)");
  if (alpha_pim->GetEntries()) alpha_pim->Write();

  theta_prot_mc->SetXTitle("theta_prot_mc (Deg)");
  if (theta_prot_mc->GetEntries()) theta_prot_mc->Write();
  theta_pip_mc->SetXTitle("theta_pip_mc (Deg)");
  if (theta_pip_mc->GetEntries()) theta_pip_mc->Write();
  theta_pim_mc->SetXTitle("theta_pim_mc (Deg)");
  if (theta_pim_mc->GetEntries()) theta_pim_mc->Write();

  Phi_gamma_mc->SetXTitle("Phi_gamma_mc (Deg)");
  if (Phi_gamma_mc->GetEntries()) Phi_gamma_mc->Write();
  Phi_prot_mc->SetXTitle("Phi_prot_mc (Deg)");
  if (Phi_prot_mc->GetEntries()) Phi_prot_mc->Write();
  Phi_pip_mc->SetXTitle("Phi_pip_mc (Deg)");
  if (Phi_pip_mc->GetEntries()) Phi_pip_mc->Write();
  Phi_pim_mc->SetXTitle("Phi_pim_mc (Deg)");
  if (Phi_pim_mc->GetEntries()) Phi_pim_mc->Write();

  alpha_prot_mc->SetXTitle("alpha_prot_mc (Deg)");
  if (alpha_prot_mc->GetEntries()) alpha_prot_mc->Write();
  alpha_pip_mc->SetXTitle("alpha_pip_mc (Deg)");
  if (alpha_pip_mc->GetEntries()) alpha_pip_mc->Write();
  alpha_pim_mc->SetXTitle("alpha_pim_mc (Deg)");
  if (alpha_pim_mc->GetEntries()) alpha_pim_mc->Write();

  theta_prot_thrown->SetXTitle("theta_prot_thrown (Deg)");
  if (theta_prot_thrown->GetEntries()) theta_prot_thrown->Write();
  theta_pip_thrown->SetXTitle("theta_pip_thrown (Deg)");
  if (theta_pip_thrown->GetEntries()) theta_pip_thrown->Write();
  theta_pim_thrown->SetXTitle("theta_pim_thrown (Deg)");
  if (theta_pim_thrown->GetEntries()) theta_pim_thrown->Write();

  Phi_gamma_thrown->SetXTitle("Phi_gamma_thrown (Deg)");
  if (Phi_gamma_thrown->GetEntries()) Phi_gamma_thrown->Write();
  Phi_prot_thrown->SetXTitle("Phi_prot_thrown (Deg)");
  if (Phi_prot_thrown->GetEntries()) Phi_prot_thrown->Write();
  Phi_pip_thrown->SetXTitle("Phi_pip_thrown (Deg)");
  if (Phi_pip_thrown->GetEntries()) Phi_pip_thrown->Write();
  Phi_pim_thrown->SetXTitle("Phi_pim_thrown (Deg)");
  if (Phi_pim_thrown->GetEntries()) Phi_pim_thrown->Write();

  alpha_prot_thrown->SetXTitle("alpha_prot_thrown (Deg)");
  if (alpha_prot_thrown->GetEntries()) alpha_prot_thrown->Write();
  alpha_pip_thrown->SetXTitle("alpha_pip_thrown (Deg)");
  if (alpha_pip_thrown->GetEntries()) alpha_pip_thrown->Write();
  alpha_pim_thrown->SetXTitle("alpha_pim_thrown (Deg)");
  if (alpha_pim_thrown->GetEntries()) alpha_pim_thrown->Write();

  W_vs_q2_twoPi->SetXTitle("W (GeV)");
  W_vs_q2_twoPi->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2_twoPi->SetOption("COLZ1");
  if (W_vs_q2_twoPi->GetEntries()) W_vs_q2_twoPi->Write();

  W_hist_twoPi->SetXTitle("W (GeV)");
  if (W_hist_twoPi->GetEntries()) W_hist_twoPi->Write();

  Q2_hist_twoPi->SetXTitle("Q^{2} (GeV^{2})");
  if (Q2_hist_twoPi->GetEntries()) Q2_hist_twoPi->Write();

  MM_twoPi->SetXTitle("MM_twoPi (GeV)");
  if (MM_twoPi->GetEntries()) MM_twoPi->Write();
  MM2_twoPi->SetXTitle("MM2_twoPi (GeV)");
  if (MM2_twoPi->GetEntries()) MM2_twoPi->Write();
  
  W_q2_15_twoPi->SetXTitle("W (GeV)");
  W_q2_15_twoPi->SetTitle("electron production #sigma (Q^{2}<1.5 )");
  if (W_q2_15_twoPi->GetEntries()) W_q2_15_twoPi->Write();
  
  W_q2_2_twoPi->SetXTitle("W (GeV)");
   W_q2_2_twoPi->SetTitle("electron production #sigma (Q^{2}<2.0 )");
  if (W_q2_2_twoPi->GetEntries()) W_q2_2_twoPi->Write();
  
  W_q2_23_twoPi->SetXTitle("W (GeV)");
   W_q2_23_twoPi->SetTitle("electron production #sigma (2.0<Q^{2}<3.0 )");
  if (W_q2_23_twoPi->GetEntries()) W_q2_23_twoPi->Write();
  
  W_q2_34_twoPi->SetXTitle("W (GeV)"); 
  W_q2_34_twoPi->SetTitle("electroproduction #sigma (3.0<Q^{2}<4.0 )");
  if (W_q2_34_twoPi->GetEntries()) W_q2_34_twoPi->Write();
  
  W_q2_15_photon_twoPi->SetXTitle("W (GeV)");
  W_q2_15_photon_twoPi->SetTitle("virtual photo-production #sigma (Q^{2}<1.5 )");
  if (W_q2_15_photon_twoPi->GetEntries()) W_q2_15_photon_twoPi->Write();
  
  W_q2_2_photon_twoPi->SetXTitle("W (GeV)");
  W_q2_2_photon_twoPi->SetTitle("virtual photo-production #sigma (Q^{2}<2.0 )");
  if (W_q2_2_photon_twoPi->GetEntries()) W_q2_2_photon_twoPi->Write();
  
  W_q2_23_photon_twoPi->SetXTitle("W (GeV)");
  W_q2_23_photon_twoPi->SetTitle("virtual photo-production #sigma (2.0<Q^{2}<3.0 )");
  if (W_q2_23_photon_twoPi->GetEntries()) W_q2_23_photon_twoPi->Write();
  
  W_q2_34_photon_twoPi->SetXTitle("W (GeV)");
  W_q2_34_photon_twoPi->SetTitle("virtual photo-production #sigma (3.0<Q^{2}<4.0 )");
  if (W_q2_34_photon_twoPi->GetEntries()) W_q2_34_photon_twoPi->Write();
  //
  //
  //modified by Jiawei
  invmass_Ppip_twoPi->SetXTitle("inv_mass_Ppi+_twoPi");
  if (invmass_Ppip_twoPi->GetEntries()) invmass_Ppip_twoPi->Write();
  
  invmass_Ppim_twoPi->SetXTitle("inv_mass_Ppi-_twoPi");
  if (invmass_Ppim_twoPi->GetEntries()) invmass_Ppim_twoPi->Write();
  
  invmass_pippim_twoPi->SetXTitle("inv_mass_pi+pi-_twoPi");
  if (invmass_pippim_twoPi->GetEntries()) invmass_pippim_twoPi->Write();
 
  
  theta_Ppip_vs_invmass_twoPi->SetXTitle("invmass_P#pi^{+}_twoPi");
  theta_Ppip_vs_invmass_twoPi->SetYTitle("theta_P#pi^{+}_twoPi");
  theta_Ppip_vs_invmass_twoPi->SetOption("COLZ1");
  theta_Ppip_vs_invmass_twoPi->Write();
  
  theta_Ppim_vs_invmass_twoPi->SetXTitle("invmass_P#pi^{-}_twoPi");
  theta_Ppim_vs_invmass_twoPi->SetYTitle("theta_P#pi^{-}_twoPi");
  theta_Ppim_vs_invmass_twoPi->SetOption("COLZ1");
  theta_Ppim_vs_invmass_twoPi->Write();
  
  theta_pippim_vs_invmass_twoPi->SetXTitle("invmass_#pi^{+}#pi^{-}_twoPi");
  theta_pippim_vs_invmass_twoPi->SetYTitle("theta_#pi^{+}#pi^{-}_twoPi");
  theta_pippim_vs_invmass_twoPi->SetOption("COLZ1");
  theta_pippim_vs_invmass_twoPi->Write(); 

  
  p_vs_phi_Prot_twoPi->SetXTitle("Phi_proton_twoPi(deg)");
  p_vs_phi_Prot_twoPi->SetYTitle("P_proton_twoPi(GeV)");
  p_vs_phi_Prot_twoPi->SetOption("COLZ1");
  p_vs_phi_Prot_twoPi->Write();

  p_vs_phi_pip_twoPi->SetXTitle("Phi_#pi^{+}_twoPi(deg)");
  p_vs_phi_pip_twoPi->SetYTitle("P_#pi^{+}_twoPi(GeV)");
  p_vs_phi_pip_twoPi->SetOption("COLZ1");
  p_vs_phi_pip_twoPi->Write();

  p_vs_phi_pim_twoPi->SetXTitle("Phi_#pi^{-}_twoPi(deg)");
  p_vs_phi_pim_twoPi->SetYTitle("P_#pi^-}_twoPi(GeV)");
  p_vs_phi_pim_twoPi->SetOption("COLZ1");
  p_vs_phi_pim_twoPi->Write();  


  p_vs_theta_Prot_twoPi->SetXTitle("Theta_proton_twoPi(deg)");
  p_vs_theta_Prot_twoPi->SetYTitle("P_proton_twoPi(GeV)");
  p_vs_theta_Prot_twoPi->SetOption("COLZ1");
  p_vs_theta_Prot_twoPi->Write();

  p_vs_theta_pip_twoPi->SetXTitle("Theta_#pi^{+}_twoPi(deg)");
  p_vs_theta_pip_twoPi->SetYTitle("P_#pi^{+}_twoPi(GeV)");
  p_vs_theta_pip_twoPi->SetOption("COLZ1");
  p_vs_theta_pip_twoPi->Write();

  p_vs_theta_pim_twoPi->SetXTitle("Theta_#pi^{-}_twoPi(deg)");
  p_vs_theta_pim_twoPi->SetYTitle("P_#pi^{-}_twoPi(GeV)");
  p_vs_theta_pim_twoPi->SetOption("COLZ1");
  p_vs_theta_pim_twoPi->Write();


  theta_vs_phi_Prot_twoPi->SetXTitle("Phi_proton_twoPi(deg)");
  theta_vs_phi_Prot_twoPi->SetYTitle("Theta_proton_twoPi(deg)");
  theta_vs_phi_Prot_twoPi->SetOption("COLZ1");
  theta_vs_phi_Prot_twoPi->Write();

  theta_vs_phi_pip_twoPi->SetXTitle("Phi_#pi^{+}_twoPi(deg)");
  theta_vs_phi_pip_twoPi->SetYTitle("Theta_#pi^{+}_twoPi(deg)");
  theta_vs_phi_pip_twoPi->SetOption("COLZ1");
  theta_vs_phi_pip_twoPi->Write();

  theta_vs_phi_pim_twoPi->SetXTitle("Phi_#pi^{-}_twoPi(deg)");
  theta_vs_phi_pim_twoPi->SetYTitle("Theta_#pi^{-}_twoPi(deg)");
  theta_vs_phi_pim_twoPi->SetOption("COLZ1");
  theta_vs_phi_pim_twoPi->Write();
  
  //
  //
  //
 
  W_vs_q2_singlePi->SetXTitle("W (GeV)");
  W_vs_q2_singlePi->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2_singlePi->SetOption("COLZ1");
  if (W_vs_q2_singlePi->GetEntries()) W_vs_q2_singlePi->Write();

  W_hist_singlePi->SetXTitle("W (GeV)");
  if (W_hist_singlePi->GetEntries()) W_hist_singlePi->Write();

  Q2_hist_singlePi->SetXTitle("Q^{2} (GeV^{2})");
  if (Q2_hist_singlePi->GetEntries()) Q2_hist_singlePi->Write();

  if (MM_neutron->GetEntries()) MM_neutron->Write();

  auto wvsq2_sec = RootOutputFile->mkdir("wvsq2_sec");
  wvsq2_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_sec[i]->SetOption("COLZ1");
    W_vs_q2_sec[i]->Write();
  }
  auto w_sec = RootOutputFile->mkdir("w_sec");
  w_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_sec[i]->SetXTitle("W (GeV)");

    W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.85, 1.05);
    gStyle->SetOptFit(01);
    W_sec[i]->Write();
  }
  Theta_vs_mom_x_mu->SetYTitle("Theta_x_mu");
  Theta_vs_mom_x_mu->SetXTitle("mom_x_mu");
  Theta_vs_mom_x_mu->SetOption("COLZ1");
  Theta_vs_mom_x_mu->Write();

  auto singlePi_sec = RootOutputFile->mkdir("singlePi_sec");
  singlePi_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_singlePi_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_singlePi_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_singlePi_sec[i]->SetOption("COLZ1");
    if (W_vs_q2_singlePi_sec[i]->GetEntries()) W_vs_q2_singlePi_sec[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_vs_MM_singlePi[i]->SetOption("COLZ1");
    W_vs_MM_singlePi[i]->SetYTitle("MM (GeV)");
    W_vs_MM_singlePi[i]->SetXTitle("W (GeV)");
    if (W_vs_MM_singlePi[i]->GetEntries()) W_vs_MM_singlePi[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_singlePi_sec[i]->SetXTitle("W (GeV)");
    if (W_singlePi_sec[i]->GetEntries()) W_singlePi_sec[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    //    if (MM_neutron_sec[i]->GetEntries()) MM_neutron_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.7, 1.1);
    MM_neutron_sec[i]->SetXTitle("Mass (GeV)");
    if (MM_neutron_sec[i]->GetEntries()) MM_neutron_sec[i]->Write();
  }

  auto twoPi_sec = RootOutputFile->mkdir("twoPi_sec");
  twoPi_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_twoPi_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_twoPi_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_twoPi_sec[i]->SetOption("COLZ1");
    if (W_vs_q2_twoPi_sec[i]->GetEntries()) W_vs_q2_twoPi_sec[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_vs_MM_twoPi[i]->SetOption("COLZ1");
    W_vs_MM_twoPi[i]->SetYTitle("MM (GeV)");
    W_vs_MM_twoPi[i]->SetXTitle("W (GeV)");
    if (W_vs_MM_twoPi[i]->GetEntries()) W_vs_MM_twoPi[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_twoPi_sec[i]->SetXTitle("W (GeV)");
    if (W_twoPi_sec[i]->GetEntries()) W_twoPi_sec[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    //  if (MM_twoPi_sec[i]->GetEntries()) MM_twoPi_sec[i]->Fit("gaus", "QMR+", "QMR+", -0.1, 0.1);
    MM_twoPi_sec[i]->SetXTitle("MisingMass (GeV)");
    MM_twoPi_sec[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    // if (MM2_twoPi_sec[i]->GetEntries()) MM2_twoPi_sec[i]->Fit("gaus", "QMR+", "QMR+", -0.1, 0.1);
    MM2_twoPi_sec[i]->SetXTitle("MisingMass (GeV)");
    MM2_twoPi_sec[i]->Write();
  }

  auto Inv_Mass_hist = RootOutputFile->mkdir("Inv_Mass_hist");
  Inv_Mass_hist->cd();
  for (short i = 0; i < w_range_num; i++) {
    for (short j = 0; j < q2_range_num; j++) {
      inv_mass_P_pip[i][j]->SetXTitle("inv_mass (GeV)");
      inv_mass_P_pip[i][j]->Write();
      inv_mass_P_pim[i][j]->SetXTitle("inv_mass (GeV)");
      inv_mass_P_pim[i][j]->Write();
      inv_mass_pip_pim[i][j]->SetXTitle("inv_mass (GeV)");
      inv_mass_pip_pim[i][j]->Write();

      theta_P_vs_mass_pip_pim[i][j]->SetYTitle("theta_P_cm");
      theta_P_vs_mass_pip_pim[i][j]->SetXTitle("inv_mass_pip_pim (GeV)");
      theta_P_vs_mass_pip_pim[i][j]->SetOption("COLZ1");
      theta_P_vs_mass_pip_pim[i][j]->Write();

      theta_pip_vs_mass_Ppim[i][j]->SetYTitle("theta_pip_cm");
      theta_pip_vs_mass_Ppim[i][j]->SetXTitle("inv_mass_Ppim (GeV)");
      theta_pip_vs_mass_Ppim[i][j]->SetOption("COLZ1");
      theta_pip_vs_mass_Ppim[i][j]->Write();

      theta_pim_vs_mass_Ppip[i][j]->SetYTitle("theta_pim_cm");
      theta_pim_vs_mass_Ppip[i][j]->SetXTitle("inv_mass_Ppip (GeV)");
      theta_pim_vs_mass_Ppip[i][j]->SetOption("COLZ1");
      theta_pim_vs_mass_Ppip[i][j]->Write();

      theta_P_lab_vs_mass_pip_pim[i][j]->SetYTitle("theta_P_cm");
      theta_P_lab_vs_mass_pip_pim[i][j]->SetXTitle("inv_mass_pip_pim (GeV)");
      theta_P_lab_vs_mass_pip_pim[i][j]->SetOption("COLZ1");
      theta_P_lab_vs_mass_pip_pim[i][j]->Write();

      theta_pip_lab_vs_mass_Ppim[i][j]->SetYTitle("theta_pip_cm");
      theta_pip_lab_vs_mass_Ppim[i][j]->SetXTitle("inv_mass_Ppim (GeV)");
      theta_pip_lab_vs_mass_Ppim[i][j]->SetOption("COLZ1");
      theta_pip_lab_vs_mass_Ppim[i][j]->Write();

      theta_pim_lab_vs_mass_Ppip[i][j]->SetYTitle("theta_pim_cm");
      theta_pim_lab_vs_mass_Ppip[i][j]->SetXTitle("inv_mass_Ppip (GeV)");
      theta_pim_lab_vs_mass_Ppip[i][j]->SetOption("COLZ1");
      theta_pim_lab_vs_mass_Ppip[i][j]->Write();
    }
  }
  auto twoPi_sec_thrown = RootOutputFile->mkdir("twoPi_sec_thrown");
  twoPi_sec_thrown->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_twoPi_sec_thrown[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_twoPi_sec_thrown[i]->SetXTitle("W (GeV)");
    W_vs_q2_twoPi_sec_thrown[i]->SetOption("COLZ1");
    if (W_vs_q2_twoPi_sec_thrown[i]->GetEntries()) W_vs_q2_twoPi_sec_thrown[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_vs_MM_twoPi_thrown[i]->SetOption("COLZ1");
    W_vs_MM_twoPi_thrown[i]->SetYTitle("MM (GeV)");
    W_vs_MM_twoPi_thrown[i]->SetXTitle("W (GeV)");
    if (W_vs_MM_twoPi[i]->GetEntries()) W_vs_MM_twoPi_thrown[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_twoPi_sec_thrown[i]->SetXTitle("W (GeV)");
    if (W_twoPi_sec_thrown[i]->GetEntries()) W_twoPi_sec_thrown[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    //  if (MM_twoPi_sec_thrown[i]->GetEntries()) MM_twoPi_sec_thrown[i]->Fit("gaus", "QMR+", "QMR+", -0.1, 0.1);
    MM_twoPi_sec_thrown[i]->SetXTitle("MisingMass (GeV)");
    MM_twoPi_sec_thrown[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    // if (MM2_twoPi_sec_thrown[i]->GetEntries()) MM2_twoPi_sec_thrown[i]->Fit("gaus", "QMR+", "QMR+", -0.1, 0.1);
    MM2_twoPi_sec_thrown[i]->SetXTitle("MisingMass (GeV)");
    MM2_twoPi_sec_thrown[i]->Write();
  }

  auto Inv_Mass_hist_thrown = RootOutputFile->mkdir("Inv_Mass_hist_thrown");
  Inv_Mass_hist_thrown->cd();
  for (short i = 0; i < w_range_num; i++) {
    for (short j = 0; j < q2_range_num; j++) {
      inv_mass_P_pip_thrown[i][j]->SetXTitle("inv_mass (GeV)");
      inv_mass_P_pip_thrown[i][j]->Write();
      inv_mass_P_pim_thrown[i][j]->SetXTitle("inv_mass (GeV)");
      inv_mass_P_pim_thrown[i][j]->Write();
      inv_mass_pip_pim_thrown[i][j]->SetXTitle("inv_mass (GeV)");
      inv_mass_pip_pim_thrown[i][j]->Write();

      theta_P_vs_mass_pip_pim_thrown[i][j]->SetYTitle("theta_P_cm");
      theta_P_vs_mass_pip_pim_thrown[i][j]->SetXTitle("inv_mass_pip_pim (GeV)");
      theta_P_vs_mass_pip_pim_thrown[i][j]->SetOption("COLZ1");
      theta_P_vs_mass_pip_pim_thrown[i][j]->Write();

      theta_pip_vs_mass_Ppim_thrown[i][j]->SetYTitle("theta_pip_cm");
      theta_pip_vs_mass_Ppim_thrown[i][j]->SetXTitle("inv_mass_Ppim (GeV)");
      theta_pip_vs_mass_Ppim_thrown[i][j]->SetOption("COLZ1");
      theta_pip_vs_mass_Ppim_thrown[i][j]->Write();

      theta_pim_vs_mass_Ppip_thrown[i][j]->SetYTitle("theta_pim_cm");
      theta_pim_vs_mass_Ppip_thrown[i][j]->SetXTitle("inv_mass_Ppip (GeV)");
      theta_pim_vs_mass_Ppip_thrown[i][j]->SetOption("COLZ1");
      theta_pim_vs_mass_Ppip_thrown[i][j]->Write();

      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetYTitle("theta_P_cm");
      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetXTitle("inv_mass_pip_pim (GeV)");
      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetOption("COLZ1");
      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->Write();

      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetYTitle("theta_pip_cm");
      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetXTitle("inv_mass_Ppim (GeV)");
      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetOption("COLZ1");
      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->Write();

      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetYTitle("theta_pim_cm");
      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetXTitle("inv_mass_Ppip (GeV)");
      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetOption("COLZ1");
      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->Write();
    }
  }
  auto Npip_sec = RootOutputFile->mkdir("Npip_sec");
  Npip_sec->cd();

  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_Npip_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_Npip_sec[i]->SetXTitle("W (GeV)");
    W_vs_q2_Npip_sec[i]->SetOption("COLZ1");
    W_vs_q2_Npip_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    W_Npip_sec[i]->SetXTitle("W (GeV)");
    W_Npip_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    MM_Npip_sec[i]->SetXTitle("Mass (GeV)");
    MM_Npip_sec[i]->Write();
  }
}
void Histogram::makeHists_electron_cuts() {
  for (auto&& cut : before_after_cut) {
    int c = cut.first;
    auto type = cut.second.c_str();
    EC_sampling_fraction[c] =
        std::make_shared<TH2D>(Form("EC_sampling_fraction%s", type), Form("EC_sampling_fraction%s", type), bins, p_min,
                               p_max, bins, zero, 0.5);
    vz_position[c] = std::make_shared<TH1D>(Form("vz_position%s", type), Form("vz_position%s", type), bins, -40, 40);
    pcal_sec[c] =
        std::make_shared<TH2D>(Form("pcal_sec%s", type), Form("pcal_sec%s", type), bins, -420, 420, bins, -420, 420);
    dcr1_sec[c] =
        std::make_shared<TH2D>(Form("dcr1_sec%s", type), Form("dcr1_sec%s", type), bins, -180, 180, bins, -180, 180);
    dcr2_sec[c] =
        std::make_shared<TH2D>(Form("dcr2_sec%s", type), Form("dcr2_sec%s", type), bins, -270, 270, bins, -270, 270);
    dcr3_sec[c] =
        std::make_shared<TH2D>(Form("dcr3_sec%s", type), Form("dcr3_sec%s", type), bins, -320, 320, bins, -320, 320);
  }
}
void Histogram::FillHists_electron_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[before_cut]->Fill(_d->vz(0));
  pcal_sec[before_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[before_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[before_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[before_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[before_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}

void Histogram::FillHists_electron_with_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[after_cut]->Fill(_d->vz(0));
  pcal_sec[after_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[after_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[after_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[after_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[after_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}

void Histogram::Write_Electron_cuts() {
  for (auto&& cut : before_after_cut) {
    int c = cut.first;
    vz_position[c]->SetXTitle("vz (GeV)");
    vz_position[c]->Fit("gaus", "QMR+", "QMR+", -7.089, 2.0);
    gStyle->SetOptFit(1111);
    if (vz_position[c]->GetEntries()) vz_position[c]->Write();
    pcal_sec[c]->SetXTitle("x (cm)");
    pcal_sec[c]->SetYTitle("y (cm)");
    pcal_sec[c]->SetOption("COLZ1");
    if (pcal_sec[c]->GetEntries()) pcal_sec[c]->Write();

    dcr1_sec[c]->SetXTitle("x (cm)");
    dcr1_sec[c]->SetYTitle("y (cm)");
    dcr1_sec[c]->SetOption("COLZ1");
    if (dcr1_sec[c]->GetEntries()) dcr1_sec[c]->Write();

    dcr2_sec[c]->SetXTitle("x (cm)");
    dcr2_sec[c]->SetYTitle("y (cm)");
    dcr2_sec[c]->SetOption("COLZ1");
    if (dcr2_sec[c]->GetEntries()) dcr2_sec[c]->Write();

    dcr3_sec[c]->SetXTitle("x (cm)");
    dcr3_sec[c]->SetYTitle("y (cm)");
    dcr3_sec[c]->SetOption("COLZ1");
    if (dcr3_sec[c]->GetEntries()) dcr3_sec[c]->Write();

    EC_sampling_fraction[c]->SetXTitle("Momentum (GeV)");
    EC_sampling_fraction[c]->SetYTitle("Sampling Fraction");
    EC_sampling_fraction[c]->SetOption("COLZ1");
    EC_sampling_fraction[c]->Write();
  }
}
void Histogram::makeHists_sector() {
  for (short i = 0; i < w_range_num; i++) {
    for (short j = 0; j < q2_range_num; j++) {
      inv_mass_P_pip[i][j] = std::make_shared<TH1D>(
          Form("inv_mass_P_pip_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("inv_mass_P_pip: %6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 1., 3.5);
      inv_mass_P_pim[i][j] = std::make_shared<TH1D>(
          Form("inv_mass_P_pim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("inv_mass_P_pim: %6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 1., 3.5);
      inv_mass_pip_pim[i][j] = std::make_shared<TH1D>(
          Form("inv_mass_pip_pim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("inv_mass_pip_pim: %6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 0.2, 2.5);
      theta_P_vs_mass_pip_pim[i][j] = std::make_shared<TH2D>(
          Form("theta_P_vs_inv_mass_pip_pim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_P_vs_inv_mass_pip_pim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins,
          0.2, 3.5, bins, zero, 180);
      theta_pim_vs_mass_Ppip[i][j] = std::make_shared<TH2D>(
          Form("theta_pim_vs_inv_mass_Ppip_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_pim_vs_inv_mass_Ppip_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 1.,
          3.5, bins, zero, 180);
      theta_pip_vs_mass_Ppim[i][j] = std::make_shared<TH2D>(
          Form("theta_pip_vs_inv_mass_Ppim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_pip_vs_inv_mass_Ppim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 1.,
          3.5, bins, zero, 180);
      theta_P_lab_vs_mass_pip_pim[i][j] = std::make_shared<TH2D>(
          Form("theta_P_lab_vs_inv_mass_pip_pim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_P_lab_vs_inv_mass_pip_pim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          bins, 0.2, 2.5, bins, zero, 60);
      theta_pim_lab_vs_mass_Ppip[i][j] = std::make_shared<TH2D>(
          Form("theta_pim_lab_vs_inv_mass_Ppip_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_pim_lab_vs_inv_mass_Ppip_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins,
          1., 3.5, bins, zero, 120);
      theta_pip_lab_vs_mass_Ppim[i][j] = std::make_shared<TH2D>(
          Form("theta_pip_lab_vs_inv_mass_Ppim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_pip_lab_vs_inv_mass_Ppim_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins,
          1., 3.5, bins, zero, 120);

      inv_mass_P_pip_thrown[i][j] = std::make_shared<TH1D>(
          Form("inv_mass_P_pip_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("inv_mass_P_pip_thrown_: %6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 1.,
          3.5);
      inv_mass_P_pim_thrown[i][j] = std::make_shared<TH1D>(
          Form("inv_mass_P_pim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("inv_mass_P_pim_thrown_: %6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 1.,
          3.5);
      inv_mass_pip_pim_thrown[i][j] = std::make_shared<TH1D>(
          Form("inv_mass_pip_pim_thrown%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("inv_mass_pip_pim_thrown: %6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()), bins, 0.2,
          2.5);
      theta_P_vs_mass_pip_pim_thrown[i][j] = std::make_shared<TH2D>(
          Form("theta_P_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_P_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          bins, 0.2, 3.5, bins, zero, 180);
      theta_pim_vs_mass_Ppip_thrown[i][j] = std::make_shared<TH2D>(
          Form("theta_pim_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_pim_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          bins, 1., 3.5, bins, zero, 180);
      theta_pip_vs_mass_Ppim_thrown[i][j] = std::make_shared<TH2D>(
          Form("theta_pip_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          Form("theta_pip_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(), q2_range_name[j].c_str()),
          bins, 1., 3.5, bins, zero, 180);
      theta_P_lab_vs_mass_pip_pim_thrown[i][j] =
          std::make_shared<TH2D>(Form("theta_P_lab_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                      q2_range_name[j].c_str()),
                                 Form("theta_P_lab_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                      q2_range_name[j].c_str()),
                                 bins, 0.2, 2.5, bins, zero, 70);
      theta_pim_lab_vs_mass_Ppip_thrown[i][j] =
          std::make_shared<TH2D>(Form("theta_pim_lab_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                      q2_range_name[j].c_str()),
                                 Form("theta_pim_lab_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                      q2_range_name[j].c_str()),
                                 bins, 1., 3.5, bins, zero, 120);
      theta_pip_lab_vs_mass_Ppim_thrown[i][j] =
          std::make_shared<TH2D>(Form("theta_pip_lab_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                      q2_range_name[j].c_str()),
                                 Form("theta_pip_lab_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                      q2_range_name[j].c_str()),
                                 bins, 1., 3.5, bins, zero, 120);
    }
  }
  for (short i = 0; i < 3; i++) {
    W_det[i] = std::make_shared<TH1D>(Form("W_det_%d", i + 1), Form("W detector: %d", i + 1), bins, zero, w_max);
    if (i == 0)
      WQ2_det[i] = std::make_shared<TH2D>(Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1), bins, zero,
                                          w_max, bins, zero, 0.5);
    else
      WQ2_det[i] = std::make_shared<TH2D>(Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1), bins, zero,
                                          w_max, bins, zero, q2_max);
  }

  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i] = std::make_shared<TH2D>(Form("wvsq2_sec_%d", i + 1), Form("W vs Q^{2} Sector: %d", i + 1), bins,
                                            zero, w_max, bins, zero, q2_max);

    W_sec[i] = std::make_shared<TH1D>(Form("w_sec_%d", i + 1), Form("W Sector: %d", i + 1), bins, zero, w_max);

    W_vs_q2_twoPi_sec[i] =
        std::make_shared<TH2D>(Form("wvsq2_sec_twoPi_%d", i + 1), Form("W vs Q^{2} W_twoPi Sector: %d", i + 1), bins,
                               zero, w_max, bins, zero, q2_max);
    W_twoPi_sec[i] =
        std::make_shared<TH1D>(Form("w_sec_twoPi_%d", i + 1), Form("W twoPi Sector: %d", i + 1), bins, zero, w_max);

    W_vs_MM_twoPi[i] = std::make_shared<TH2D>(Form("W_vs_MM_twoPi_%d", i + 1), Form("W_vs_MM_twoPi_%d", i + 1), bins,
                                              zero, w_max, bins, -q2_max, q2_max);

    MM_twoPi_sec[i] =
        std::make_shared<TH1D>(Form("MM_sec_%d", i + 1), Form("MM twoPi Sector: %d", i + 1), bins, -0.4, 0.4);
    MM2_twoPi_sec[i] =
        std::make_shared<TH1D>(Form("MM_SQ_sec_%d", i + 1), Form("MM_SQ twoPi Sector: %d", i + 1), bins, -0.2, 0.2);

    W_vs_MM_twoPi_thrown[i] =
        std::make_shared<TH2D>(Form("W_vs_MM_twoPi_thrown_%d", i + 1), Form("W_vs_MM_twoPi_thrown_%d", i + 1), bins,
                               zero, w_max, bins, -q2_max, q2_max);

    W_vs_q2_twoPi_sec_thrown[i] = std::make_shared<TH2D>(Form("wvsq2_sec_twoPi_thrown_%d", i + 1),
                                                         Form("W vs Q^{2} W_twoPi Sector_thrown: %d", i + 1), bins,
                                                         zero, w_max, bins, zero, q2_max);
    W_twoPi_sec_thrown[i] = std::make_shared<TH1D>(Form("w_sec_twoPi_thrown_%d", i + 1),
                                                   Form("W twoPi Sector_thrown: %d", i + 1), bins, zero, w_max);

    MM_twoPi_sec_thrown[i] = std::make_shared<TH1D>(Form("MM_sec_thrown_%d", i + 1),
                                                    Form("MM twoPi Sector_thrown: %d", i + 1), bins, -0.4, 0.4);
    MM2_twoPi_sec_thrown[i] = std::make_shared<TH1D>(Form("MM_SQ_sec_thrown_%d", i + 1),
                                                     Form("MM_SQ twoPi Sector_thrown: %d", i + 1), bins, -0.2, 0.2);

    W_vs_q2_singlePi_sec[i] =
        std::make_shared<TH2D>(Form("wvsq2_sec_singlePi_%d", i + 1), Form("W vs Q^{2} W_singlePi Sector: %d", i + 1),
                               bins, zero, w_max, bins, zero, q2_max);

    W_singlePi_sec[i] = std::make_shared<TH1D>(Form("w_sec_singlePi_%d", i + 1), Form("W singlePi Sector: %d", i + 1),
                                               bins, zero, w_max);

    W_vs_q2_Npip_sec[i] =
        std::make_shared<TH2D>(Form("wvsq2_sec_Npip_%d", i + 1), Form("W vs Q^{2} W_Npip Sector: %d", i + 1), bins,
                               zero, w_max, bins, zero, q2_max);

    W_Npip_sec[i] =
        std::make_shared<TH1D>(Form("w_sec_Npip_%d", i + 1), Form("W Npip Sector: %d", i + 1), bins, zero, w_max);

    MM_neutron_sec[i] =
        std::make_shared<TH1D>(Form("MM_Sec_%d", i + 1), Form("MM neutron Sector: %d", i + 1), bins, zero, 4.0);

    MM_Npip_sec[i] = std::make_shared<TH1D>(Form("MM_Npip_Sec_%d", i + 1), Form("MM^{2} neutron pip Sector: %d", i + 1),
                                            bins, -2.0, 2.0);

    W_vs_MM_singlePi[i] = std::make_shared<TH2D>(Form("W_vs_MM_singlePi_%d", i + 1), Form("W_vs_MM_singlePi_%d", i + 1),
                                                 bins, zero, w_max, bins, -q2_max, q2_max);
  }
}

void Histogram::makeHists_deltat() {
  std::string tof = "";
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        tof = "ftof";
        delta_t_hist[p][c][i][0] =
            std::make_shared<TH2D>(Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   bins, p_min, p_max, bins, Dt_min, Dt_max);

        tof = "ctof";
        delta_t_hist[p][c][i][1] =
            std::make_shared<TH2D>(Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   bins, 0, 3.0, bins, -6.0, 6.0);
      }
    }
  }
}

void Histogram::Fill_deltat_pi(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt, int part) {
  auto _cuts = std::make_unique<Cuts>(data, dt);
  int charge = data->charge(part);
  bool fc = dt->ctof();
  int pid = data->pid(part);
  float mom = data->p(part);
  float time = NAN;
  if (fc)
    time = dt->dt_ctof_Pi();
  else
    time = dt->dt_Pi();

  if (charge == 1) {
    delta_t_hist[1][0][0][fc]->Fill(mom, time);
    if (_cuts->IsPip(part))
      delta_t_hist[1][0][1][fc]->Fill(mom, time);
    else
      delta_t_hist[1][0][2][fc]->Fill(mom, time);
  } else if (charge == -1) {
    delta_t_hist[1][1][0][fc]->Fill(mom, time);
    if (_cuts->IsPim(part))
      delta_t_hist[1][1][1][fc]->Fill(mom, time);
    else
      delta_t_hist[1][1][2][fc]->Fill(mom, time);
  }
}

void Histogram::Fill_deltat_prot(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt,
                                 int part) {
  auto _cuts = std::make_unique<Cuts>(data, dt);
  int charge = data->charge(part);
  bool fc = dt->ctof();
  int pid = data->pid(part);
  float mom = data->p(part);
  float time = NAN;

  if (fc)
    time = dt->dt_ctof_P();
  else
    time = dt->dt_P();

  if (charge == 1) {
    delta_t_hist[2][0][0][fc]->Fill(mom, time);
    if (_cuts->IsProton(part))
      delta_t_hist[2][0][1][fc]->Fill(mom, time);
    else
      delta_t_hist[2][0][2][fc]->Fill(mom, time);

    delta_t_hist[2][1][0][fc]->Fill(mom, time);
    if (pid == PROTON)
      delta_t_hist[2][1][1][fc]->Fill(mom, time);
    else
      delta_t_hist[2][1][2][fc]->Fill(mom, time);
  }
}

// do
// pid == 11 at first;
// skim vs pid at 0 is 11. compare.
void Histogram::Write_deltat() {
  TDirectory* ftof_folder = RootOutputFile->mkdir("ftof");
  ftof_folder->cd();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i][0]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i][0]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i][0]->SetOption("COLZ1");
        if (delta_t_hist[p][c][i][0]->GetEntries() > 1) delta_t_hist[p][c][i][0]->Write();
      }
    }
  }
  TDirectory* ctof_folder = RootOutputFile->mkdir("ctof");
  ctof_folder->cd();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i][1]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i][1]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i][1]->SetOption("COLZ1");
        if (delta_t_hist[p][c][i][1]->GetEntries() > 1) delta_t_hist[p][c][i][1]->Write();
      }
    }
  }
}

void Histogram::makeHists_MomVsBeta() {
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        momvsbeta_hist[p][c][i] = std::make_shared<TH2D>(
            Form("mom_vs_beta_%s_%s_%s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
            Form("Momentum vs #beta %s %s %s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
            bins, p_min, p_max, bins, zero, 1.2);
      }
    }
  }
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Branches12>& data, int part) {
  int good_ID = 0;
  float beta = data->beta(part);
  float mom = data->p(part);
  int charge = data->charge(part);
  int pid = data->pid(part);
  if (beta != 0) {
    momentum->Fill(mom);
    for (short p = 0; p < particle_num; p++) {
      switch (p) {
        case 0:
          good_ID = ELECTRON;
          break;
        case 1:
          good_ID = PIP;
          break;
        case 2:
          good_ID = PROTON;
          break;
        case 3:
          good_ID = KP;
          break;
      }

      momvsbeta_hist[p][0][0]->Fill(mom, beta);
      if (good_ID == abs(pid)) {
        momvsbeta_hist[p][0][1]->Fill(mom, beta);
      } else {
        momvsbeta_hist[p][0][2]->Fill(mom, beta);
      }

      if (charge == -1) {
        momvsbeta_hist[p][2][0]->Fill(mom, beta);
        if (-good_ID == pid) {
          momvsbeta_hist[p][2][1]->Fill(mom, beta);
        } else {
          momvsbeta_hist[p][2][2]->Fill(mom, beta);
        }
      } else if (charge == 1) {
        momvsbeta_hist[p][1][0]->Fill(mom, beta);
        if (good_ID == pid) {
          momvsbeta_hist[p][1][1]->Fill(mom, beta);
        } else {
          momvsbeta_hist[p][1][2]->Fill(mom, beta);
        }
      }
    }
  }
}

void Histogram::Write_MomVsBeta() {
  momentum->SetXTitle("Momentum (GeV)");
  momentum->Write();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        momvsbeta_hist[p][c][i]->SetYTitle("#beta");
        momvsbeta_hist[p][c][i]->SetOption("COLZ1");
        momvsbeta_hist[p][c][i]->Write();
      }
    }
  }
}
