#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <Rtypes.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TStopwatch.h>

// Define max number of particles in the event
#define MAX_TRACKS 15000
#define HARMONICNUMBER 2
#define WEIGHT 1

Float_t Qy(Float_t phi) {
	return WEIGHT * sin(HARMONICNUMBER * phi);
}

Float_t Qx(Float_t phi) {
	return WEIGHT * cos(HARMONICNUMBER * phi);
}

Float_t PHI(Float_t momy, Float_t momx) {
	return TMath::ATan2(momy, momx);
}

Float_t PSI_EP(Float_t Qx, Float_t Qy) {
	return TMath::ATan2(Qy, Qx) / HARMONICNUMBER;
}

Float_t vnobs(std::vector<Float_t> angles, Float_t ep_angle) {

	TH1F *htmp = new TH1F("h1", "h1 title", 100, -10, 10);
	Float_t tmp;
	Float_t mean = 0.0;

	for (Int_t i = 0; i < angles.size(); i++) {
			tmp = cos(HARMONICNUMBER * (angles[i] - ep_angle));
			htmp->Fill(tmp);
		}
		mean = htmp->GetMean();
		delete htmp;
		return mean;
}


std::vector<Float_t> vec_Cos_N_Phi_Psi(std::vector<Float_t> angles, Float_t ep_angle) {

	for (Int_t i = 0; i < angles.size(); i++) {
			angles[i] = cos(HARMONICNUMBER * (angles[i] - ep_angle));
		}

		return angles;
}


void readTree(TString inListName, TString outFileName) {
		// Start timer
		TStopwatch timer;
		timer.Start();

		// Read input from the file list
		std::ifstream file(inListName.Data());
		std::string line;
		TChain *chain = new TChain("mctree");
		Int_t nFiles = 0;
		while (std::getline(file, line)) {
			chain->Add(line.c_str());
			nFiles++;
		}

		// Init variables to read to
		Float_t         bimp; // - impact parameter
  	Float_t         phi2;
  	Float_t         phi3;
  	Float_t         ecc2;
  	Float_t         ecc3;
  	Int_t           npart;
  	Int_t           nh; // - number of particles in the event
  	Float_t         momx[MAX_TRACKS];   //[nh] - px
  	Float_t         momy[MAX_TRACKS];   //[nh] - py
  	Float_t         momz[MAX_TRACKS];   //[nh] - pz
  	Float_t         ene[MAX_TRACKS];   //[nh] - Energy
  	Int_t           hid[MAX_TRACKS];   //[nh]
  	Int_t           pdg[MAX_TRACKS];   //[nh] - PDG code (211 - pion, 321 - kaon, 2212 - proton)
  	Short_t         charge[MAX_TRACKS];   //[nh]	- charge

		chain->SetBranchAddress("bimp", &bimp);
  	chain->SetBranchAddress("phi2", &phi2);
  	chain->SetBranchAddress("phi3", &phi3);
  	chain->SetBranchAddress("ecc2", &ecc2);
  	chain->SetBranchAddress("ecc3", &ecc3);
  	chain->SetBranchAddress("npart", &npart);
  	chain->SetBranchAddress("nh", &nh);
  	chain->SetBranchAddress("momx", momx);
  	chain->SetBranchAddress("momy", momy);
  	chain->SetBranchAddress("momz", momz);
 		chain->SetBranchAddress("ene", ene);
  	chain->SetBranchAddress("hid", hid);
  	chain->SetBranchAddress("pdg", pdg);
  	chain->SetBranchAddress("charge", charge);

		 // Manage output
		 TFile *fo = new TFile(outFileName.Data(), "recreate");

		 /*Init my own variables*/

		 Int_t nbins = 10;

		 Float_t bimpValues[11] = {
       	 0.00, 4.17,
       	 6.03, 7.40,
       	 8.56, 9.63,
       	 10.59, 11.46,
       	 12.29, 13.30,
       	 14.94
			 };


    	Float_t Qx_tpcE = 0.0, Qx_tpcW = 0.0;
    	Float_t Qy_tpcE = 0.0, Qy_tpcW = 0.0;
    	Float_t phi = 0.0;
    	// Float_t qW = 0.0;
    	Float_t Psi_EP_tpcW = 0.0;
    	Float_t Psi_EP_tpcE = 0.0;

    	Float_t error = 0.0;
    	Float_t etta = 0.0;
    	Float_t momModule = 0.0;
    	Float_t momTransverse = 0.0;
    	Float_t pre_Resolution = 0.0;
			Float_t Resolution = 0.0;
			Float_t vnobs_tpcW = 0.0;
			Float_t vnobs_tpcE = 0.0;

    	Bool_t flag = false;

			char name[20];
    	char title[20];

			/* Add parameters od mom transverse to create h[bin][npTbins]*/

			Int_t npTbins = 10; // 0.5 - 3.6 GeV/c - number of npTbins bins
			Float_t bin_w[11] = {
				0.2, 0.4,
				0.6, 0.8,
				1.0, 1.2,
				1.4, 1.8,
				2.3, 2.8,
				4.0
			};

			Float_t maxpT = 4.0; // max pT
			Float_t minpT = 0.2; // min pT


			TH1F *hRes[nbins];
			TH1F *hFlowtpcE[nbins][npTbins];
			TH1F *hFlowtpcW[nbins][npTbins];

			Float_t Psi_EP_tpcEmat[nbins][npTbins];
			Float_t Psi_EP_tpcWmat[nbins][npTbins];

			std::vector<Float_t> tmpE;
			std::vector<Float_t> tmpW;

			std::vector<Float_t> tpcE_phi[nbins][npTbins];
			std::vector<Float_t> tpcW_phi[nbins][npTbins];

			/* End of init */

			// Init output histograms, profile, etc.

			// Reading events

			Long64_t nEvents = chain->GetEntries();
			std::cout << "Read " << nFiles << " files. " << nEvents << " events total." << std::endl;

			for (Int_t j = 0; j < npTbins; j++) {
				for ( Int_t i = 0; i < nbins; i++ ) {
					hFlowtpcE[i][j] = new TH1F(Form("tpcW: nBin: %d, npTbin: %d", i + 1, j + 1),
																 Form("v_{2}tpcW[%d][%d]", i + 1, j + 1),
																 400, -2.0, 2.0);

				  hFlowtpcW[i][j] = new TH1F(Form("tpcE: nBin: %d, npTbin: %d", i + 1, j + 1),
																 Form("v_{2}tpcE[%d][%d]", i + 1, j + 1),
																 400, -2.0, 2.0);
				}
			}

			for (Int_t i = 0; i < nbins; i++) {
				sprintf(name, "nBin: %d", i + 1); // IP - impact parameter
				sprintf(title, "Resolution, %d", i + 1);
				hRes[i] = new TH1F(name, title, 200, 0, 10);
			}


			// Collisions
			for (Long64_t ievent = 0; ievent < nEvents; ievent++) {

				if (ievent % 1000 == 0) std::cout << "Event [" << ievent << "/" << nEvents << "]" << std::endl;
				chain->GetEntry(ievent);

				// Do event-wise stuff (fill histograms, etc.)
				for (Int_t j = 0; j < npTbins; j++) {

							for (Int_t i = 0; i < nbins; i++) {

								if(bimpValues[i] < bimp &&  bimp < bimpValues[i + 1]) {

									// Reading particles from one collision
									for (int itrack = 0; itrack < nh; itrack++) {

										// Do particle-wise stuff (fill histograms, etc.)

										/* Begin of conditions */
										if (charge[itrack] != 0) {

											    momModule = sqrt(pow(momx[itrack], 2) + pow(momy[itrack], 2) + pow(momz[itrack], 2));
											    etta = 0.5 * log((momModule + momz[itrack]) / (momModule - momz[itrack]));
											    momTransverse = sqrt(pow(momx[itrack], 2) + pow(momy[itrack], 2));

													if (bin_w[j] < momTransverse &&  momTransverse < bin_w[j + 1]) {

																// TPC(E)
														    if (- 1.5 < etta && etta < -0.1) {
																    phi = PHI(momy[itrack], momx[itrack]);
																    Qx_tpcE += Qx(phi);
																    Qy_tpcE += Qy(phi);
																    // qW += WEIGHT;
																    flag = true;
																		tpcE_phi[i][j].push_back(phi);
														    }

																// TPC(W)
														    if (0.1 < etta && etta < 1.5) {
																    phi = PHI(momy[itrack], momx[itrack]);
																    Qx_tpcW += Qx(phi);
																    Qy_tpcW += Qy(phi);
																    // qW += WEIGHT;
																    flag = true;
																		tpcW_phi[i][j].push_back(phi);
														    }

												  }

											}
											etta = 0.0;
											momModule = 0.0;
											momTransverse = 0.0;
											phi = 0.0;
											/* End of conditions */
										}
										if (flag == false) { continue;}


										Psi_EP_tpcE = PSI_EP(Qy_tpcE, Qx_tpcE);
										Psi_EP_tpcW = PSI_EP(Qy_tpcW, Qx_tpcW);

										Psi_EP_tpcEmat[i][j] = Psi_EP_tpcE;
										Psi_EP_tpcWmat[i][j] = Psi_EP_tpcW;
								    pre_Resolution = cos(HARMONICNUMBER * (Psi_EP_tpcE - Psi_EP_tpcW));
								    //pre_Resolution /= qW;
								    hRes[i]->Fill(pre_Resolution);

								    error = 0.0;
								    flag = false;
								    // qW = 0.0;
								    Qx_tpcE = 0.0;
								    Qy_tpcE = 0.0;
								    Qx_tpcW = 0.0;
								    Qy_tpcW = 0.0;
								    Psi_EP_tpcE = 0.0;
								    Psi_EP_tpcW = 0.0;
								    pre_Resolution = 0.0;
								   //if (Cut(ientry) < 0) continue;
								}

							}

					 }

				}

			// Fill histograms here:


			for (Int_t j = 0; j < npTbins; j++) {
				for (Int_t i = 0; i < nbins; i++ ) {
					Resolution = sqrt(hRes[i]->GetMean());
					tmpE = vec_Cos_N_Phi_Psi(tpcE_phi[i][j], Psi_EP_tpcWmat[i][j]);
					tmpW = vec_Cos_N_Phi_Psi(tpcW_phi[i][j], Psi_EP_tpcEmat[i][j]);

					for (Int_t k = 0; k < tmpE.size(); k++) {
						hFlowtpcE[i][j]->Fill(tmpE[k] / Resolution);
					}

					for (Int_t k = 0; k < tmpW.size(); k++) {
						hFlowtpcW[i][j]->Fill(tmpW[k] / Resolution);
					}

					hFlowtpcE[i][j]->Write();
					hFlowtpcW[i][j]->Write();

					Resolution = 0.0;
					tmpE.clear();
					tmpW.clear();
				}
			}

			// Save output to the file
			std::cout << "Save output information to the file: " << outFileName.Data() << std::endl;
			fo->cd();
			fo->Close();

			std::cout << "Program is finished successfully!" << std::endl;

			// Print out timer info
			timer.Stop();
			timer.Print();


}


