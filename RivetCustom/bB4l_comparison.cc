// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/Vector4.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  struct bB4l_comparison_Plots { };



  /// Top pair production with central jet veto
  class bB4l_comparison : public Analysis {
  public:

    /// Constructor
    bB4l_comparison()
      : Analysis("bB4l_comparison"),

    _jet_ntag(0)
    //_hMap(),
    //_histLimit(6)
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 100);
      addProjection(fs, "ALL_FS");

      /// Get electrons from truth record
      IdentifiedFinalState elec_fs(Cuts::abseta < 100);
      elec_fs.acceptIdPair(PID::ELECTRON);
      addProjection(elec_fs, "ELEC_FS");

      /// Get muons which pass the initial kinematic cuts:
      IdentifiedFinalState muon_fs(Cuts::abseta < 100);
      muon_fs.acceptIdPair(PID::MUON);
      addProjection(muon_fs, "MUON_FS");

      /// Get all neutrinos. These will not be used to form jets.
      /// We'll use the highest 2 pT neutrinos to calculate the MET
      IdentifiedFinalState neutrino_fs(Cuts::abseta < 100);
      neutrino_fs.acceptNeutrinos();
      addProjection(neutrino_fs, "NEUTRINO_FS");

      // Final state used as input for jet-finding.
      // We include everything except the muons and neutrinos
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      jet_input.addVetoPairId(PID::MUON);
      addProjection(jet_input, "JET_INPUT");

      // Get the jets
      FastJets jets(jet_input, FastJets::ANTIKT, 0.5);
      addProjection(jets, "JETS");

/*
      for (unsigned int ihist = 0; ihist < _histLimit ; ihist++) {
        const unsigned int threshLimit = _thresholdLimit(ihist);
        for (unsigned int ithres = 0; ithres < threshLimit; ithres++) {
          _histogram(ihist, ithres); // Create all histograms
        }
      }
*/

	// Create histograms
	_histPt1 = bookHisto1D("Pt1", 25, 0, 400);
	_histEta1 = bookHisto1D("Eta1", 20, -2.5, 2.5);
	_histPhiemu = bookHisto1D("Phiemu", 20, 0, 2*PI);
	_histdeltaRb = bookHisto1D("deltaRb", 20, 0, 5);
	_histdeltaRl = bookHisto1D("deltaRl", 20, 0, 5);
	_histPtMiss = bookHisto1D("PtMiss", 25, 0, 400);
	_histHT = bookHisto1D("HT", 20, 0, 1200);
//	_histMlb = bookHisto1D("Mlb", 50, 0, 350);

	_histmWjB = bookHisto1D("mWjB", 40, 150, 200);
	_histmljB = bookHisto1D("mljB", 50, 0, 350);
	_histmjB  = bookHisto1D("mjB", 20, 5, 50);
	_histdeltaR = bookHisto1D("deltaR", 100, 0, 2);
	_histxB   = bookHisto1D("xB", 20, 0, 1);
	_histpTbDec = bookHisto1D("pTbDec", 15, 0, 30);


	_histPt1_T = bookHisto1D("Pt1_T", 25, 0, 400);
        _histEta1_T = bookHisto1D("Eta1_T", 20, -2.5, 2.5);
        _histPhiemu_T = bookHisto1D("Phiemu_T", 20, 0, 2*PI);
        _histdeltaRb_T = bookHisto1D("deltaRb_T", 20, 0, 5);
        _histdeltaRl_T = bookHisto1D("deltaRl_T", 20, 0, 5);
        _histPtMiss_T = bookHisto1D("PtMiss_T", 25, 0, 400);
        _histHT_T = bookHisto1D("HT_T", 20, 0, 1200);
//      _histMlb = bookHisto1D("Mlb", 50, 0, 350);

        _histmWjB_T = bookHisto1D("mWjB_T", 40, 150, 200);
        _histmljB_T = bookHisto1D("mljB_T", 50, 0, 350);
        _histmjB_T  = bookHisto1D("mjB_T", 20, 5, 50);
        _histdeltaR_T = bookHisto1D("deltaR_T", 100, 0, 2);
        _histxB_T   = bookHisto1D("xB_T", 20, 0, 1);
        _histpTbDec_T = bookHisto1D("pTbDec_T", 15, 0, 30);
//	_histmlb = bookHisto1D("mlb", 50, 0, 350);
//	_histmLeptonsjB = bookHisto1D("mLeptonsjB", 50, 0, 350);

	_hist_selected = bookHisto1D("selected", 5, 0, 5);
	_hist_ptlepton = bookHisto1D("ptlepton", 10, 0, 400);
	_nevents = 0;
	_nbjets = 0;
	_nleptons = 0;
	
   }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      bool hadronization = true;

      _nevents += 1;

      double weight = event.weight();
//      const double normWW = 1.0/3.9465;
//      weight *= normWW;

      /// Get the various sets of final state particles
      Particles elecFS = applyProjection<IdentifiedFinalState>(event, "ELEC_FS").particlesByPt();
      Particles muonFS = applyProjection<IdentifiedFinalState>(event, "MUON_FS").particlesByPt();
      const Particles& neutrinoFS = applyProjection<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();

      // Get all jets with pT > 30 GeV
      const Jets& jets = applyProjection<FastJets>(event, "JETS").jetsByPt();

      // Keep any jets that pass the initial rapidity cut
      vector<const Jet*> central_jets;
      foreach(const Jet& j, jets) {
        if (j.absrap() < 2.5) central_jets.push_back(&j);
      }


      // Get b hadrons with pT > 5 GeV
      /// @todo This is a hack -- replace with UnstableFinalState
      GenParticle B_hadrons;
      GenParticle B_bbar_hadrons;
      vector<GenParticle const *> allParticles = particles(event.genEvent());
      double maxpTHad = 0;
      double maxpTHadBar = 0;
      bool foundB=false;
      bool foundBbar=false;
      for (size_t i = 0; i < allParticles.size(); i++) {
        const GenParticle* p = allParticles[i];
        if (!PID::isHadron(p->pdg_id()) || !PID::hasBottom(p->pdg_id())) continue;
        //if (p->momentum().perp() < 5*GeV) continue;
	if (Particle(*p).hasAncestor(5) && Particle(*p).pT()>maxpTHad){
	  B_hadrons = *p;
	  maxpTHad = Particle(*p).pT();
	  foundB=true;
	}
	if (Particle(*p).hasAncestor(-5) && Particle(*p).pT()>maxpTHadBar){
	  B_bbar_hadrons = *p;
	  maxpTHadBar = Particle(*p).pT();
	  foundBbar=true;
	}

	if(foundB) B_hadrons.print();
	if(foundBbar) B_bbar_hadrons.print();


//	cout << (p->pdg_id()>0 && Particle(*p).pT()>maxpTHad) << "\n";
//	cout << i << ": " << Particle(*p).pT() << "\n";
      }

      vector<const Jet*> b_jets;
      foreach (const Jet* j, central_jets) {
//        foreach (const GenParticle* b, B_hadrons) {
        if (j->containsParticle(B_hadrons)) b_jets.push_back(j);
	else if (!hadronization && j->containsBottom()) b_jets.push_back(j);
//        }
      }

      vector<const Jet*> b_bar_jets;
      foreach (const Jet* j, central_jets) {
	if (j->containsParticle(B_bbar_hadrons)) b_bar_jets.push_back(j);
	else if (!hadronization && j->containsBottom()) b_bar_jets.push_back(j);
      }

      vector<const Jet*> b_jets_T;
      foreach (const Jet* j, central_jets) {
//        foreach (const GenParticle* b, B_hadrons) {
        if (j->containsParticleId(5)) b_jets_T.push_back(j);
//        }
      }

      vector<const Jet*> b_bar_jets_T;
      foreach (const Jet* j, central_jets) {
        if (j->containsParticleId(-5)) b_bar_jets_T.push_back(j);
      }
 
      if (!b_jets_T.empty() && !b_bar_jets_T.empty()) _nbjets += 1;

      vector<const GenParticle*> WPlus;
      vector<const GenParticle*> WMinus;
      foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
           //cout << p->pdg_id() << " ";
          if (p->pdg_id() == 24){ WPlus.push_back(p);}
	  if (p->pdg_id() == -24){ WMinus.push_back(p);}
           //const GenVertex* pv = p->production_vertex();
      }

 
      vector<const GenParticle*> leptonsP, leptonsM, bs, bbars;
      if (!WPlus.empty() && !WMinus.empty()){
      for (HepMC::GenVertex::particle_iterator iter = WPlus[0]->end_vertex()->particles_begin(HepMC::descendants); iter != WPlus[0]->end_vertex()->particles_end(HepMC::descendants); ++iter) {
	if ((*iter)->pdg_id()==-11 || (*iter)->pdg_id()==-13) {leptonsP.push_back(*iter);}
     }
      for (HepMC::GenVertex::particle_iterator iter = WMinus[0]->end_vertex()->particles_begin(HepMC::descendants); iter != WMinus[0]->end_vertex()->particles_end(HepMC::descendants); ++iter) {
        if ((*iter)->pdg_id()==11 || (*iter)->pdg_id()==13) {leptonsM.push_back(*iter);}
     }

      for (HepMC::GenVertex::particle_iterator iter = WPlus[0]->production_vertex()->particles_begin(HepMC::descendants); iter != WPlus[0]->production_vertex()->particles_end(HepMC::descendants); ++iter) {
        if ((*iter)->pdg_id()==5) {bs.push_back(*iter);}
     }
      for (HepMC::GenVertex::particle_iterator iter = WMinus[0]->production_vertex()->particles_begin(HepMC::descendants); iter != WMinus[0]->production_vertex()->particles_end(HepMC::descendants); ++iter) {
        if ((*iter)->pdg_id()==-5) {bbars.push_back(*iter);}
     }}


// TRUTH

     vector<const GenParticle*> leptonsP_T;
     vector<const GenParticle*> leptonsM_T;
     vector<const GenParticle*> bs_T;
     vector<const GenParticle*> bbars_T;
     bool foundM = false, foundP = false, foundb = false, foundbbar = false;
      foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
	p->print();
        if (foundM == false && p->pdg_id()==13) {leptonsM_T.push_back(p); foundM=true;}
	if (foundP == false && p->pdg_id()==-11) {leptonsP_T.push_back(p); foundP=true;}
	if (foundb == false && p->pdg_id()==5) {bs_T.push_back(p); foundb=true;}
	if (foundbbar == false && p->pdg_id()==-5) {bbars_T.push_back(p); foundbbar=true;}
     }

     if (!leptonsM_T.empty() && !leptonsP_T.empty()) _nleptons += 1;

      // Get b hadrons with pT > 5 GeV
      //       /// @todo This is a hack -- replace with UnstableFinalState
      vector<GenParticle const *> B_hadrons_old;
      for (size_t i = 0; i < allParticles.size(); i++) {
         const GenParticle* p2 = allParticles[i];
         if (!PID::isHadron(p2->pdg_id()) || !PID::hasBottom(p2->pdg_id())) continue;
                                //if (p->momentum().perp() < 5*GeV) continue;
         B_hadrons_old.push_back(p2);
      }
      //
      // For each of the good jets, check whether any are b-jets (via dR matching)
      vector<const Jet*> b_jets_old;
      foreach (const Jet* j, central_jets) {
        bool isbJet = false;
        foreach (const GenParticle* b, B_hadrons_old) {
          if (deltaR(j->momentum(), FourMomentum(b->momentum())) < 0.4) isbJet = true;
        }
        if (isbJet) b_jets_old.push_back(j);
      }

      // Get the MET by taking the vector sum of all neutrinos
      /// @todo Use MissingMomentum instead?
      double MET = 0;
      FourMomentum p_MET;
      foreach (const Particle& p, neutrinoFS) {
        p_MET = p_MET + p.momentum();
      }
      MET = p_MET.pT();


      // Get the electrons and muons if there aren't any in the final state
      //(event.genEvent())->print();
      if(elecFS.empty() || muonFS.empty()) {
         foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
	   //cout << p->pdg_id() << " ";
           if (p->pdg_id() != -11 && p->pdg_id() != 13) continue;
           //const GenVertex* pv = p->production_vertex();
           bool passed = true;
           //if (pv) {
           //  foreach (const GenParticle* pp, particles_in(pv)) {
           //    if ( p->pdg_id() == pp->pdg_id() ) {
           //      passed = false;
           //      break;
           //    }
           //  }
           //}
      	if (passed && p->pdg_id() == -11 && Particle(*p).hasAncestor(24)) {
	 if (!elecFS.empty()) {
		if (Particle(*p).pT() > elecFS[0].pT()) elecFS[0] = Particle(*p);
	 }
	 else elecFS.push_back(Particle(*p));
	}
      	if (passed && p->pdg_id() == 13 && Particle(*p).hasAncestor(-24)) {
	 if (!muonFS.empty()) {
		if (Particle(*p).pT() > muonFS[0].pT()) muonFS[0] = Particle(*p);
	 }
	 else muonFS.push_back(Particle(*p));
	}
      }
      }

      if(!WPlus.empty()) WPlus[0]->print();
      if(!WMinus.empty()) WMinus[0]->print();
      if(!leptonsP.empty()) leptonsP[0]->print();
      if(!leptonsM.empty()) leptonsM[0]->print();

      cout << leptonsP_T.size() + leptonsM_T.size() << "\n";

      // Finally, the same again with the emu channel
      //if (elecFS.size() == 1 && muonFS.size() == 1) {
        // With the desired charge signs
        //if (charge(elecFS[0]) >= 0 && charge(muonFS[0]) <= 0) {
          // Calculate HT: scalar sum of the pTs of the leptons and all good jets
/*          double HT = 0;
	  if(elecFS.size() >= 1 && muonFS.size() >= 1) {
          HT += elecFS[0].pT();
          HT += muonFS[0].pT();
	  }
          foreach (const Jet* j, central_jets)
            HT += fabs(j->pT());
          // Keep events with HT > 130 GeV
          if (HT > 130.0*GeV) {
            // And again we want 2 or more b-jets
            if (b_jets.size() > 1 && elecFS.size() >= 1 && muonFS.size() >= 1) {
		 if (MET >= 20.0*GeV) {
              		passed_emu = true;

			// Additional requirements for the invariant mass event selection
			if(elecFS.size() >= 1 && muonFS.size() >= 1) {
			if(elecFS[0].pT() >= 20.*GeV && muonFS[0].pT() >= 20.*GeV && fabs(elecFS[0].eta()) < 2.5 && fabs(muonFS[0].eta()) < 2.5) {

			if(deltaR(b_jets[0]->momentum(), elecFS[0].momentum()) >= 0.4 && deltaR(b_jets[0]->momentum(), muonFS[0].momentum()) >= 0.4
			&& deltaR(b_jets[1]->momentum(), elecFS[0].momentum()) >= 0.4 && deltaR(b_jets[1]->momentum(), muonFS[0].momentum()) >= 0.4) {
				if(MET >= 60*GeV && (elecFS[0].momentum() + muonFS[0].momentum()).mass() >= 15*GeV) {
					if(abs((elecFS[0].momentum() + muonFS[0].momentum()).mass() - 91*GeV) >= 10*GeV) {
						passed_Mlb = true;
					}
				}
			}
			}
			}
		}
            }
	}
	//}
	//}

      if (passed_emu == true) {
*/
/*
	double pTcut[4] = {30.,40.,60.,80.};
	// Count the jet multiplicity for 30, 40, 60 and 80GeV
      	unsigned int ithres, jet_n[4] = {0,0,0,0};
      	for (ithres = 0; ithres < 4; ithres++) {
		foreach(const Jet* j, central_jets) {
        		if (j->pT() > pTcut[ithres]) jet_n[ithres]++;
		}
      	unsigned int ncutoff[4] = {8,7,6,5};
      	if (jet_n[ithres] > ncutoff[ithres]) jet_n[ithres] = ncutoff[ithres];
      }


      // Fill histograms
      for (unsigned int ihist = 0; ihist < 6; ihist++) {
        unsigned int threshLimit = _thresholdLimit(ihist);
        for (ithres = 0; ithres < threshLimit; ithres++) {
          if (jet_n[ithres] < 2) continue; // 2 or more jets for ljets
          // Fill
          if (ihist == 0) _histogram(ihist, ithres)->fill(jet_n[ithres], weight); // njets
          else if (ihist == 1) _histogram(ihist, ithres)->fill(central_jets[0]->pT(), weight); // leading jet pT
          else if (ihist == 2) _histogram(ihist, ithres)->fill(central_jets[1]->pT(), weight); // 2nd jet pT
          else if (ihist == 3 && jet_n[ithres] >= 3) _histogram(ihist, ithres)->fill(central_jets[2]->pT(), weight); // 3rd jet pT
          else if (ihist == 4 && jet_n[ithres] >= 4) _histogram(ihist, ithres)->fill(central_jets[3]->pT(), weight); // 4th jet pT
          else if (ihist == 5 && jet_n[ithres] >= 5) _histogram(ihist, ithres)->fill(central_jets[4]->pT(), weight); // 5th jet pT
        }
      }
*/
    if(b_jets.size()>=1 && b_bar_jets.size()>=1) {
	_histPt1->fill(0.5*(b_jets[0]->pT()+b_bar_jets[0]->pT()), weight);
	_histEta1->fill(0.5*(b_jets[0]->eta()+b_bar_jets[0]->eta()), weight);
	_histdeltaRb->fill(deltaR(b_jets[0]->momentum(), b_bar_jets[0]->momentum()), weight);
		//_histPhiemu->fill(deltaPhi(elecFS[0], muonFS[0]), weight);
		//_histPhiemu->fill(mapAngle0To2Pi(elecFS[0].momentum().phi() - muonFS[0].momentum().phi()), weight);
		//_histdeltaRl->fill(deltaR(elecFS[0], muonFS[0]), weight);
	if(!WPlus.empty() && !WMinus.empty()) _histmWjB->fill(((Particle(WPlus[0]).momentum()+b_jets[0]->momentum()).mass() + (Particle(WMinus[0]).momentum()+b_bar_jets[0]->momentum()).mass())/2, weight);
	_histmjB->fill((b_jets[0]->momentum().mass() + b_bar_jets[0]->momentum().mass())/2, weight);
	if(!leptonsP.empty() && !leptonsM.empty()) {
		 _histmljB->fill(((Particle(leptonsP[0]).momentum()+b_jets[0]->momentum()).mass() + (Particle(leptonsM[0]).momentum()+b_bar_jets[0]->momentum()).mass())/2, weight);
}

	
     }
	_histPtMiss->fill(MET, weight);

    if(b_jets_T.size()>=1 && b_bar_jets_T.size()>=1) {
        _histPt1_T->fill(0.5*(b_jets_T[0]->pT()+b_bar_jets_T[0]->pT()), weight);
        _histEta1_T->fill(0.5*(b_jets_T[0]->eta()+b_bar_jets_T[0]->eta()), weight);
        _histdeltaRb_T->fill(deltaR(b_jets_T[0]->momentum(), b_bar_jets_T[0]->momentum()), weight);
        if (!WPlus.empty() && !WMinus.empty()) _histmWjB_T->fill(((Particle(WPlus[0]).momentum()+b_jets[0]->momentum()).mass() + (Particle(WMinus[0]).momentum()+b_bar_jets[0]->momentum()).mass())/2, weight);
        _histmjB_T->fill((b_jets_T[0]->momentum().mass() + b_bar_jets_T[0]->momentum().mass())/2, weight);
        if(!leptonsP_T.empty() && !leptonsM_T.empty()) {
                 _histmljB_T->fill(((Particle(leptonsP_T[0]).momentum()+b_jets_T[0]->momentum()).mass() + (Particle(leptonsM_T[0]).momentum()+b_bar_jets_T[0]->momentum()).mass())/2, weight);

}


     }
        _histPtMiss_T->fill(MET, weight);

	if(!leptonsP_T.empty() && !leptonsM_T.empty()){
        //_histPhiemu_T->fill(deltaPhi(leptonsP_T[0], leptonsM_T[0]), weight);
        _histPhiemu_T->fill(mapAngle0To2Pi(Particle(leptonsP_T[0]).momentum().phi() - Particle(leptonsM_T[0]).momentum().phi()), weight);
        _histdeltaRl->fill(deltaR(Particle(leptonsP_T[0]), Particle(leptonsM_T[0])), weight);
	
	_hist_ptlepton->fill(1/2.0*(Particle(leptonsP_T[0]).pT()+Particle(leptonsM_T[0]).pT()), weight);
}

	//_histHT->fill(HT, weight);
	//if(passed_Mlb == true) {
//        if(b_jets_old.size() >= 2 && !elecFS.empty() && !muonFS.empty()) {
//		if((elecFS[0].momentum() + b_jets_old[0]->momentum()).mass() + (muonFS[0].momentum() + b_jets_old[1]->momentum()).mass()
//		> (elecFS[0].momentum() + b_jets_old[1]->momentum()).mass() + (muonFS[0].momentum() + b_jets_old[0]->momentum()).mass()) {
//		_histMlb->fill(((elecFS[0].momentum() + b_jets_old[1]->momentum()).mass() + (muonFS[0].momentum() + b_jets_old[0]->momentum()).mass())/2, weight);
//		}
//		else _histMlb->fill(((elecFS[0].momentum() + b_jets_old[0]->momentum()).mass() + (muonFS[0].momentum() + b_jets_old[1]->momentum()).mass())/2, weight);
//	}
//	if(b_jets.size()>=1 && b_bar_jets.size()>=1 && !elecFS.empty() && !muonFS.empty()){
//	  _histmLeptonsjB->fill(((elecFS[0].momentum()+b_jets[0]->momentum()).mass() + (muonFS[0].momentum()+b_bar_jets[0]->momentum()).mass())/2, weight);
//	}
//	if(bs.size()>=1 && bbars.size()>=1 && !elecFS.empty() && !muonFS.empty()){
//	   _histmlb->fill(((elecFS[0].momentum()+Particle(bs[0]).momentum()).mass() + (muonFS[0].momentum()+Particle(bbars[0]).momentum()).mass())/2, weight);
//	}
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {
//	normalize(_histPt1);
//	normalize(_histEta1);
//	normalize(_histPhiemu);
//	normalize(_histdeltaRb);
//	normalize(_histdeltaRl);
//	normalize(_histPtMiss);
//	normalize(_histHT);
//	normalize(_histMlb);

	const double norm = crossSection()/sumOfWeights();
	scale(_histPt1, norm);
	scale(_histEta1, norm);
	scale(_histPhiemu, norm);
	scale(_histdeltaRb, norm);
	scale(_histdeltaRl, norm);
	scale(_histPtMiss, norm);
	scale(_histHT, norm);
//	scale(_histMlb, norm);
	scale(_histmWjB, norm);
	scale(_histmljB, norm);
	scale(_histmjB, norm);

	scale(_histPt1_T, norm);
        scale(_histEta1_T, norm);
        scale(_histPhiemu_T, norm);
        scale(_histdeltaRb_T, norm);
        scale(_histdeltaRl_T, norm);
        scale(_histPtMiss_T, norm);
        scale(_histHT_T, norm);
//      scale(_histMlb, norm);
        scale(_histmWjB_T, norm);
        scale(_histmljB_T, norm);
        scale(_histmjB_T, norm);

	scale(_hist_ptlepton, norm);
//	scale(_histmLeptonsjB, norm);
//        scale(_histmlb, norm);
      //const double norm = crossSection()/sumOfWeights();
      //typedef map<unsigned int, Histo1DPtr>::value_type IDtoHisto1DPtr; ///< @todo Remove when C++11 allowed
      //foreach (IDtoHisto1DPtr ihpair, _hMap) scale(ihpair.second, norm); ///< @todo Use normalize(ihpair.second, crossSection())

	_hist_selected->fill(0.5, _nevents);
	_hist_selected->fill(1.5, _nbjets);
	_hist_selected->fill(2.5, _nleptons);
    }



  private:

    /// @name Histogram helper functions
    //@{

	Histo1DPtr _histPt1;
	Histo1DPtr _histEta1;
	Histo1DPtr _histPhiemu;
	Histo1DPtr _histdeltaRb;
	Histo1DPtr _histdeltaRl;
	Histo1DPtr _histPtMiss;
	Histo1DPtr _histHT;
//	Histo1DPtr _histMlb;

	Histo1DPtr _histPt1_T;
        Histo1DPtr _histEta1_T;
        Histo1DPtr _histPhiemu_T;
        Histo1DPtr _histdeltaRb_T;
        Histo1DPtr _histdeltaRl_T;
        Histo1DPtr _histPtMiss_T;
        Histo1DPtr _histHT_T;
//      Histo1DPtr _histMlb;

//	Histo1DPtr _histmlb;
//        Histo1DPtr _histmLeptonsjB;

	Histo1DPtr _histmWjB;
        Histo1DPtr _histmljB;
        Histo1DPtr _histmjB;
        Histo1DPtr _histdeltaR;
        Histo1DPtr _histxB;
        Histo1DPtr _histpTbDec;

	Histo1DPtr _histmWjB_T;
        Histo1DPtr _histmljB_T;
        Histo1DPtr _histmjB_T;
        Histo1DPtr _histdeltaR_T;
        Histo1DPtr _histxB_T;
        Histo1DPtr _histpTbDec_T;

	Histo1DPtr _hist_selected;
	Histo1DPtr _hist_ptlepton;
	int _nevents;
	int _nbjets;
	int _nleptons;

/*
    unsigned int _thresholdLimit(unsigned int histId) {
      if (histId == 0) return 4;
      return 1;
    }

    Histo1DPtr _histogram(unsigned int histId, unsigned int thresholdId) {
      assert(histId < _histLimit);
      assert(thresholdId < _thresholdLimit(histId));

      const unsigned int hInd = (histId == 0) ? thresholdId : (_thresholdLimit(0) + (histId-1) + thresholdId);
      if (_hMap.find(hInd) != _hMap.end()) return _hMap[hInd];

      if (histId == 0) _hMap.insert(make_pair(hInd,bookHisto1D(1,thresholdId+1,1)));
      else _hMap.insert(make_pair(hInd,bookHisto1D(2,histId,1)));
      return _hMap[hInd];
    }
*/

private:

    unsigned int _jet_ntag;

    //map<unsigned int, Histo1DPtr> _hMap;
    //unsigned int _histLimit;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(bB4l_comparison);

}
