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

namespace Rivet {


  struct ATLAS_LUDOVIC3_Plots { };



  /// Top pair production with central jet veto
  class ATLAS_LUDOVIC3 : public Analysis {
  public:

    /// Constructor
    ATLAS_LUDOVIC3()
      : Analysis("ATLAS_LUDOVIC3"),

    _jet_ntag(0),
    _hMap(),
    _histLimit(6)
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 5);
      addProjection(fs, "ALL_FS");

      /// Get electrons from truth record
      IdentifiedFinalState elec_fs(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      elec_fs.acceptIdPair(PID::ELECTRON);
      addProjection(elec_fs, "ELEC_FS");

      /// Get muons which pass the initial kinematic cuts:
      IdentifiedFinalState muon_fs(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      muon_fs.acceptIdPair(PID::MUON);
      addProjection(muon_fs, "MUON_FS");

      /// Get all neutrinos. These will not be used to form jets.
      /// We'll use the highest 2 pT neutrinos to calculate the MET
      IdentifiedFinalState neutrino_fs(Cuts::abseta < 5.0);
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

      for (unsigned int ihist = 0; ihist < _histLimit ; ihist++) {
        const unsigned int threshLimit = _thresholdLimit(ihist);
        for (unsigned int ithres = 0; ithres < threshLimit; ithres++) {
          _histogram(ihist, ithres); // Create all histograms
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      /// Get the various sets of final state particles
      const Particles& elecFS = applyProjection<IdentifiedFinalState>(event, "ELEC_FS").particlesByPt();
      const Particles& muonFS = applyProjection<IdentifiedFinalState>(event, "MUON_FS").particlesByPt();
      const Particles& neutrinoFS = applyProjection<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();

      // Get all jets with pT > 30 GeV
      const Jets& jets = applyProjection<FastJets>(event, "JETS").jetsByPt(30.0*GeV);

      // Keep any jets that pass the initial rapidity cut
      vector<const Jet*> central_jets;
      foreach(const Jet& j, jets) {
        if (j.absrap() < 2.5) central_jets.push_back(&j);
      }


      // Get b hadrons with pT > 5 GeV
      /// @todo This is a hack -- replace with UnstableFinalState
      vector<GenParticle const *> B_hadrons;
      vector<GenParticle const *> allParticles = particles(event.genEvent());
      for (size_t i = 0; i < allParticles.size(); i++) {
        const GenParticle* p = allParticles[i];
        if (!PID::isHadron(p->pdg_id()) || !PID::hasBottom(p->pdg_id())) continue;
        //if (p->momentum().perp() < 5*GeV) continue;
        B_hadrons.push_back(p);
      }

      // For each of the good jets, check whether any are b-jets (via dR matching)
      vector<const Jet*> b_jets;
      foreach (const Jet* j, central_jets) {
        bool isbJet = false;
        foreach (const GenParticle* b, B_hadrons) {
          if (deltaR(j->momentum(), FourMomentum(b->momentum())) < 0.35) isbJet = true;
        }
        if (isbJet) b_jets.push_back(j);
      }

      // Get the MET by taking the vector sum of all neutrinos
      /// @todo Use MissingMomentum instead?
      double MET = 0;
      FourMomentum p_MET;
      foreach (const Particle& p, neutrinoFS) {
        p_MET = p_MET + p.momentum();
      }
      MET = p_MET.pT();

      bool passed_emu = false;
      // Finally, the same again with the emu channel
      //if (elecFS.size() >= 1 && muonFS.size() >= 1) {
        // With the desired charge signs
        //if (charge(elecFS[0]) >= 0 && charge(muonFS[0]) <= 0) {
          // Calculate HT: scalar sum of the pTs of the leptons and all good jets
          double HT = 0;
          if(elecFS.size() >= 1 && muonFS.size() >= 1) {
	  HT += elecFS[0].pT();
          HT += muonFS[0].pT();
	  }
          foreach (const Jet* j, central_jets)
            HT += fabs(j->pT());
          // Keep events with HT > 130 GeV
          if (HT > 130.0*GeV) {
            // And again we want 2 or more b-jets
            if (b_jets.size() > 1) {
		 if (MET >= 20.0*GeV) {
              		passed_emu = true;
		}
            }
	}
	//}
	//}

      if (passed_emu == true) {

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

      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = crossSection()/sumOfWeights();
      typedef map<unsigned int, Histo1DPtr>::value_type IDtoHisto1DPtr; ///< @todo Remove when C++11 allowed
      foreach (IDtoHisto1DPtr ihpair, _hMap) scale(ihpair.second, norm); ///< @todo Use normalize(ihpair.second, crossSection())
    }



  private:

    /// @name Histogram helper functions
    //@{

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

private:

    unsigned int _jet_ntag;

    map<unsigned int, Histo1DPtr> _hMap;
    unsigned int _histLimit;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_LUDOVIC3);

}
