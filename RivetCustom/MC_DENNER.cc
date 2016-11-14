// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Tools/RivetMT2.hh"
#include "Rivet/Tools/mt2_bisect.hh"

namespace Rivet {


  class MC_DENNER : public Analysis {

#include "NLOHisto1D.cc"

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_DENNER()
      : Analysis("MC_DENNER")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    void init() {

       VetoedFinalState fs;
       fs.addVetoPairId(PID::ELECTRON);
       fs.addVetoPairId(PID::MUON);
       fs.addVetoPairId(PID::NU_E);
       fs.addVetoPairId(PID::NU_MU);
       addProjection(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");

       IdentifiedFinalState elecs;
       elecs.acceptIdPair(PID::ELECTRON);
       addProjection(elecs, "elecs");

       IdentifiedFinalState muons;
       muons.acceptIdPair(PID::MUON);
       addProjection(muons, "muons");

       addProjection(VisibleFinalState(-5, 5),"vfs");


       _h_xs = bookNLOHisto1D("xs", 1, 0.5, 1.5);
       _h_xscut = bookNLOHisto1D("xscut", 1, 0.5, 1.5);
       _h_pTbmin = bookNLOHisto1D("pTbmin", 25, 0., 400.);
       _h_pTbmax = bookNLOHisto1D("pTbmax", 25, 0., 400.);
       _h_pTmu = bookNLOHisto1D("pTmu", 40, 0., 1000.);
       _h_etamu = bookNLOHisto1D("etamu", 50, -5, 5);
       _h_mll = bookNLOHisto1D("mll", 40, 0., 400.);
       _h_mlb = bookNLOHisto1D("mlb", 25, 0., 200.);
       _h_etdphi = bookNLOHisto1D("etdphi", 40, 0., 400.);
       _h_mt2 = bookNLOHisto1D("mt2", 50, 0., 500.);
       _h_malt = bookNLOHisto1D("malt", 40, 0., 400.);
       _h_pTlb = bookNLOHisto1D("pTlb", 40, 0., 800.);
       _h_pTl = bookNLOHisto1D("pTl", 40, 0., 800.);
       _h_pTb = bookNLOHisto1D("pTb", 40, 0., 800.);
    }


    void analyze(const Event& event) {
_h_xs->fill(1 , event);

//beginning of cuts ---------------------------------------
	// TODO: Change for hadron level running. More than one electron.
       ParticleVector elec = applyProjection<IdentifiedFinalState>(event, "elecs").particlesByPt();
	if (elec.size() !=1){
		vetoEvent;
	}
	if (elec[0].momentum().pT()<20 || fabs(elec[0].momentum().eta())>2.5){
		vetoEvent;
	}
       ParticleVector muon = applyProjection<IdentifiedFinalState>(event, "muons").particlesByPt();
	if (muon.size() !=1){
		vetoEvent;
	}
	if (muon[0].momentum().pT()<20 || fabs(muon[0].momentum().eta())>2.5){
		vetoEvent;
	}

	Jets bjets;
	const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
	const Jets alljets = jetpro.jetsByPt();
//cuts for ttbarjet-----------------------------------------
//if (jetpro.jetsByPt(15.*GeV).size()<3){
//	vetoEvent;
//}
//cuts for ttbarjet-----------------------------------------

	double jetpt = 0;
	foreach (const Jet& jet, alljets) {
		jetpt += jet.momentum().pT(); //for H_T observable
		if ( fabs( jet.momentum().eta() ) < 2.5 && jet.momentum().pT() > 30  ) {
			if (jet.containsBottom() ) {
				bjets.push_back(jet);
			}
		}
	}

	if (bjets.size() < 2) {
		vetoEvent;
	}

	if(elec[0].momentum().pT() + muon[0].momentum().pT() + bjets[0].momentum().pT() + bjets[1].momentum().pT()  < 130*GeV)
	{
		vetoEvent;
	}

	ParticleVector vfs_particles = applyProjection<VisibleFinalState>(event, "vfs").particles();
	FourMomentum pTmiss;
	foreach ( const Particle & p, vfs_particles ) {
		pTmiss -= p.momentum();
	}

	if(pTmiss.pT()  < 20*GeV)
	{
		vetoEvent;
	}

//end of cuts ---------------------------------------


	if( bjets[0].momentum().pT() > bjets[1].momentum().pT() ){
		_h_pTbmax->fill(bjets[0].momentum().pT(), event);
		_h_pTbmin->fill(bjets[1].momentum().pT(), event);
	}else{
		_h_pTbmax->fill(bjets[1].momentum().pT(), event);
		_h_pTbmin->fill(bjets[0].momentum().pT(), event);
	}
	_h_pTmu->fill(muon[0].momentum().pT(), event);
	_h_etamu->fill(muon[0].momentum().eta(), event);
	_h_mll->fill((muon[0].momentum()+elec[0].momentum()).mass(), event);

	_h_pTl->fill(elec[0].momentum().pT(), event);
	_h_pTl->fill(muon[0].momentum().pT(), event);
	_h_pTb->fill(bjets[0].momentum().pT(), event);
	_h_pTb->fill(bjets[1].momentum().pT(), event);

	FourMomentum p00 = elec[0].momentum() + bjets[0].momentum();
	FourMomentum p10 = muon[0].momentum() + bjets[0].momentum();
	FourMomentum p01 = elec[0].momentum() + bjets[1].momentum();
	FourMomentum p11 = muon[0].momentum() + bjets[1].momentum();
	double mlb;
	double etdphi;
	double m_T2;
	if((p00.mass()+p11.mass())<(p01.mass()+p10.mass())){
		mlb =  (p00.mass() + p11.mass())/2; 
		double pa[3]    = {p00.mass() , p00.x() , p00.y() };
		double pb[3]    = {p11.mass() , p11.x() , p11.y() };
		double pmiss[3] = {-999.999 , pTmiss.x() , pTmiss.y() };
		mt2_bisect::mt2 mt2_event;
		mt2_event.set_momenta( pa, pb, pmiss );
		mt2_event.set_mn( 0.0 );
		m_T2 = mt2_event.get_mt2();
		etdphi=(elec[0].momentum().Et()*mapAngle0ToPi(elec[0].momentum().phi()-bjets[0].momentum().phi()) + muon[0].momentum().Et()*mapAngle0ToPi(muon[0].momentum().phi()-bjets[1].momentum().phi()))/2;
		_h_pTlb->fill(p00.pT(), event);
		_h_pTlb->fill(p11.pT(), event);
	} else {
		mlb =  (p10.mass() + p01.mass())/2;  
		double pa[3]    = {p10.mass() , p10.x() , p10.y() };
		double pb[3]    = {p01.mass() , p01.x() , p01.y() };
		double pmiss[3] = {-999.999 , pTmiss.x() , pTmiss.y() };
		mt2_bisect::mt2 mt2_event;
		mt2_event.set_momenta( pa, pb, pmiss );
		mt2_event.set_mn( 0.0 );
		m_T2 = mt2_event.get_mt2();
		etdphi=(elec[0].momentum().Et()*mapAngle0ToPi(elec[0].momentum().phi()-bjets[1].momentum().phi()) + muon[0].momentum().Et()*mapAngle0ToPi(muon[0].momentum().phi()-bjets[0].momentum().phi()))/2;
		_h_pTlb->fill(p01.pT(), event);
		_h_pTlb->fill(p10.pT(), event);
	}
	_h_mlb->fill(mlb , event);
	_h_mt2->fill(m_T2 , event);
	_h_etdphi->fill(etdphi , event);

	double malt;
	if(elec[0].momentum().vector3().unit().dot(bjets[0].momentum().vector3().unit())+muon[0].momentum().vector3().unit().dot(bjets[1].momentum().vector3().unit()) > elec[0].momentum().vector3().unit().dot(bjets[1].momentum().vector3().unit())+muon[0].momentum().vector3().unit().dot(bjets[0].momentum().vector3().unit())){
		malt=sqrt((elec[0].momentum().E()-elec[0].momentum().vector3().dot(bjets[0].momentum().vector3().unit()))*(muon[0].momentum().E()-muon[0].momentum().vector3().dot(bjets[1].momentum().vector3().unit())));
	} else {
		malt=sqrt((elec[0].momentum().E()-elec[0].momentum().vector3().dot(bjets[1].momentum().vector3().unit()))*(muon[0].momentum().E()-muon[0].momentum().vector3().dot(bjets[0].momentum().vector3().unit())));
	}

_h_malt->fill(malt,event);
_h_xscut->fill(1 , event);

    }

    void finalize() {

double factor = crossSection()/sumOfWeights();

 _h_xs->finalize();
scale(_h_xs, factor);
 _h_xscut->finalize();
scale(_h_xscut, factor);
 _h_pTbmin->finalize();
scale(_h_pTbmin, factor);
 _h_pTbmax->finalize();
scale(_h_pTbmax, factor);
 _h_pTmu->finalize();
scale(_h_pTmu, factor);
 _h_etamu->finalize();
scale(_h_etamu, factor);
 _h_mll->finalize();
scale(_h_mll, factor);
 _h_mlb->finalize();
scale(_h_mlb, factor);
 _h_etdphi->finalize();
scale(_h_etdphi, factor);
 _h_mt2->finalize();
scale(_h_mt2, factor);
 _h_malt->finalize();
scale(_h_malt, factor);
 _h_pTlb->finalize();
scale(_h_pTlb, factor);
 _h_pTl->finalize();
scale(_h_pTl, factor);
 _h_pTb->finalize();
scale(_h_pTb, factor);
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    //@{
  NLOHisto1DPtr  _h_xs;
  NLOHisto1DPtr  _h_xscut;
  NLOHisto1DPtr  _h_pTbmin;
  NLOHisto1DPtr  _h_pTbmax;
  NLOHisto1DPtr  _h_pTmu;
  NLOHisto1DPtr  _h_etamu;
  NLOHisto1DPtr  _h_mll;
  NLOHisto1DPtr  _h_mlb;
  NLOHisto1DPtr  _h_etdphi;
  NLOHisto1DPtr  _h_mt2;
  NLOHisto1DPtr  _h_malt;
  NLOHisto1DPtr  _h_pTlb;
  NLOHisto1DPtr  _h_pTl;
  NLOHisto1DPtr  _h_pTb;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_DENNER);

}
