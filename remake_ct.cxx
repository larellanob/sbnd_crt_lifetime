double linear_eq(double z,double y1, double z1, double y2, double z2)
  {
    double m = (y2-y1)/(z2-z1);
    double b = (y1*z2-z1*y2)/(z2-z1);
    return m*z+b;
 }

double get_vdiff(std::pair<TVector3,TVector3> crt, std::pair<TVector3,TVector3> tpc)
{
  double z = 200.0; // calculate y diff at z = 200cm
  double c_crt = linear_eq(z, crt.first.Y(), crt.first.Z(),crt.second.Y(),crt.second.Z());
  double c_tpc = linear_eq(z, tpc.first.Y(), tpc.first.Z(),tpc.second.Y(),tpc.second.Z());
  return std::abs(c_crt-c_tpc);
}

std::pair<TVector3,TVector3> x_correct(std::pair<TVector3,TVector3> crt,
				       std::pair<TVector3,TVector3> tpc,
				       double thetaxz
				       )
{
  // crate TPC track to match CRT track at z = 200
  // call only after final matching
  double z = 200.0;
  double x_crt = linear_eq(z, crt.first.X(), crt.first.Z(),crt.second.X(),crt.second.Z());
  double thetaxz_rad = thetaxz*TMath::DegToRad();
  TVector3 corrected_tpc1 {
    x_crt-((z-tpc.first.Z())*tan(thetaxz_rad)),
    tpc.first.Y(),
    tpc.first.Z()
  };
  TVector3 corrected_tpc2 {
    x_crt-((tpc.second.Z()-z)*tan(thetaxz_rad)),
    tpc.second.Y(),
    tpc.second.Z()
  }; 
  return std::make_pair(corrected_tpc1, corrected_tpc2);
}


void remake_ct()
{
  gROOT->SetStyle("uboone_sty");

  // output file for corrected muon tracks
  TFile *of = new TFile("CORRECTED_muon_hitdumper_NS.root","recreate");
  of->mkdir("hitdumper");
  of->cd("hitdumper");
  TTree *outtree = new TTree("hitdumpertree_corr","Corrected TPC x-coordinate and matched to CRT");
  std::vector<double> v_crt_x1;
  std::vector<double> v_crt_y1;
  std::vector<double> v_crt_z1;
  std::vector<double> v_crt_x2;
  std::vector<double> v_crt_y2;
  std::vector<double> v_crt_z2;
  std::vector<double> v_tpc_x1;
  std::vector<double> v_tpc_y1;
  std::vector<double> v_tpc_z1;
  std::vector<double> v_tpc_x2;
  std::vector<double> v_tpc_y2;
  std::vector<double> v_tpc_z2;
  std::vector<double> v_tpc_uncorrected_x1;
  std::vector<double> v_tpc_uncorrected_y1;
  std::vector<double> v_tpc_uncorrected_z1;
  std::vector<double> v_tpc_uncorrected_x2;
  std::vector<double> v_tpc_uncorrected_y2;
  std::vector<double> v_tpc_uncorrected_z2;
  std::vector<double> crt_thetaxz;
  std::vector<double> crt_thetayz;
  std::vector<double> tpc_thetaxz;
  std::vector<double> tpc_thetayz;
  outtree->Branch("ct_x1",&v_crt_x1);
  outtree->Branch("ct_y1",&v_crt_y1);
  outtree->Branch("ct_z1",&v_crt_z1);
  outtree->Branch("ct_x2",&v_crt_x2);
  outtree->Branch("ct_y2",&v_crt_y2);
  outtree->Branch("ct_z2",&v_crt_z2);
  outtree->Branch("muontrk_x1",&v_tpc_x1);
  outtree->Branch("muontrk_y1",&v_tpc_y1);
  outtree->Branch("muontrk_z1",&v_tpc_z1);
  outtree->Branch("muontrk_x2",&v_tpc_x2);
  outtree->Branch("muontrk_y2",&v_tpc_y2);
  outtree->Branch("muontrk_z2",&v_tpc_z2);
  outtree->Branch("muontrk_uncorrected_x1",&v_tpc_uncorrected_x1);
  outtree->Branch("muontrk_uncorrected_y1",&v_tpc_uncorrected_y1);
  outtree->Branch("muontrk_uncorrected_z1",&v_tpc_uncorrected_z1);
  outtree->Branch("muontrk_uncorrected_x2",&v_tpc_uncorrected_x2);
  outtree->Branch("muontrk_uncorrected_y2",&v_tpc_uncorrected_y2);
  outtree->Branch("muontrk_uncorrected_z2",&v_tpc_uncorrected_z2);
  outtree->Branch("crt_thetaxz",&crt_thetaxz);
  outtree->Branch("crt_thetayz",&crt_thetayz);
  outtree->Branch("tpc_thetaxz",&tpc_thetaxz);
  outtree->Branch("tpc_thetayz",&tpc_thetayz);
  
  // input file
  TString filename {"muon_hitdumper_NS.root"};
  TFile *f = new TFile(filename);
  TTreeReader reader("hitdumper/hitdumpertree",f);

  // evt run subrun
  TTreeReaderValue<int> run (reader,"run");
  TTreeReaderValue<int> subrun (reader,"subrun");
  TTreeReaderValue<int> evt (reader,"event");

  // crt hits
  TTreeReaderArray<double> chit_x(reader,"chit_x");
  TTreeReaderArray<double> chit_y(reader,"chit_y");
  TTreeReaderArray<double> chit_z(reader,"chit_z");
  TTreeReaderArray<double> chit_t(reader,"chit_time");
  TTreeReaderArray<int> chit_plane(reader,"chit_plane");

  // crt tracks
  /*
  TTreeReaderArray<double> ct_z1(reader,"ct_z1");
  TTreeReaderArray<double> ct_z2(reader,"ct_z2");
  TTreeReaderArray<double> ct_t(reader,"ct_time");
  */
  
  // muon tracks
  TTreeReaderArray<float> muontrk_x1(reader,"muontrk_x1");
  TTreeReaderArray<float> muontrk_y1(reader,"muontrk_y1");
  TTreeReaderArray<float> muontrk_z1(reader,"muontrk_z1");
  TTreeReaderArray<float> muontrk_x2(reader,"muontrk_x2");
  TTreeReaderArray<float> muontrk_y2(reader,"muontrk_y2");
  TTreeReaderArray<float> muontrk_z2(reader,"muontrk_z2");
  TTreeReaderArray<double> muontrk_t0(reader,"muontrk_t0");
  TTreeReaderArray<float> muontrk_th_xz(reader,"muontrk_theta_xz");
  TTreeReaderArray<float> muontrk_th_yz(reader,"muontrk_theta_yz");
  

  TH1F *h_ud_theta = new TH1F("h_ud_theta",";CRT Track #theta (deg); Entries",100,0,60);
  TH1F *h_ud_phi = new TH1F("h_ud_phi",";CRT Track #phi (deg);Entries",100,-180,180);
  TH2F *h_crt_typeall = new TH2F("h_crt_typeall",
				 ";CRT track distance travelled (cm); CRT track travel time (ns)",
				 100,0,2000,100,0,100);
  TH2F *h_crt_type0 = new TH2F("h_crt_type0",
			    ";CRT track distance travelled (U-D) (cm); CRT track travel time (ns)",
			    100,900,1500,50,0.0,100);
  TH2F *h_crt_type1 = new TH2F("h_crt_type1",
			    ";CRT track distance travelled (U-D) (cm); CRT track travel time (ns)",
			    100,900,1500,50,0.0,100);
  TH2F *h_crt_type2 = new TH2F("h_crt_type2",
			    ";CRT track distance travelled (U-D) (cm); CRT track travel time (ns)",
			    100,900,1500,50,0.0,100);
  TH2F *h_crt_type3 = new TH2F("h_crt_type3",
			    ";CRT track distance travelled (U-D) (cm); CRT track travel time (ns)",
			    100,900,1500,50,0.0,100);
  TH2F *h_crt_type4 = new TH2F("h_crt_type4",
			    ";CRT track distance travelled (U-D) (cm); CRT track travel time (ns)",
			    100,900,1500,50,0.0,100);
  TH2F *h_crt_type5 = new TH2F("h_crt_type5",
			    ";CRT track distance travelled (U-D) (cm); CRT track travel time (ns)",
			    100,900,1500,50,0.0,100);
  TH2F *h_match = new TH2F("h_match",
			   ";CRT-muon |#Delta#theta| (deg);CRT-muon |#Delta#phi| (deg)",
			   100,0,10,100,0,10);
  TH2F *h_match2 = new TH2F("h_match2",
			   ";CRT-muon |#Delta#theta_{XZ}| (deg);CRT-muon |#Delta#theta_{YZ}| (deg)",
			   100,0,10,100,0,10);

  TH1F * h_types = new TH1F("h_types",";Track type for remade CRT tracks; Entries",6,0,6);
  TH1F * h_nct = new TH1F("h_nct",";Number of remade CRT tracks; Entries",20,0,20);
  gStyle->SetOptStat(1);
  int events = -1;
  int n_total_muontrks = 0;
  int n_total_crttrks  = 0;
  int n_total_matched  = 0;
  int n_total_matched2  = 0;
  int n_eve_matched2  = 0;

  outtree->Branch("event",&events);
  outtree->Branch("n_matched",&n_eve_matched2);


  while ( reader.Next() ) {
    n_eve_matched2 = 0;
    v_crt_x1.clear();
    v_crt_y1.clear();
    v_crt_z1.clear();
    v_crt_x2.clear();
    v_crt_y2.clear();
    v_crt_z2.clear();
    v_tpc_x1.clear();
    v_tpc_y1.clear();
    v_tpc_z1.clear();
    v_tpc_x2.clear();
    v_tpc_y2.clear();
    v_tpc_z2.clear();
    v_tpc_uncorrected_x1.clear();
    v_tpc_uncorrected_y1.clear();
    v_tpc_uncorrected_z1.clear();
    v_tpc_uncorrected_x2.clear();
    v_tpc_uncorrected_y2.clear();
    v_tpc_uncorrected_z2.clear();
    crt_thetaxz.clear();
    crt_thetayz.clear();
    tpc_thetaxz.clear();
    tpc_thetayz.clear();
    
    events++;
    //if ( events < 308 ) continue;
    //if ( events > 312 ) break;
    /*
    if ( events < 9100 ) {
      events++;
      continue;
    }
    */

    // muontrks loop
    std::vector<std::pair<TVector3,TVector3>> v_muontrks;
    std::vector<double> v_muontrks_t;
    std::vector<std::pair<double,double>> v_muontrks_theta_xz_yz;
    for ( int i = 0; i < muontrk_x1.GetSize(); i++ ) {
      TVector3 v_mtrk = (
			 TVector3{muontrk_x1[i],muontrk_y1[i],muontrk_z1[i]}
			 - TVector3{muontrk_x2[i],muontrk_y2[i],muontrk_z2[i]}
			 );
      auto v_mtrk1 =  TVector3{muontrk_x1[i],muontrk_y1[i],muontrk_z1[i]};
      auto v_mtrk2 =  TVector3{muontrk_x2[i],muontrk_y2[i],muontrk_z2[i]};
      v_muontrks.push_back(std::make_pair(v_mtrk1,v_mtrk2));
      v_muontrks_t.push_back(muontrk_t0[i]);
      v_muontrks_theta_xz_yz.push_back(std::make_pair(muontrk_th_xz[i],muontrk_th_yz[i]));
      n_total_muontrks++;
    }
    
    // crt hits loop
    // we will be looking for up-downstream crossing tracks

    std::vector<std::pair<TVector3,TVector3>> v_crttrks;
    std::vector<double> v_crttrks_t;
    std::vector<int> v_crttrks_type;

    int n_ct = 0;
    for ( int i = 0; i < chit_z.GetSize(); i++ ) {
      // will match this "i" chit to a "j" one
      bool count_crttrks = true; // will not count same i for different js
      /**
	 plane 2 is downstream. hits seem to be stored in same order
	 in commissioning ntuples: 0-4-2-1-6-5-3. Loop over hits until
	 you find upstream plane
      */
      for ( int j = i+1; j < chit_z.GetSize(); j++ ) {
	if ( chit_plane[j] == chit_plane[i] ) continue;
	double delta_t = std::abs(chit_t[j] - chit_t[i]);
	double avg_t = (chit_t[j]+chit_t[i])/2.0;

	// ctime is in us. we demand <100 ns coincidence
	if ( delta_t > 0.1 || delta_t < 0.00) { 
	  continue;
	}
	auto v1 = TVector3{chit_x[i],chit_y[i],chit_z[i]};
	auto v2 = TVector3{chit_x[j],chit_y[j],chit_z[j]};
	TVector3 v_ud = (
			 TVector3{chit_x[i],chit_y[i],chit_z[i]}
			 - TVector3{chit_x[j],chit_y[j],chit_z[j]}
			 );
	h_crt_typeall->Fill(v_ud.Mag(), delta_t*1000.);
	// TYPE 0: Anode-cathode crosser
	// ...
	// TYPE 4: UP DOWNSTREAM (only this one implemented)
	int crt_muon_type = -1;
	if ( ( chit_plane[i] == 2 && chit_plane[j] == 1 )
	     || ( chit_plane[i] == 1 && chit_plane[j] == 2 )
	     ) {
	  n_ct++;
	  h_crt_type4->Fill(v_ud.Mag(), delta_t*1000.);
	  // make it so always second point is further downstream
	  if ( v2.Z() < v1.Z() ) {
	    TVector3 vtemp;
	    vtemp = v2;
	    v2 = v1;
	    v1 = vtemp;
	  }
	  v_crttrks.push_back(std::make_pair(v1,v2));
	  v_crttrks_t.push_back(avg_t);
	  v_crttrks_type.push_back(4);
	  h_ud_theta->Fill(v_ud.Theta()*TMath::RadToDeg());
	  h_ud_phi->Fill(v_ud.Phi()*TMath::RadToDeg());
	  if ( count_crttrks == true ) {
	    n_total_crttrks++;
	    count_crttrks = false;
	  }
	}
      } // end of j looping
    } // end of i crt hits loop
    h_nct->Fill(n_ct);

    // muontrk / crttrk matching
    for ( int m = 0; m < v_muontrks.size(); m++ ) {
      bool b_event_match2 = false;
      auto muontrk1 = v_muontrks[m].first;
      auto muontrk2 = v_muontrks[m].second;

      double best_match   = 1000.0;
      double best_match_i = -1;
      double best_match2   = 1000.0;
      double best_match_i2 = -1;
      double best_crt_thetaxz = -999;
      double best_crt_thetayz = -999;
      for ( int c = 0; c < v_crttrks.size(); c++ ) {
	auto crttrk1 = v_crttrks[c].first;
	auto crttrk2 = v_crttrks[c].second;
	auto crttrk = crttrk2 - crttrk1;

	/** 
	    instead of using ROOTs theta and phi, calculate and use
	    theta_xz and theta_yz

	    very important: make sure you construct crttrks with
	    v2.Z() > v1.Z()
	*/
	/*
	crttrk1.Print();
	crttrk2.Print();
	std::cout << "---\n";
	*/
	TVector3 crttrk_xz = {
	  crttrk2.X() - crttrk1.X(),
	  0.0,
	  crttrk2.Z() - crttrk1.Z()
	};

	TVector3 crttrk_yz = {
	  0.0,
	  crttrk2.Y() - crttrk1.Y(),
	  crttrk2.Z() - crttrk1.Z()
	};
	double th_xz = TMath::ACos(crttrk_xz.Z()/crttrk_xz.Mag())*TMath::RadToDeg();
	double th_yz = TMath::ACos(crttrk_yz.Z()/crttrk_yz.Mag())*TMath::RadToDeg();
	if ( ( crttrk2.Y() - crttrk1.Y() ) < 0.0 ) th_yz = -th_yz;

	double deltaxz = std::abs(v_muontrks_theta_xz_yz[m].first-th_xz);
	double deltayz = std::abs(v_muontrks_theta_xz_yz[m].second-th_yz);
	if ( ( deltaxz < 2 ) && (deltayz < 2 ) ) {
	  // match using theta xz and theta yz
	  if ( sqrt(deltaxz*deltaxz + deltayz*deltayz) < best_match2 ) {
	    //if ( b_event_match2 == true ) std::cout << "better match " << std::endl;
	    best_match_i2 = c;
	    best_match2 = sqrt(deltaxz*deltaxz + deltayz*deltayz);
	    best_crt_thetaxz = th_xz;
	    best_crt_thetayz = th_yz;
	    /*
	    std::cout << "MATCHHH CRT, original TPC, corrected TPC:\n";
	    v_crttrks[c].first.Print();
	    v_crttrks[c].second.Print();
	    std::cout << "---" << std::endl;
	    v_muontrks[m].first.Print();
	    v_muontrks[m].second.Print();
	    std::cout << "---" << std::endl;
	    std::cout << v_muontrks_theta_xz_yz[m].first << " " << v_muontrks_theta_xz_yz[m].second << std::endl;
	    std::cout << th_xz << " " << th_yz << std::endl;
	    std::cout << "ctttrkxz, yz" << std::endl;
	    crttrk_xz.Print();
	    crttrk_yz.Print();
	    */
	  }
	  if ( b_event_match2 == false ) {
	    b_event_match2 = true;
	    //std::cout << "first match " << std::endl;
	  }

	}
      } // end loop over crttracks

	/**
	   at this point we should have matches based on angle, create
	   new tvector3 objects and store those and all other muoninfo
	   as "corrected_muontrk" or something
	   
	   the indices of muon track and best-match crt track are m
	   and c.
	 
	   calculate center in x direction for both crt and tpc
	   
	   i actually can't do it directly for tpc, because endpoints
	   are degenerated. need to first of all get true endpoints by
	   using the orientation and the "supposedly" non-degenerated
	   y and z endpoints
	*/
      if ( b_event_match2 ) {
	//std::cout << "Event, nmatched: " << events << " " << n_total_matched2 << std::endl;
	int c = best_match_i2;
	double vdiff = get_vdiff(v_crttrks[c],v_muontrks[m]);
	//std::cout << vdiff << std::endl;
	if ( vdiff < 2.0 ) { // cm
	  n_total_matched2++;
	  n_eve_matched2++;
	  
	  /*
	  std::cout << "CRT track: " << std::endl;
	  v_crttrks[c].first.Print();
	  v_crttrks[c].second.Print();
	  std::cout << "TPC track: " << std::endl;
	  v_muontrks[m].first.Print();
	  v_muontrks[m].second.Print();
	  std::cout << "-----------" << std::endl;
	  std::cout << "CRT xz yz: " << std::endl;
	  std::cout << best_crt_thetaxz << " " << best_crt_thetayz << std::endl;
	  std::cout << "TPC xz yz: " << std::endl;
	  std::cout << v_muontrks_theta_xz_yz[m].first << " " <<
	    v_muontrks_theta_xz_yz[m].second << std::endl;
	    
	  std::cout << "vdiff at z=200: " << std::endl;
	  std::cout << vdiff <<  tpc_ycenter << std::endl;
	  */
	  
	  // do not comment this line
	  h_match2->Fill(std::abs(v_muontrks_theta_xz_yz[m].first  - best_crt_thetaxz),
			 std::abs(v_muontrks_theta_xz_yz[m].second - best_crt_thetayz));
	  auto corrected_tpc = x_correct(v_crttrks[c],v_muontrks[m],v_muontrks_theta_xz_yz[m].first);	  
	  if ( events < 10 || 1 ) {
	    /*
	    std::cout << "CRT, original TPC, corrected TPC:\n";
	    v_crttrks[c].first.Print();
	    v_crttrks[c].second.Print();
	    std::cout << "---" << std::endl;
	    v_muontrks[m].first.Print();
	    v_muontrks[m].second.Print();
	    std::cout << "---" << std::endl;

	    corrected_tpc.first.Print();
	    corrected_tpc.second.Print();
	    std::cout << "ANGLES thetaxz thetayz CRT, TPC: " << std::endl;
	    std::cout << best_crt_thetaxz  << " " << best_crt_thetayz << std::endl;
	    std::cout << v_muontrks_theta_xz_yz[m].first << " " << v_muontrks_theta_xz_yz[m].second << std::endl;
	    std::cout << "-----------------------" << std::endl;
	    */
	  }
	  v_crt_x1.push_back(v_crttrks[c].first.X());
	  v_crt_y1.push_back(v_crttrks[c].first.Y());
	  v_crt_z1.push_back(v_crttrks[c].first.Z());
	  v_crt_x2.push_back(v_crttrks[c].second.X());
	  v_crt_y2.push_back(v_crttrks[c].second.Y());
	  v_crt_z2.push_back(v_crttrks[c].second.Z());
	  v_tpc_x1.push_back(corrected_tpc.first.X());
	  v_tpc_y1.push_back(corrected_tpc.first.Y());
	  v_tpc_z1.push_back(corrected_tpc.first.Z());
	  v_tpc_x2.push_back(corrected_tpc.second.X());
	  v_tpc_y2.push_back(corrected_tpc.second.Y());
	  v_tpc_z2.push_back(corrected_tpc.second.Z());
	  v_tpc_uncorrected_x1.push_back(v_muontrks[m].first.X());
	  v_tpc_uncorrected_y1.push_back(v_muontrks[m].first.Y());
	  v_tpc_uncorrected_z1.push_back(v_muontrks[m].first.Z());
	  v_tpc_uncorrected_x2.push_back(v_muontrks[m].second.X());
	  v_tpc_uncorrected_y2.push_back(v_muontrks[m].second.Y());
	  v_tpc_uncorrected_z2.push_back(v_muontrks[m].second.Z());
	  crt_thetaxz.push_back(best_crt_thetaxz);
	  crt_thetayz.push_back(best_crt_thetayz);
	  tpc_thetaxz.push_back(v_muontrks_theta_xz_yz[m].first);
	  tpc_thetayz.push_back(v_muontrks_theta_xz_yz[m].second);
	}
	/*
	auto crttrk1 = v_crttrks[c].first;
	auto crttrk2 = v_crttrks[c].second;
	auto crttrk = crttrk2 - crttrk1;
	TVector3 crttrk2_xz = {crttrk2.X(),crttrk1.Y(),crttrk2.Z()};
	TVector3 crttrk_xz = crttrk1 - crttrk2_xz;
	double th_yz = TMath::ACos((crttrk*crttrk_xz)/(crttrk.Mag()*crttrk_xz.Mag()))*TMath::RadToDeg();
	TVector3 unitz = {0,0,1};
	double th_xz = TMath::ACos((crttrk_xz*unitz)/(crttrk_xz.Mag()*unitz.Mag()))*TMath::RadToDeg();
	
	std::cout << v_muontrks_theta_xz_yz[m].first << " " <<  th_xz << std::endl;
	std::cout << v_muontrks_theta_xz_yz[m].second << " " <<  th_yz << std::endl;
	std::cout << "-----" << std::endl;
	*/
	

	/*
	if (    std::abs(v_muontrks_theta_xz_yz[m].first  - th_xz) < 2.0
	     && std::abs(v_muontrks_theta_xz_yz[m].second - th_yz) < 1.0
	     ) {
	  
	}
	*/  
	     
	
      }
      
    } // end loop over muontrks
    if ( v_crt_x1.size() > 0 ) {
      outtree->Fill();
    }

  }
  of->Write();


  std::cout << Form("Looped over %i events", events) << std::endl;
  std::cout << Form("Total:\nMuontrk: %i\nCRTtrk : %i\nMatched: %i\nMatched2: %i\n",
		    n_total_muontrks,
		    n_total_crttrks,
		    n_total_matched,
		    n_total_matched2);

  // 1d histos
  auto *c1 = new TCanvas();
  h_ud_theta->Draw();
  c1->SaveAs("img/ct_remake/ud_theta.png");
  c1->SaveAs("img/ct_remake/ud_theta.pdf");
  h_ud_phi->Draw();
  c1->SaveAs("img/ct_remake/ud_phi.png");
  c1->SaveAs("img/ct_remake/ud_phi.pdf");
  h_nct->Draw();
  c1->SaveAs("img/ct_remake/ud_nct.png");
  c1->SaveAs("img/ct_remake/ud_nct.pdf");
  delete c1;
  
  gROOT->GetStyle("uboone_sty_colz");
  uboone_sty_colz->cd();
  //uboone_sty_colz->SetPadRightMargin(0.2);
  auto *c2 = new TCanvas();
  // 2d histos
  gStyle->SetOptStat(1);
  h_crt_type4->Draw("colz");
  c2->SaveAs("img/ct_remake/h_crt_type4.png");
  c2->SaveAs("img/ct_remake/h_crt_type4.pdf");
  h_crt_typeall->Draw("colz");
  c2->SaveAs("img/ct_remake/h_crt_typeall.png");
  c2->SaveAs("img/ct_remake/h_crt_typeall.pdf");
  h_match->Draw("colz");
  c2->SaveAs("img/ct_remake/ud_match.pdf");
  c2->SaveAs("img/ct_remake/ud_match.png");
  h_match2->Draw("colz");
  c2->SaveAs("img/ct_remake/ud_match2.pdf");
  c2->SaveAs("img/ct_remake/ud_match2.png");
  
}
