double linear_eq(double h,double v1, double h1, double v2, double h2)
  {
    // think of v: vertical axis, h: horizontal axis
    double m = (v2-v1)/(h2-h1);       // slope
    double b = (v1*h2-h1*v2)/(h2-h1); // intercept
    return m*h+b;  // value of vertical at horizontal = h
 }

double get_vdiff(std::pair<TVector3,TVector3> crt, std::pair<TVector3,TVector3> tpc)
{
  double z = 200.0; // calculate y diff at z = 200cm
  double c_crt = linear_eq(z, crt.first.Y(), crt.first.Z(),crt.second.Y(),crt.second.Z());
  double c_tpc = linear_eq(z, tpc.first.Y(), tpc.first.Z(),tpc.second.Y(),tpc.second.Z());
  return std::abs(c_crt-c_tpc);
}

std::pair<TVector3,TVector3> x_correct(std::pair<TVector3,TVector3> crt_trk,
				       std::pair<TVector3,TVector3> tpc_trk,
				       double thetaxz,
				       double z_tpctrk_center
				       )
{
  // creates TPC track to match CRT track at z = 200
  // call only after final matching
  // returns std::pair of new endpoints (as TVector3)

  
  // first get the x position of the crt track at z=200, we will move the
  // matched tpc track in x by this amount
  double x_crt_trk = linear_eq(z_tpctrk_center, crt_trk.first.X(), crt_trk.first.Z(),crt_trk.second.X(),crt_trk.second.Z());
  // check that at that point crt track is in the correct tpc
  double thetaxz_rad = thetaxz*TMath::DegToRad();
  TVector3 corrected_tpc_endpoint1 {
    x_crt_trk-((z_tpctrk_center-tpc_trk.first.Z())*tan(thetaxz_rad)),
    tpc_trk.first.Y(),
    tpc_trk.first.Z()
  };
  TVector3 corrected_tpc_endpoint2 {
    x_crt_trk+((tpc_trk.second.Z()-z_tpctrk_center)*tan(thetaxz_rad)),
    tpc_trk.second.Y(),
    tpc_trk.second.Z()
  };
  return std::make_pair(corrected_tpc_endpoint1, corrected_tpc_endpoint2);
}


void remake_ct(TString modifier = "")
{
  gROOT->SetStyle("uboone_sty");
  // modify the angles of the muontracks and make angle range -90 to
  // 90, with z2>z1
  bool modify_angles = true;
  bool remake_ct = false; // remake crt tracks from crt hits

  // output file for corrected muon tracks
  if ( modifier != "" ) {
    modifier = "_"+modifier;
  }
  TFile *of = new TFile(Form("trees/muon_hitdumper_NS%s_corrected.root",modifier.Data()),"recreate");
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
  std::vector<double> v_mhit_charge;
  std::vector<int> v_mhit_wire;
  int out_run;
  int out_subrun;
  int out_evt;
  // is the crt track restricted to just tpc 0 or tpc 1? or does it
  // cross from -x to x or x to -x?
  std::vector<int> v_tpc_tpc;
  outtree->Branch("run",&out_run);
  outtree->Branch("subrun",&out_subrun);
  outtree->Branch("evt",&out_evt);
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
  outtree->Branch("muontrk_tpc",&v_tpc_tpc);
  outtree->Branch("mhit_charge",&v_mhit_charge);
  outtree->Branch("mhit_wire",&v_mhit_wire);
  
  // input file
  TString filename {Form("trees/muon_hitdumper_NS%s.root",modifier.Data())};
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

  TTreeReaderArray<double> ct_x1(reader,"ct_x1");
  TTreeReaderArray<double> ct_y1(reader,"ct_y1");
  TTreeReaderArray<double> ct_z1(reader,"ct_z1");
  TTreeReaderArray<double> ct_x2(reader,"ct_x2");
  TTreeReaderArray<double> ct_y2(reader,"ct_y2");
  TTreeReaderArray<double> ct_z2(reader,"ct_z2");
  TTreeReaderArray<double> ct_t(reader,"ct_time");

  
  // muon tracks
  TTreeReaderArray<float> muontrk_x1(reader,"muontrk_x1");
  TTreeReaderArray<float> muontrk_y1(reader,"muontrk_y1");
  TTreeReaderArray<float> muontrk_z1(reader,"muontrk_z1");
  TTreeReaderArray<float> muontrk_x2(reader,"muontrk_x2");
  TTreeReaderArray<float> muontrk_y2(reader,"muontrk_y2");
  TTreeReaderArray<float> muontrk_z2(reader,"muontrk_z2");
  TTreeReaderArray<double> muontrk_t0(reader,"muontrk_t0");
  TTreeReaderArray<int> muontrk_tpc(reader,"muontrk_tpc");
  TTreeReaderArray<float> muontrk_th_xz(reader,"muontrk_theta_xz");
  TTreeReaderArray<float> muontrk_th_yz(reader,"muontrk_theta_yz");

  // muon hits
  /**
     these are the hits from which the muontracks were originally
     calculated, it's what we use in the python EVD to make the
     (wire,tick) event displays. they are very long vectors (average
     of 3k hits per muon track) and have an entry called mhit_trk that
     should match the index of the muon track we are matching. We also
     need to match to save only from the collection plane using the
     variable mhit_plane. 

     We will definitely save the mhit_charge for such hit. Other
     variables we might save are mhit_wire and mhit_peakT.

     mhit_peakT should be worthless as we do time matching anyway, but
     i'm not 100% sure.

     mhit_wire can be scaled easily to match the z coordinate.  We
     could use to get the other coordinates (x and y). Bear in mind
     that track hits are sometimes saved out of order (e.g. a track
     that crosses wires 400 to 1300 could store first the hits from
     1000 to 1300 and then go backwards from 1000 to 400). The
     mhit_wire variable becomes absolutely necessary if we want to
     make our landau+gauss histograms from tracks that cross different
     x-bins. If we just use angles which cross a single x-bin then we
     can ignore the position information altogether.
  */
  TTreeReaderArray<int> mhit_trk(reader,"mhit_trk");
  TTreeReaderArray<int> mhit_plane(reader,"mhit_plane");
  TTreeReaderArray<int> mhit_wire(reader,"mhit_wire");
  TTreeReaderArray<double> mhit_charge(reader,"mhit_charge");
  

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
  TH2F *h_match2 = new TH2F("h_match2",
			   ";CRT-muon |#Delta#theta_{XZ}| (deg);CRT-muon |#Delta#theta_{YZ}| (deg)",
			   100,0,10,100,0,10);

  TH1F * h_types = new TH1F("h_types",";Track type for remade CRT tracks; Entries",6,0,6);
  TH1F * h_nct = new TH1F("h_nct",";Number of remade CRT tracks; Entries",20,0,20);

  // angle validation
  TH1F *h_tpc_thetaxz = new TH1F("h_tpc_thetaxz",";TPC #theta_{XZ};Entries",60,-180,180);
  TH1F *h_tpc_thetayz = new TH1F("h_tpc_thetayz",";TPC #theta_{YZ};Entries",60,-180,180);
  TH1F *h_crt_thetaxz = new TH1F("h_crt_thetaxz",";CRT #theta_{XZ};Entries",60,-180,180);
  TH1F *h_crt_thetayz = new TH1F("h_crt_thetayz",";CRT #theta_{YZ};Entries",60,-180,180);
  TH1F *h_n_matched = new TH1F("h_n_matched",";Number of matches per event;Entries",5,-0.5,4.5);
  gStyle->SetOptStat(1);
  int events = -1;
  int n_total_muontrks = 0;
  int n_total_crttrks  = 0;
  int n_total_matched  = 0;
  int n_eve_matched  = 0;

  outtree->Branch("event",&events);
  outtree->Branch("n_matched",&n_eve_matched);


  std::cout << "Reading file" << std::endl;
  while ( reader.Next() ) {
    out_run = *run;
    out_subrun = *subrun;
    out_evt = *evt;
    n_eve_matched = 0;
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
    v_tpc_tpc.clear();
    v_mhit_charge.clear();
    v_mhit_wire.clear();
    
    events++;
    //if ( events > 409 ) break;

    // muontrks loop
    std::vector<std::pair<TVector3,TVector3>> v_muontrks;
    std::vector<std::pair<double,double>> v_muontrks_theta_xz_yz;
    std::vector<int> v_muontrks_tpc;
    for ( int i = 0; i < muontrk_x1.GetSize(); i++ ) {
      auto v_mtrk1 =  TVector3{muontrk_x1[i],muontrk_y1[i],muontrk_z1[i]};
      auto v_mtrk2 =  TVector3{muontrk_x2[i],muontrk_y2[i],muontrk_z2[i]};

      double modified_thetaxz = muontrk_th_xz[i];
      double modified_thetayz = muontrk_th_yz[i];
      if ( modify_angles ) {
	// make z2 > z1 always
	if ( v_mtrk2.Z() < v_mtrk1.Z() ) {
	  TVector3 v_temp;
	  v_temp = v_mtrk2;
	  v_mtrk2 = v_mtrk1;
	  v_mtrk1 = v_temp;
	}
	// modify angles
	if ( muontrk_th_xz[i] > 90 && muontrk_tpc[i] == 0 ) {
	  modified_thetaxz = muontrk_th_xz[i] - 180.;
	} else if ( muontrk_th_xz[i] > 90 && muontrk_tpc[i] == 1 ) {
	  modified_thetaxz = 180. - muontrk_th_xz[i];
	} else if ( muontrk_th_xz[i] < 90 && muontrk_tpc[i] == 1 ) {
	  modified_thetaxz = -muontrk_th_xz[i];
	}
	if ( muontrk_th_yz[i] < -90 ) {
	  modified_thetayz = muontrk_th_yz[i] + 180.;
	}
	if ( muontrk_th_yz[i] > 90 ) {
	  modified_thetayz = muontrk_th_yz[i] - 180.;
	}
      }
      v_muontrks.push_back(std::make_pair(v_mtrk1,v_mtrk2));
      v_muontrks_theta_xz_yz.push_back(std::make_pair(modified_thetaxz,modified_thetayz));
      v_muontrks_tpc.push_back(muontrk_tpc[i]);
      n_total_muontrks++;
    }
    
    // crt hits loop
    // we will be looking for up-downstream crossing tracks

    std::vector<std::pair<TVector3,TVector3>> v_crttrks;
    std::vector<double> v_crttrks_t;
    std::vector<int> v_crttrks_type;

    int n_ct = 0;
    if ( !remake_ct ) {
      // use crt tracks stored in the file
      for ( int i = 0 ; i < ct_x1.GetSize(); i++ ) {
	TVector3 ct1 {ct_x1[i],ct_y1[i],ct_z1[i]};
	TVector3 ct2 {ct_x2[i],ct_y2[i],ct_z2[i]};
	// make second one always downstream
	if ( ct_z2[i] < ct_z1[i] ) {
	  TVector3 cttemp;
	  cttemp = ct2;
	  ct2 = ct1;
	  ct1 = cttemp;
	}
	// upstream hit
	if ( ct1.Z() < -179.06 || ct1.Z() > -179.04 ) {
	  continue;
	}
	// downstream hit
	if ( ct2.Z() < 770.94 || ct2.Z() > 770.96 ) {
	  continue;
	}
	v_crttrks.push_back(std::make_pair(ct1,ct2));
	v_crttrks_t.push_back(ct_t[i]);
	v_crttrks_type.push_back(4);
	n_ct++;
	n_total_crttrks++;
      }
    } else { // (remake_ct == true)
      // remake crt tracks
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
	    v_crttrks_type.push_back(4); // because hitplanes 2 and 1
	    if ( count_crttrks == true ) {
	      n_total_crttrks++;
	      count_crttrks = false;
	    }
	  }
	} // end of j looping
      } // end of i crt hits loop
    }
    h_nct->Fill(n_ct);

    // muontrk / crttrk matching
    for ( int m = 0; m < v_muontrks_theta_xz_yz.size(); m++ ) {
      bool b_event_match2 = false;

      double best_match   = 1000.0;
      double best_match_i = -1;
      double best_match2   = 1000.0;
      double best_match_i2 = -1;
      double best_crt_thetaxz = -999;
      double best_crt_thetayz = -999;
      for ( int c = 0; c < v_crttrks.size(); c++ ) {
	// focus on matching "tpcs" for now
	
	auto crttrk1 = v_crttrks[c].first;
	auto crttrk2 = v_crttrks[c].second;
	auto crttrk = crttrk2 - crttrk1;

	// check that at the middle z point of the muon track, the crt
	// track is in the same TPC as the muon track
	double z_middle_tpctrk = (v_muontrks[m].first.Z()+v_muontrks[m].second.Z())/2.0;
	double x_of_crt_at_tpc_middle = linear_eq(z_middle_tpctrk,
						  crttrk1.X(),
						  crttrk1.Z(),
						  crttrk2.X(),
						  crttrk2.Z()
						  );
	int crt_tpc = -1;
	if ( x_of_crt_at_tpc_middle > 0 ) {
	  crt_tpc = 1;
	} else if ( x_of_crt_at_tpc_middle < 0 ) {
	  crt_tpc = 0;
	}
	if ( crt_tpc != v_muontrks_tpc[m] ) {
	  continue;
	}

	/** previous method, which doesn't match what's on the muontrackproducer
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
	*/
	double dx = crttrk2.X() - crttrk1.X();
	double dy = crttrk2.Y() - crttrk1.Y();
	double dz = crttrk2.Z() - crttrk1.Z();
	double th_xz = atan2(dx,dz)*TMath::RadToDeg();
	double th_yz = atan2(dy,dz)*TMath::RadToDeg();
	//if ( th_xz < 0.0 ) th_xz = -th_xz;
	h_crt_thetaxz->Fill(th_xz);
	h_crt_thetayz->Fill(th_yz);
	if ( !modify_angles && ( crttrk2.Y() - crttrk1.Y() ) < 0.0 ) th_yz = -th_yz;

	double deltaxz = std::abs(v_muontrks_theta_xz_yz[m].first-th_xz);
	double deltayz = std::abs(v_muontrks_theta_xz_yz[m].second-th_yz);
	if ( ( deltaxz < 2 ) && (deltayz < 2 ) ) {
	  // match using theta xz and theta yz
	  if ( sqrt(deltaxz*deltaxz + deltayz*deltayz) < best_match2 ) {
	    best_match_i2 = c;
	    best_match2 = sqrt(deltaxz*deltaxz + deltayz*deltayz);
	    best_crt_thetaxz = th_xz;
	    best_crt_thetayz = th_yz;
	  }
	  if ( b_event_match2 == false ) {
	    b_event_match2 = true;
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
	int c = best_match_i2;
	double vdiff = get_vdiff(v_crttrks[c],v_muontrks[m]);
	if ( vdiff < 2.0 ) { // cm
	  n_total_matched++;
	  n_eve_matched++;
	  
	  // do not comment this line
	  h_match2->Fill(std::abs(v_muontrks_theta_xz_yz[m].first  - best_crt_thetaxz),
			 std::abs(v_muontrks_theta_xz_yz[m].second - best_crt_thetayz));
	  double z_tpctrk_center = (v_muontrks[m].first.Z()+v_muontrks[m].second.Z())/2.0;
	  auto corrected_tpc = x_correct(v_crttrks[c],v_muontrks[m],v_muontrks_theta_xz_yz[m].first,z_tpctrk_center);
	  h_tpc_thetaxz->Fill(v_muontrks_theta_xz_yz[m].first);
	  h_tpc_thetayz->Fill(v_muontrks_theta_xz_yz[m].second);
	  //h_tpc_thetayz->Fill(muontrk_th_yz[i]);

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
	  v_tpc_tpc.push_back(v_muontrks_tpc[m]);
	  // sadly need to loop over the mhits now
	  for ( int mh = 0; mh < mhit_trk.GetSize(); mh++ ) {
	    if ( mhit_trk[mh] != m ) { // match current hit
	      continue;
	    }
	    if ( mhit_plane[mh] != 2 ) { // only save collection plane
	      continue;
	    }
	    v_mhit_charge.push_back(mhit_charge[mh]);
	    v_mhit_wire.push_back(mhit_wire[mh]);
	    
	  }
	}
      }
      
    } // end loop over muontrks
    if ( v_crt_x1.size() > 0 ) {
      h_n_matched->Fill(n_eve_matched);
      outtree->Fill();
    }

  }
  of->Write();


  std::cout << Form("Looped over %i events", events) << std::endl;
  std::cout << Form("Total:\nMuontrk: %i\nCRTtrk : %i\nMatched: %i\n",
		    n_total_muontrks,
		    n_total_crttrks,
		    n_total_matched);
  // remove leading '_' from filename modifier
  if ( modifier[0] == '_' ) {
    modifier.Remove(0,1); // remove 1 character, starting from char 0
  }
  if ( modifier == "" ) {
    modifier = "default";
  }
  // mkdir ( name, recursive)
  gSystem->mkdir(Form("img/ct_remake/%s",modifier.Data()),true);
  
  // 1d histos
  auto *c1 = new TCanvas();
  h_nct->Draw();
  c1->SaveAs(Form("img/ct_remake/%s/ud_nct.png",modifier.Data()));
  c1->SaveAs(Form("img/ct_remake/%s/ud_nct.pdf",modifier.Data()));
  delete c1;

  TStyle *ts = gROOT->GetStyle("uboone_sty_colz");
  if ( ts != nullptr ) ts->cd();
  
  auto *c2 = new TCanvas();
  // 2d histos
  gStyle->SetOptStat(1);
  h_crt_type4->Draw("colz");
  c2->SaveAs(Form("img/ct_remake/%s/h_crt_type4.png",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/h_crt_type4.pdf",modifier.Data()));
  h_crt_typeall->Draw("colz");
  c2->SaveAs(Form("img/ct_remake/%s/h_crt_typeall.png",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/h_crt_typeall.pdf",modifier.Data()));
  h_match2->Draw("colz");
  c2->SaveAs(Form("img/ct_remake/%s/ud_match2.pdf",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/ud_match2.png",modifier.Data()));

  // angle validation histos
  h_tpc_thetaxz->Draw();
  c2->SaveAs(Form("img/ct_remake/%s/tpc_thetaxz.pdf",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/tpc_thetaxz.png",modifier.Data()));
  h_tpc_thetayz->Draw();
  c2->SaveAs(Form("img/ct_remake/%s/tpc_thetayz.pdf",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/tpc_thetayz.png",modifier.Data()));
  
  h_crt_thetaxz->Draw();
  c2->SaveAs(Form("img/ct_remake/%s/crt_thetaxz.pdf",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/crt_thetaxz.png",modifier.Data()));
  h_crt_thetayz->Draw();
  c2->SaveAs(Form("img/ct_remake/%s/crt_thetayz.pdf",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/crt_thetayz.png",modifier.Data()));
  // number of matches histo
  c2->SetLogy();
  h_n_matched->Draw();
  c2->SaveAs(Form("img/ct_remake/%s/n_matches.pdf",modifier.Data()));
  c2->SaveAs(Form("img/ct_remake/%s/n_matches.png",modifier.Data()));
  
}
