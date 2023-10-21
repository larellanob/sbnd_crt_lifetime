void draw_muontype_text(TLatex *tl_muontype, int type_to_check) {
  if ( type_to_check == 0 ) {
    tl_muontype->DrawLatexNDC(0.2,0.83,"#splitline{Type 0}{Anode-cathode crosser}");
  } else if ( type_to_check == 1 ) {
    tl_muontype->DrawLatexNDC(0.2,0.83,"#splitline{Type 1}{Anode piercer}");
  } else if ( type_to_check == 2 ) {
    tl_muontype->DrawLatexNDC(0.2,0.83,"#splitline{Type 2}{Cathode piercer}");
  } else if ( type_to_check == 3 ) {
    tl_muontype->DrawLatexNDC(0.2,0.83,"#splitline{Type 3}{Top-bottom crosser}");
  } else if ( type_to_check == 4 ) {
    tl_muontype->DrawLatexNDC(0.2,0.83,"#splitline{Type 4}{Up-downstream crosser}");
  } else if ( type_to_check == 5 ) {
    tl_muontype->DrawLatexNDC(0.2,0.83,"#splitline{Type 5}{Other}");
  }

}

void validate_thetaxz_thetayz(int type_to_check = 1, float match_threshold = 2.0 /*degrees*/ )
{
  // if you want to create the muon tracks sorted in z with z2>z1,
  // then you have to make some changes to the angle
  bool sort_in_z = false;
  // give me totally uncorrected angle, just raw atan2
  // this will not match angles in ntuple, set to false if you want to match
  bool uncorrected = false;
  // modify the angles to match your own vision instead of this
  // tpc-dependent angle
  bool modify_angles = false;
  if ( modify_angles ) sort_in_z = true;
  /**
     the purpose of this macro is to see if I understand the thetaxz
     and thetayz angles which are stored in the commissioning ntuples,
     my first guess is that it goes from x1y1z1 to x2y2z2, and so
     similarly oriented tracks can either be close to zero or close to
     180 degrees.

     we will work with non-degenerate tracks, that is type-0, type-1
     or type-2
   */

  /**
     results so far (09-oct-2023):

     - conclusion: my implementation of theta_xz and theta_yz is
       correct, i.e. it matches the theta_xz and theta_yz branches on
       the ntuple, which was *presumably* stored before any x-axis
       corrections to the muon tracks. there's only a few <2 deg
       discrepancies for type-0

       - CORRECTION TO BELOW: i've made many more matches, now
         everything matches within 2 degrees, and only 19 don't match
         for type-0 within 1 deg. I needed to set the orientation of
         the tracks to always increase in z, and that fixed all
         unfixed matches.

     - i've made my own calculation of the theta_xz and theta_yz from
       the endpoints we're given for tpc-tracks. originally wanted to
       just look at types 0, 1 and 2, found that there is perfect
       match at all levels for type-1 (anode piercer) but not perfect
       for any other, which implies that even if the track doesn't
       make it to the anode when being x-corrected, its endpoints are
       still affected. presumably type-1 doesn't need any correction
       in the x-axis, this makes sense as it's the only type which
       1. one endpoint touches the anode
       2. the other endpoint touches a different CRT-plane

       for example I thought type-0 (anode-cathode crosser) wouldn't
       need a correction, because it touches the anode, but the other
       end of the track pierces the cathode, which doesn't have a
       CRT-plane in it
   */

  gROOT->SetStyle("uboone_sty");
  TString filename {"muon_hitdumper_NS.root"};
  TFile *f = new TFile(filename);
  TTreeReader reader("hitdumper/hitdumpertree",f);

  // evt run subrun
  TTreeReaderValue<int> run (reader,"run");
  TTreeReaderValue<int> subrun (reader,"subrun");
  TTreeReaderValue<int> evt (reader,"event");

  
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
  TTreeReaderArray<int> muontrk_type(reader,"muontrk_type");
  TTreeReaderArray<int> muontrk_tpc(reader,"muontrk_tpc");

  TH1F *h_ntuple_thetaxz = new TH1F("h_ntuple_thetaxz",";#theta_{XZ};Entries",60,-180,180);
  TH1F *h_ntuple_thetayz = new TH1F("h_ntuple_thetayz",";#theta_{YZ};Entries",60,-180,180);
  TH1F *h_own_thetaxz = new TH1F("h_own_thetaxz",";#theta_{XZ};Entries",60,-180,180);
  TH1F *h_own_thetayz = new TH1F("h_own_thetayz",";#theta_{YZ};Entries",60,-180,180);
  TH1F *h_modified_thetaxz = new TH1F("h_modified_thetaxz",";#theta_{XZ};Entries",60,-180,180);
  TH1F *h_modified_thetayz = new TH1F("h_modified_thetayz",";#theta_{YZ};Entries",60,-180,180);
  
  int good_match = 0;
  int bad_match = 0;
  int events = -1;
  while ( reader.Next() ) {
    events++;
    //if ( events > 100 ) break;
    for ( int i = 0; i < muontrk_x1.GetSize(); i++ ) {
      //if ( muontrk_type[i] != 0 and muontrk_type[i] != 1 and muontrk_type[i] != 2 ) continue;
      if (  muontrk_type[i] != type_to_check ) continue;
      h_ntuple_thetaxz->Fill(muontrk_th_xz[i]);
      h_ntuple_thetayz->Fill(muontrk_th_yz[i]);
      // tvector3 pointing from x1y1z1 to x2y2z2
      TVector3 v1 {muontrk_x1[i],muontrk_y1[i],muontrk_z1[i]}; // first endpoint
      TVector3 v2 {muontrk_x2[i],muontrk_y2[i],muontrk_z2[i]}; // second endpoint
      TVector3 vtemp;

      if ( sort_in_z && (v2.Z() < v1.Z()) ) {
	vtemp = v2;
	v2 = v1;
	v1 = vtemp;
      }

      /**
	 now we calculate our own thetaxz and thetayz
      */

      double dx = v2.X() - v1.X();
      double dy = v2.Y() - v1.Y();
      double dz = v2.Z() - v1.Z();
      double theta_xz = atan2(dx,dz)*TMath::RadToDeg();
      double theta_yz = atan2(dy,dz)*TMath::RadToDeg();

      
      if ( sort_in_z && !modify_angles ) {
	if ( muontrk_tpc[i] == 0 && theta_xz < 0 ) {
	  theta_xz = 180+theta_xz;
	}
	if ( muontrk_tpc[i] == 1 ) {
	  if ( theta_xz > 0 ) theta_xz = 180-theta_xz;
	  if ( theta_xz < 0 ) theta_xz = -theta_xz;
	}
      }
      double modified_thetaxz = muontrk_th_xz[i];
      double modified_thetayz = muontrk_th_yz[i];
      if ( sort_in_z && modify_angles ) {
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
	h_modified_thetaxz->Fill(modified_thetaxz);
	h_modified_thetayz->Fill(modified_thetayz);
	std::cout << "regular and modified" << std::endl;
	std::cout << muontrk_th_xz[i] << " " << modified_thetaxz << std::endl;
      }
      
      if ( !sort_in_z && !uncorrected ) {
	if ( muontrk_tpc[i] == 1 ) theta_xz = -theta_xz;
      }
      //if ( theta_xz < 0 ) theta_xz = -theta_xz; // effectively same as previous line
      h_own_thetayz->Fill(theta_yz);
      h_own_thetaxz->Fill(theta_xz);
      // they should be matched at this point
      // BAD MATCHES
      if ( 
	   std::abs(theta_xz - muontrk_th_xz[i]) > match_threshold
	   || std::abs(theta_yz - muontrk_th_yz[i]) > match_threshold
	  ) {

	/*
	  std::cout << Form("event %i we have muontrack type %i: \n",events,muontrk_type[i]);
	  std::cout << muontrk_x1[i] << "\t" << muontrk_y1[i] << "\t" << muontrk_z1[i] << std::endl;
	  std::cout << muontrk_x2[i] << "\t" << muontrk_y2[i] << "\t" << muontrk_z2[i] << std::endl;
	  std::cout << "theta xz, yz own : " << theta_xz << " " << theta_yz << std::endl;
	  std::cout << "theta xz, yz ntu : " << muontrk_th_xz[i] << " " << muontrk_th_yz[i] << std::endl;
	  std::cout << "-------" << std::endl;
	*/
	bad_match++;
	//std::cout << "bad: " << theta_xz << " " << muontrk_th_xz[i] << std::endl;
	//std::cout << "bad: " << theta_yz << " " << muontrk_th_yz[i] << std::endl;
	//v1.Print();
	//v2.Print();

      } else { // GOOD MAATCHES
	good_match++;
	//std::cout << theta_xz << " " << muontrk_th_xz[i] << std::endl;
	//std::cout << theta_yz << " " << muontrk_th_yz[i] << std::endl;

      }
    }
  }
  std::cout << "Muon track type: " << type_to_check << std::endl;
  std::cout << "Match threshold: " << match_threshold << " (deg)" << std::endl;
  std::cout << "Good matches: " << good_match << std::endl;
  std::cout << "Bad matches: " << bad_match << std::endl;

  // draw histograms
  TCanvas *c1 = new TCanvas();
  TLatex *tl_muontype = new TLatex();
  tl_muontype->SetTextFont(132);
    
  TLegend *myleg = new TLegend(0.5,0.7,0.95,0.95);
  h_own_thetaxz->SetMarkerStyle(kFullSquare);
  h_own_thetayz->SetMarkerStyle(kFullSquare);
  h_ntuple_thetaxz->SetLineColor(kBlue);
  h_ntuple_thetayz->SetLineColor(kBlue);
  h_ntuple_thetaxz->SetLineWidth(2);
  h_ntuple_thetayz->SetLineWidth(2);
  
  
  double max = std::max(h_own_thetaxz->GetMaximum(),h_ntuple_thetaxz->GetMaximum());
  TString filename_modifier;
  h_ntuple_thetaxz->SetMaximum(max*1.5);
  h_ntuple_thetaxz->Draw();
  h_own_thetaxz->Draw("SAME P");
  if ( modify_angles ) {
    h_modified_thetaxz->Draw("SAME HIST");
    h_modified_thetaxz->SetLineColor(kRed);
    h_modified_thetaxz->SetLineWidth(2);
    myleg->AddEntry(h_modified_thetaxz,Form("Modified muontracks (%.0f)",h_modified_thetaxz->GetEntries()));
    filename_modifier = "_modified";
  }
  myleg->AddEntry(h_own_thetaxz, Form("Own implementation (%.0f)",h_own_thetaxz->GetEntries()));
  myleg->AddEntry(h_ntuple_thetaxz, Form("Commissioning ntuples (%.0f)",h_ntuple_thetaxz->GetEntries()));
  myleg->Draw("same");
  draw_muontype_text(tl_muontype,type_to_check);
  c1->SaveAs(Form("img/match_validation/h_thetaxz_type%i%s.png",type_to_check,filename_modifier.Data()));
  c1->SaveAs(Form("img/match_validation/h_thetaxz_type%i%s.pdf",type_to_check,filename_modifier.Data()));

  myleg->Clear();
  max = std::max(h_own_thetayz->GetMaximum(),h_ntuple_thetayz->GetMaximum());
  h_ntuple_thetayz->SetMaximum(max*1.5);
  h_ntuple_thetayz->Draw();
  h_own_thetayz->Draw("SAME P");
  if ( modify_angles ) {
    h_modified_thetayz->Draw("SAME HIST");
    h_modified_thetayz->SetLineColor(kRed);
    h_modified_thetayz->SetLineWidth(2);
    myleg->AddEntry(h_modified_thetayz,Form("Modified muontracks (%.0f)",h_modified_thetayz->GetEntries()));
    filename_modifier = "_modified";
  }
  myleg->AddEntry(h_own_thetayz, Form("Own implementation (%.0f)",h_own_thetayz->GetEntries()));
  myleg->AddEntry(h_ntuple_thetayz, Form("Commissioning ntuples (%.0f)",h_ntuple_thetayz->GetEntries()));
  myleg->Draw("same");
  draw_muontype_text(tl_muontype,type_to_check);
  c1->SaveAs(Form("img/match_validation/h_thetayz_type%i%s.png",type_to_check,filename_modifier.Data()));
  c1->SaveAs(Form("img/match_validation/h_thetayz_type%i%s.pdf",type_to_check,filename_modifier.Data()));
  
}
