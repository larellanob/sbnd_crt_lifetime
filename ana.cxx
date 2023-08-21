void ana()
{
  TString filename {"muon_hitdumper_NS.root"};
  TFile *f = new TFile(filename);
  TTreeReader reader("hitdumper/hitdumpertree",f);

  // evt run subrun
  TTreeReaderValue<int> run (reader,"run");
  TTreeReaderValue<int> subrun (reader,"subrun");
  TTreeReaderValue<int> evt (reader,"event");

  
  // muon tracks
  TTreeReaderArray<float> muon_z1(reader,"muontrk_z1");
  TTreeReaderArray<float> muon_z2(reader,"muontrk_z2");
  TTreeReaderArray<int> muon_type(reader,"muontrk_type");
  TTreeReaderArray<float> muon_th_xz(reader,"muontrk_theta_xz");
  TTreeReaderArray<float> muon_th_yz(reader,"muontrk_theta_yz");
  TTreeReaderArray<int> muon_tpc(reader,"muontrk_tpc");

  // muon hits
  TTreeReaderArray<int> mhit_trk(reader,"mhit_trk"); // compare to nmuontrks
  TTreeReaderArray<int> mhit_tpc(reader,"mhit_tpc");
  TTreeReaderArray<int> mhit_plane(reader,"mhit_plane");
  TTreeReaderArray<int> mhit_wire(reader,"mhit_wire");
  TTreeReaderArray<double> mhit_peakT(reader,"mhit_peakT");
  TTreeReaderArray<double> mhit_charge(reader,"mhit_charge");
  

  // crt hits
  TTreeReaderArray<double> chit_z(reader,"chit_z");
  TTreeReaderArray<double> chit_t(reader,"chit_time");
  TTreeReaderArray<int> chit_plane(reader,"chit_plane");

  // crt tracks
  TTreeReaderArray<double> ct_z1(reader,"ct_z1");
  TTreeReaderArray<double> ct_z2(reader,"ct_z2");
  TTreeReaderArray<double> ct_t(reader,"ct_time");
  
  int events = 0;
  int crossing_counter = 0;

  double eff = 0;
  double pur = 0;

  int TP = 0;
  int FP = 0;
  int FN = 0;

  int total_matches = 0;

  TH1F h_xzangle("h_xzangle",";Angle #theta_{XZ} (deg); Entries",30, 0, 30);
  TH1F h_yzangle("h_yzangle",";Angle #theta_{YZ} (deg); Entries",360, -180, 180);
  TH2F h_time_charge_tpc0("h_time_charge_tpc0",
			  ";Median drift time (TPC 0);Median charge",
			  20,0,4000,100,0,1500);
  TH2F h_time_charge_tpc1("h_time_charge_tpc1",
			  ";Median drift time (TPC 1);Median charge",
			  20,0,4000,100,0,1500);
  
  while ( reader.Next() ) {


    bool b_crt = false;
    bool b_muon = false;

    std::vector<int> i_match;
    
    // crt tracks loop
    for ( int i = 0; i < ct_z1.GetSize(); i++ ) {
      double ct_up   = std::min(ct_z1[i],ct_z2[i]); // upstream
      double ct_down = std::max(ct_z1[i],ct_z2[i]); // downstream
      
      if ( ct_up < -170.0 && ct_down > 770.0 ) {
	//std::cout << "CRT HITS: " << ct_up << " " << ct_down << std::endl;
	b_crt = true;
      }
    } // end crt tracks loop
    
    
    // muon tracks loop
    bool found_type_5 = false;
    for ( int i = 0; i < muon_z1.GetSize(); i++ ) {
      /*
	type 4 muon tracks (up-downstream)
	
	we can use this type of track as our "truth", we would ideally
	identify all of these using our crt selection, so we have two
	metrics:

	1. efficiency

 	   eff = (TP)/(TP+FN)

	 2. purity

	   pur = (TP)/(TP+FP)

	 TP: True positives: events which we have both crt track and muon track
	 FN: False negatives: we have muon track but no crt track
	 FP: False positive: we have crt track but no muon track
       */
      if ( muon_type[i] == 4 ) {
	b_muon = true;
	//std::cout << Form("Muon type 4: run %i subrun %i event %i",*run,*subrun,*evt) << std::endl;
	if ( b_crt ) {
	  i_match.push_back(i);
	  total_matches++;
	  //std::cout << "match for i: " << i << std::endl;
	  //std::cout << "total matches: " << total_matches << std::endl;
	}
      }
      
      /*
	type 5 muon tracks (other)
	
	question we want to answer is: what fraction of type-5 muon
	tracks hit the edges and thus could potentially confuse our
	crt selection
       */
      /*
      if ( muon_type[i] == 5 ) {
	if ( !found_type_5 ) {
	  std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
	  std::cout << "Event: " << events << std::endl;
	  std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
	  found_type_5 = true;
	}
	std::cout << muon_z1[i] << " " << muon_z2[i] << std::endl;
	crossing_counter++;
      }
      */
      
    } // end muon tracks loop

    // now we've looped we can calculate eff and pur
    if ( b_crt and b_muon ) {
      TP++;
    } else if ( b_crt and !b_muon ) {
      FP++;
    } else if ( !b_crt and b_muon ) {
      FN++;
    }
    
    if ( crossing_counter > 50 ) break;
    events++;


    // now if it's a match we can enter a loop with the i_match index

    if ( b_crt and b_muon ) {
      // loop over matched i_match vector
      for ( int m = 0; m < i_match.size(); m++ ) {
	//std::cout << "Event, muon: " << *evt << " " <<  i_match[m] << std::endl;
	int mi = i_match[m]; // muon index
	float xzangle = muon_th_xz[mi];
	float yzangle = muon_th_yz[mi];
	if ( xzangle > 90 ) {
	  xzangle = 180. - xzangle;
	}
	//std::cout << Form("Run, Event: %i, %i. Angle xz for match: %.3f",*run,*evt,xzangle) << std::endl;
	h_xzangle.Fill(xzangle);
	h_yzangle.Fill(yzangle);
	if (xzangle > 5 ) {
	  //std::cout << "Discarding due to large angle" << std::endl;
	  continue;
	}
	// loop over muon hits and select those with matching i_match[mi]
	TH2F h_wire_charge_tpc0("h_wire_charge_tpc0",
				"Run 2, Event 1806;Wire (plane 2);Charge",
				500,0,2000,2000,0,5000);
	TH1F h_charge_tpc0("h_charge_tpc0",
			   "Run 2, Event 1806;Charge;Entries",
			   1000,0,5000);
	TH1F h_time_tpc0("h_time_tpc0",
			 "Run 2, Event 5577;Peak Time;Entries",
			 100,0,4000);
	TH2F h_wire_charge_tpc1("h_wire_charge_tpc1",
				"Run 2, Event 5577;Wire (plane 2);Charge",
				500,0,2000,2000,0,5000);
	TH1F h_charge_tpc1("h_charge_tpc1",
			   "Run 2, Event 5577;Charge;Entries",
			   1000,0,5000);
	TH1F h_time_tpc1("h_time_tpc1",
			 "Run 2, Event 5577;Peak Time;Entries",
			 100,0,4000);
	for ( int mh = 0; mh < mhit_trk.GetSize(); mh++ ) {
	  //select those with matching i_match[m]
	  if ( mhit_trk[mh] != mi ) {
	    continue;
	  }
	  if ( mhit_plane[mh] != 2 ) {
	    continue;
	  }
	  if ( mhit_tpc[mh] == 0 ) {
	    h_wire_charge_tpc0.Fill(mhit_wire[mh],mhit_charge[mh]);
	    h_charge_tpc0.Fill(mhit_charge[mh]);
	    h_time_tpc0.Fill(mhit_peakT[mh]);
	  }
	  if ( mhit_tpc[mh] == 1 ) {
	    h_wire_charge_tpc1.Fill(mhit_wire[mh],mhit_charge[mh]);
	    h_charge_tpc1.Fill(mhit_charge[mh]);
	    h_time_tpc1.Fill(mhit_peakT[mh]);
	  }
	}
	// median
	double q;
	double median_charge_tpc0 = 0;
	double median_time_tpc0   = 0;
	double median_charge_tpc1 = 0;
	double median_time_tpc1   = 0;
	q = 0.5; // 0.5 quantile for median
	if ( h_charge_tpc0.Integral() != 0 ) {
	  h_charge_tpc0.GetQuantiles(1,&median_charge_tpc0,&q);
	}
	if ( h_time_tpc0.Integral() != 0 ) {
	  h_time_tpc0.GetQuantiles(1,&median_time_tpc0,&q);
	}
	if ( h_charge_tpc1.Integral() != 0 ) {
	  h_charge_tpc1.GetQuantiles(1,&median_charge_tpc1,&q);
	}
	if ( h_time_tpc1.Integral() != 0 ) {
	  h_time_tpc1.GetQuantiles(1,&median_time_tpc1,&q);
	}
	if ( median_charge_tpc0 > 800 ) {
	  std::cout << "-------------------------------" << std::endl;
	  std::cout << Form("Run, event, median time: %i, %i, %.3f",*run,*evt,median_time_tpc0) << std::endl;
	  std::cout << "large median charge (TPC0): " << median_charge_tpc0 << std::endl;
	}

	if ( median_time_tpc0 != 0 && median_charge_tpc0 != 0 ) {
	  h_time_charge_tpc0.Fill(median_time_tpc0,median_charge_tpc0);
	}
	if ( median_time_tpc1 != 0 && median_charge_tpc1 != 0 ) {
	  h_time_charge_tpc1.Fill(median_time_tpc1,median_charge_tpc1);
	}
	// particular examples
	if ( *run == 2 && *evt == 5577 ) {
	  TLatex * text_median = new TLatex();
	  text_median->SetTextSize(0.06);
	  text_median->SetTextFont(132);
	  text_median->SetTextAlign(33);
	  auto * c2 = new TCanvas();

	  // TPC1
	  h_wire_charge_tpc1.Draw("colz");
	  c2->SaveAs("img/dqdx/2d_wire_charge_tpc1.pdf");
	  c2->SaveAs("img/dqdx/2d_wire_charge_tpc1.png");
	  h_charge_tpc1.Draw();
	  text_median->DrawLatexNDC(.8,.8,Form("Median: %.3f",median_charge_tpc1));
	  c2->SaveAs("img/dqdx/1d_charge_tpc1.pdf");
	  c2->SaveAs("img/dqdx/1d_charge_tpc1.png");
	  h_time_tpc1.Draw();
	  text_median->DrawLatexNDC(.8,.8,Form("Median: %.3f",median_time_tpc1));
	  c2->SaveAs("img/dqdx/1d_time_tpc1.pdf");
	  c2->SaveAs("img/dqdx/1d_time_tpc1.png");


	  delete c2;
	  delete text_median;
	}
	if ( *run == 2 && *evt == 1806 ) {
	  TLatex * text_median = new TLatex();
	  text_median->SetTextSize(0.06);
	  text_median->SetTextFont(132);
	  text_median->SetTextAlign(33);
	  auto * c2 = new TCanvas();
	  
	  // TPC0
	  h_wire_charge_tpc0.Draw("colz");
	  c2->SaveAs("img/dqdx/2d_wire_charge_tpc0.pdf");
	  c2->SaveAs("img/dqdx/2d_wire_charge_tpc0.png");
	  h_charge_tpc0.Draw();
	  text_median->DrawLatexNDC(.8,.8,Form("Median: %.3f",median_charge_tpc0));
	  c2->SaveAs("img/dqdx/1d_charge_tpc0.pdf");
	  c2->SaveAs("img/dqdx/1d_charge_tpc0.png");
	  h_time_tpc0.Draw();
	  text_median->DrawLatexNDC(.8,.8,Form("Median: %.3f",median_time_tpc0));
	  c2->SaveAs("img/dqdx/1d_time_tpc0.pdf");
	  c2->SaveAs("img/dqdx/1d_time_tpc0.png");

	  delete c2;
	  delete text_median;
	}

	
      }
    }
  }
  eff = (double)TP/( (double)TP + (double)FN );
  pur = (double)TP/( (double)TP + (double)FP );

  std::cout << Form("Looped over %i events", events) << std::endl;
  std::cout << Form("TP: %i \t FP: %i \t FN: %i",TP,FP,FN) << std::endl;
  std::cout << Form("EFF: %.3f, PUR: %.3f",eff,pur) << std::endl;

  auto * c1 = new TCanvas();
  h_xzangle.Draw();
  c1->SaveAs("img/dqdx/xzangle.pdf");
  c1->SaveAs("img/dqdx/xzangle.png");
  
  h_time_charge_tpc0.Draw("colz");
  c1->SaveAs("img/dqdx/2d_time_charge_tpc0.pdf");
  c1->SaveAs("img/dqdx/2d_time_charge_tpc0.png");
  h_time_charge_tpc1.Draw("colz");
  c1->SaveAs("img/dqdx/2d_time_charge_tpc1.pdf");
  c1->SaveAs("img/dqdx/2d_time_charge_tpc1.png");
  //delete c3;
  
}
