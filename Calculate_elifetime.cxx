void Calculate_elifetime()
{
  // open file
  TString filename {"CORRECTED_muon_hitdumper_NS.root"};
  TFile *f = new TFile(filename);

  TTreeReader reader("hitdumper/hitdumpertree_corr",f);
  // load only the corrected muon tracks and muon hits
  // muon tracks
  TTreeReaderArray<double> muontrk_x1(reader,"muontrk_x1");
  TTreeReaderArray<double> muontrk_y1(reader,"muontrk_y1");
  TTreeReaderArray<double> muontrk_z1(reader,"muontrk_z1");
  TTreeReaderArray<double> muontrk_x2(reader,"muontrk_x2");
  TTreeReaderArray<double> muontrk_y2(reader,"muontrk_y2");
  TTreeReaderArray<double> muontrk_z2(reader,"muontrk_z2");
  //muon hits
  TTreeReaderArray<double> mhit_charge(reader,"mhit_charge");
  TTreeReaderArray<int> mhit_wire(reader,"mhit_wire");
  
  // angles
  TTreeReaderArray<double> theta_xz(reader,"tpc_thetaxz");
  TTreeReaderArray<double> theta_yz(reader,"tpc_thetayz");
  // misc
  TTreeReaderValue<int> ev(reader,"event");

  // figure out number of front-back crossing muon tracks per bin for
  // different x-binnings
  int nbins1 = 10;
  int nbins2 = 20;
  int nbins3 = 30;
  int nbins4 = 40;
  int nbins5 = 50;

  double safe_angle1 = 10.0;
  double safe_angle2 = 5.0;
  double safe_angle3 = 3.0;
  double safe_angle4 = 2.5;
  double safe_angle5 = 2.0;

  std::vector<std::pair<int,double>> v_binning;
  v_binning.push_back(std::make_pair(nbins1,safe_angle1));
  v_binning.push_back(std::make_pair(nbins2,safe_angle2));
  v_binning.push_back(std::make_pair(nbins3,safe_angle3));
  v_binning.push_back(std::make_pair(nbins4,safe_angle4));
  v_binning.push_back(std::make_pair(nbins5,safe_angle5));

  int binning = 2; // binning 2 -> 30 bins
  int nbins = v_binning[binning].first;

  double bin_width_x = (202.05-(-202.05))/nbins;
  double drift_vel = 160.; // 0.16 cm/us = 160 cm/ms
  double bin_width_t = ((202.05-(-202.05))/drift_vel)/nbins;
  double dz = 0.3; // 3mm wire pitch
  
  TH1F *h_binning1 = new TH1F("h_binning1",
			      Form("%i bins, angle selection;x-coordinate;Number of matched tracks",nbins),
			      nbins,
			      -202.05,
			      202.05);
  
  std::cout << "Bin width in x: " << bin_width_x << " cm" << std::endl;
  std::cout << "Bin width in t: " << bin_width_t << " ms" << std::endl;

  std::vector<TH1F*> v_charge_per_bin;
  std::vector<TH1F*> v_dQdx;
  std::vector<int> v_tracks_per_bin;
  for ( int i = 0; i < nbins; i++ ) {
    v_charge_per_bin.push_back(new TH1F(Form("h_xbin%i",i+1),
					Form("x bin %i (%.2f < x < %.2f);Charge per hit;Number of hits",
					     i+1,-202.05+i*bin_width_x, -202.05+(i+1)*bin_width_x),
					1000,
					0,
					5000)
			       );
    v_dQdx.push_back(new TH1F(Form("h_dqdx_xbin%i",i+1),
			      Form("x bin %i (%.2f < x < %.2f);dQ/dx (ADC/cm);Number of hits",
				   i+1,-202.05+i*bin_width_x, -202.05+(i+1)*bin_width_x),
			      100,
			      0,
			      5000)
		     );
    
    v_tracks_per_bin.push_back(0);
  }
  

  // read the file
  while ( reader.Next() ) {
    // angle cut
    for ( int i = 0; i < muontrk_x1.GetSize(); i++ ) {
      // filter out tracks not straight enough
      if ( theta_xz[i] < -v_binning[binning].second || theta_xz[i] > v_binning[binning].second ) {
	continue;
      }
      // filter out tracks not crossing enough
      /*
	if ( muontrk_z1[i] > 15 || muontrk_z2[i] < 480 ) {
	continue;
	}
      */
      // figure out bin
      // root's Fill actually returns the bin number
      int bin_filled = h_binning1->Fill((muontrk_x1[i]+muontrk_x2[i])/2.0);
      //std::cout << bin_filled << std::endl;
      if ( bin_filled < 1 || bin_filled > nbins ) continue;
      v_tracks_per_bin[bin_filled-1] = v_tracks_per_bin[bin_filled-1]+=1;
      // dQ/dx, let's call it dQ/dr to avoid confusion
      double dr = dz*sqrt(1.0 +
			  tan(theta_xz[i]*TMath::DegToRad())*tan(theta_xz[i]*TMath::DegToRad()) +
			  tan(theta_yz[i]*TMath::DegToRad())*tan(theta_yz[i]*TMath::DegToRad())
			  );
      //std::cout << "dz, theta_xz, theta_yz: " << dz << " " << theta_xz << " " << theta_yz << std::endl;
      std::cout << "dr = " << dr << std::endl;
      for ( int ch = 0; ch < mhit_charge.GetSize(); ch++ ) {
	v_charge_per_bin[bin_filled-1]->Fill(mhit_charge[ch]);
	v_dQdx[bin_filled-1]->Fill(mhit_charge[ch]/dr);
	
      }
    } // end looping over event muon tracks
  } // end looping over events


  // I have the histograms. let's draw, fit, and so on

  
  TCanvas *c1 = new TCanvas();
  h_binning1->Draw();
  c1->SaveAs(Form("img/elifetime/binning%i.pdf",binning));
  c1->SaveAs(Form("img/elifetime/binning%i.png",binning));

  TLatex * l_ntracks = new TLatex();
  TLatex * l_lpg = new TLatex();
  TLatex * l_par = new TLatex();
  TH1F *h_elifetime_tpc0 = new TH1F("h_elifetime_tpc0",
				    "Drift time",
				    nbins/2,
				    0,
				    bin_width_t*nbins/2);
  TH1F *h_elifetime_tpc1 = new TH1F("h_elifetime_tpc1",
				    "Drift time",
				    nbins/2,
				    0,
				    bin_width_t*nbins/2);
  int tpc = 0;
  double drift_time = -999.;
  TLegend *myleg = new TLegend(0.65,0.3,0.9,0.6);
				   
  for ( int i = 0; i < nbins; i++ ) {
  //for ( int i = 0; i < 5; i++ ) {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Bin: " << i+1 << std::endl;
    std::cout  << std::endl;
    std::cout  << std::endl;
    if ( i >= 15 ) tpc = 1;
    double par[6] = {200.0,500.0,80.0,400.0,500.0,50.0};
    // tf1convolution ("f1", "f2", xmin, xmax, bool useFFT)


    double center_ini = 1500;

    double low_conv = 500.0;
    double high_conv = 2250.0;

    double low_gaus = 2250.0;
    double high_gaus= 3500.0;
    
    TF1Convolution * LconvG = new TF1Convolution("landau","gaus",low_conv,high_conv,false);
    TF1 * g = new TF1("g","gaus",low_gaus,high_gaus);
    TF1 * f_LconvG = new TF1("f_LconvG",*LconvG,low_conv,high_conv,LconvG->GetNpar());
    //g->SetParameters(50,center_ini*2,50);
    //TF1 *fit_full = new TF1("fit_full","*LconvG+g",low_conv,high_conv);
    /*
    auto lpg = new TF1("lpg","landau(0)+gaus(3)",low_conv,800);
    lpg->SetParameters(par);
    v_charge_per_bin[i]->Fit("lpg","","",low_conv,high_conv);
    */
    double pars[6] = {100.0,center_ini,50.0,
		      100.0,center_ini,50.0};
    double pars2[3] = { 100.0,center_ini*2,50.0
    };
    //fit_full->SetParameters(pars);
    /*
    f_LconvG->SetParameters(pars);
    g->SetParameters(pars2);
    f_LconvG->SetLineColor(kGreen+1);
    g->SetLineColor(kBlue);
    f_LconvG->SetLineWidth(3);
    g->SetLineWidth(3);
    v_charge_per_bin[i]->Draw();
    v_charge_per_bin[i]->Fit("f_LconvG","","",low_conv,high_conv);
    v_charge_per_bin[i]->Fit("g","","",low_gaus,high_gaus);
    f_LconvG->Draw("same");
    //LconvG->Draw();
    //g->Draw();
    //l_lpg->DrawLatexNDC(0.5,0.6,"G(C_{G},#mu_{G},#sigma_{G})+L(C_{L},#mu_{L},#sigma_{L})");
    l_ntracks->DrawLatexNDC(0.35,0.85,Form("N tracks: %i",v_tracks_per_bin[i]));
    l_lpg->DrawLatexNDC(0.35,0.7,"C_{G}exp#left(#frac{x-#mu_{G}}{#sigma_{G}}#right)+C_{L} * L(#mu_{L},#sigma_{L})");
    myleg->Clear();
    myleg->AddEntry(f_LconvG,"L(x) #otimes G_{1}(x)");
    myleg->AddEntry(g,"G_{2}(x)");
    myleg->Draw("same");
    */
    /*
    l_par->DrawLatexNDC(0.35,0.5,
			Form("#splitline{C_{G}: %.2f, #mu_{G}: %.2f, #sigma_{G}: %.2f}{C_{L}:%.2f, #mu_{L}: %.2f, #sigma_{L}: %.2f}",
			     lpg->GetParameter(0),
			     lpg->GetParameter(1),
			     lpg->GetParameter(2),
			     lpg->GetParameter(3),
			     lpg->GetParameter(4),
			     lpg->GetParameter(5)));
    */
    /*
    l_par->DrawLatexNDC(0.35,0.5,
			Form("#splitline{C_{G}: %.2f, #mu_{G}: %.2f, #sigma_{G}: %.2f}{C_{L}:%.2f, #mu_{L}: %.2f, #sigma_{L}: %.2f}",
			     fit_full->GetParameter(0),
			     fit_full->GetParameter(1),
			     fit_full->GetParameter(2),
			     fit_full->GetParameter(3),
			     fit_full->GetParameter(4),
			     fit_full->GetParameter(5)));
    */

    // moved to dQdx, no dQ/d wire now
    //c1->SaveAs(Form("img/elifetime/charge/charge_binning%i_bin%i.pdf",binning,i+1));
    //c1->SaveAs(Form("img/elifetime/charge/charge_binning%i_bin%i.png",binning,i+1));


    // dQdx proper
    v_dQdx[i]->Draw();
    l_ntracks->DrawLatexNDC(0.4,0.8,Form("N tracks: %i",v_tracks_per_bin[i]));

    center_ini = 1500;
    auto lpg = new TF1("lpg","gaus",low_conv,high_conv);
    double dQdx_landau[3] = { 1000.0, center_ini, 150.0};
    lpg->SetParameters(dQdx_landau);
    /*    double dQdx_pars[6] = {1000.0,center_ini,150.0,
			   150.0,center_ini,250.0};
    */
     double dQdx_pars[5] = {1000.0,center_ini,150.0,
			   center_ini,250.0};

    double dQdx_pars2[3] = { 100.0,center_ini*2,50.0
    };
    f_LconvG->SetParameters(dQdx_pars);
    g->SetParameters(dQdx_pars2);
    f_LconvG->SetLineColor(kGreen+1);
    g->SetLineColor(kBlue);
    f_LconvG->SetLineWidth(3);
    g->SetLineWidth(3);
    //v_dQdx[i]->Fit("lpg","","",low_conv,high_conv);
    v_dQdx[i]->Fit("f_LconvG","","",low_conv,high_conv);
    v_dQdx[i]->Fit("g","","",low_gaus,high_gaus);
    //lpg->Draw("same");
    f_LconvG->Draw("same");
    
    myleg->Clear();
    myleg->AddEntry(f_LconvG,"L(x) #otimes G_{1}(x)");
    myleg->AddEntry(g,"G_{2}(x)");
    myleg->Draw("same");
    l_par->DrawLatexNDC(0.4,0.6,Form("#splitline{#mu_{L} = %.2f, #mu_{G_{1}} = %.2f}{#mu_{L}+#mu_{G_{1}} = %.2f}",
				       f_LconvG->GetParameter(1),
				       f_LconvG->GetParameter(3),
				       f_LconvG->GetParameter(1)+f_LconvG->GetParameter(3)
				       )
			);
    
    c1->SaveAs(Form("img/elifetime/charge/dQdx_binning%i_bin%i.pdf",binning,i+1));
    c1->SaveAs(Form("img/elifetime/charge/dQdx_binning%i_bin%i.png",binning,i+1));
    

    
    /*
    if ( tpc == 0 ) {
      drift_time = i*bin_width_t;
      h_elifetime_tpc0->SetBinContent(i,fit_full->GetParameter(4));
    }
    if ( tpc == 1 ) {
      drift_time = (nbins*bin_width_t)-(nbins-i)*bin_width_t;
      std::cout << drift_time << std::endl;
      h_elifetime_tpc1->SetBinContent(i,fit_full->GetParameter(4));
    }
    */
  }
  gStyle->SetOptStat(0);
  h_elifetime_tpc0->Draw();
  c1->SaveAs("img/elifetime/elifetime_tpc0.pdf");
  c1->SaveAs("img/elifetime/elifetime_tpc0.png");
  h_elifetime_tpc1->Draw();
  c1->SaveAs("img/elifetime/elifetime_tpc1.pdf");
  c1->SaveAs("img/elifetime/elifetime_tpc1.png");
  
}
