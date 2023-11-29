double langaufun(double *x, double *par) {
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

void draw_x_bins(int index, std::vector<std::pair<TVector3,TVector3>> tracks, TH1F* binning, TString modifier = "")
{
  // histo. bins start from 1, vector index starts from 0
  int bin = index+1;

  // base histogram SBND (z,x) coordinates
  TH2F h_tpc_xz_trackbins("h_tpc_xz_trackbins",
			  Form("Bin %i (%.2f < #it{x} < %.2f), %i tracks;#it{z} coordinate (cm);#it{x} coordinate (cm)",
			       bin,binning->GetBinLowEdge(bin),binning->GetBinLowEdge(bin+1),(int)tracks.size()),
			  100,5.25,504.15,
			  30,-202.05,202.05);

  // draw "tpc" 2d histo and bin lines
  TCanvas ctemp("c","c",1200,1000);
  ctemp.SetTopMargin(0.1);
  h_tpc_xz_trackbins.Draw();
  int nbins = binning->GetNbinsX();
  for ( int i = 0; i < nbins; i++ ) {
    double bin_lowedge = binning->GetBinLowEdge(i+1);
    double bin_upedge = binning->GetBinLowEdge(i+2);
    TLine *l = new TLine(5.25,bin_lowedge,504.15,bin_lowedge);
    if ( bin_lowedge == 0.0 ) {
      l->SetLineWidth(3);
    }
    l->Draw("same");
  }

  // highlight current bin using tbox with transparency
  double bin_lowedge = binning->GetBinLowEdge(bin);
  double bin_upedge = binning->GetBinLowEdge(bin+1);
  TBox *b = new TBox(5.25,bin_lowedge,504.15,bin_upedge);
  b->SetFillColorAlpha(kBlue, 0.25);
  b->Draw("same");

  // draw track lines
  for ( int i = 0; i < tracks.size(); i++ ) {
    TLine *l = new TLine(tracks[i].first.Z(),tracks[i].first.X(),
			 tracks[i].second.Z(),tracks[i].second.X());
    l->Draw("same");
    l->SetLineColor(kRed);
    l->SetLineWidth(2);
  }

  // save canvas to file
  gSystem->mkdir(Form("img/elifetime/%s/tracks_per_xbin",modifier.Data()),true);
  ctemp.SaveAs(Form("img/elifetime/%s/tracks_per_xbin/tracks_bin%i.png",modifier.Data(),bin));
  ctemp.SaveAs(Form("img/elifetime/%s/tracks_per_xbin/tracks_bin%i.pdf",modifier.Data(),bin));
}
  

void Calculate_elifetime(TString modifier = "")
{
  gROOT->SetStyle("uboone_sty");
  // open file
  if ( modifier != "" ) {
    modifier = "_"+modifier;
  }
  TString filename {Form("trees/muon_hitdumper_NS%s_corrected.root",modifier.Data())};
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
  double drift_vel = 160.0; // 0.16 cm/us = 160 cm/ms
  double bin_width_t = ((202.05-(-202.05))/drift_vel)/nbins;
  double dz = 0.3; // 3mm wire pitch
  
  TH1F *h_binning1 = new TH1F("h_binning1",
			      Form("%i bins, |#theta_{XZ}| < %.2f deg;#it{x}-coordinate;Selected matched tracks",
				   nbins,v_binning[binning].second),
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
					Form(";Charge per hit;Number of hits"),
					1000,
					0,
					5000)
			       );
    v_dQdx.push_back(new TH1F(Form("h_dqdx_xbin%i",i+1),
			      Form(";dQ/dx (ADC/cm);Number of hits"),
			      100,
			      0,
			      5000)
		     );
    
    v_tracks_per_bin.push_back(0);
  }
  

  // read the file
  std::vector<std::pair<TVector3,TVector3>> v_binned_muontrks[nbins];
  while ( reader.Next() ) {
    // array of vectors of tracks that go in each bin division
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
      for ( int ch = 0; ch < mhit_charge.GetSize(); ch++ ) {
	v_charge_per_bin[bin_filled-1]->Fill(mhit_charge[ch]);
	v_dQdx[bin_filled-1]->Fill(mhit_charge[ch]/dr);
	
      }
      // draw binned tracks
      TVector3 endpoint1 {muontrk_x1[i],muontrk_y1[i],muontrk_z1[i]};
      TVector3 endpoint2 {muontrk_x2[i],muontrk_y2[i],muontrk_z2[i]};
      TVector3 muontrack = (
			    TVector3{muontrk_x2[i],muontrk_y2[i],muontrk_z2[i]}
			    - TVector3{ muontrk_x1[i],muontrk_y1[i],muontrk_z1[i]});
      v_binned_muontrks[bin_filled-1].push_back(std::make_pair(endpoint1,endpoint2));
    } // end looping over event muon tracks
    
  } // end looping over events

  // I have the histograms. let's draw, fit, and so on

  // prepare image output directory
  // remove leading '_' from filename modifier
  if ( modifier[0] == '_' ) {
    modifier.Remove(0,1); // remove 1 character, starting from char 0
  }
  // mkdir ( name, recursive)
  if ( modifier == "" ) {
    modifier = "default";
  }
  gSystem->mkdir(Form("img/elifetime/%s/charge/raw",modifier.Data()),true);

  // loop over bins and draw tracks
  for ( int i = 0; i < nbins; i++ ) {
    draw_x_bins(i,v_binned_muontrks[i],h_binning1,modifier);
  }
  TCanvas *c1 = new TCanvas();
  h_binning1->Draw();
  c1->SaveAs(Form("img/elifetime/%s/binning%i.pdf",modifier.Data(),binning));
  c1->SaveAs(Form("img/elifetime/%s/binning%i.png",modifier.Data(),binning));

  TLatex * l_xbin = new TLatex();
  TLatex * l_ntracks = new TLatex();
  TLatex * l_par = new TLatex();
  l_xbin->SetTextFont(132);
  l_ntracks->SetTextFont(132);
  l_par->SetTextFont(132);
  TH1F *h_elifetime_tpc0 = new TH1F("h_elifetime_tpc0",
				    ";Drift time (ms); dQ/dx MPV (ADC/cm)",
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
    //if ( i != 7 ) continue;
    //for ( int i = 0; i < 5; i++ ) {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Bin: " << i+1 << std::endl;
    std::cout  << std::endl;
    std::cout  << std::endl;
    if ( i >= 15 ) tpc = 1;
    double par[6] = {200.0,500.0,80.0,400.0,500.0,50.0};
    // tf1convolution ("f1", "f2", xmin, xmax, bool useFFT)


    // look for the peak and grab the x position of it
    double conv_center = v_dQdx[i]->GetBinCenter(v_dQdx[i]->GetMaximumBin());
    double gaus_center = conv_center*2.0;
    double conv_gaus_intersection = conv_center +0.75*(gaus_center-conv_center);


    double low_gaus = conv_gaus_intersection;
    double high_gaus = gaus_center + (gaus_center-conv_gaus_intersection);
    double high_conv = conv_gaus_intersection;
    double low_conv = std::max(0.0,conv_center-0.5*0.75*(gaus_center-conv_center));
    
    /*
    double low_conv = 200.0;
    double high_conv = 2000.0;

    double low_gaus = 2250.0;
    double high_gaus= 3500.0;
    */
    
    TF1 * g = new TF1("g","gaus",low_gaus,high_gaus);

    // dQdx proper
    v_dQdx[i]->SetMaximum(v_dQdx[i]->GetMaximum()*1.2);
    v_dQdx[i]->Draw();


    /** 
	We want the integral of the number of hits, to get a good
	first approximation of the scale of our fit, so we compute the
	integral of the histogram between the fitting range and
	scale. The hardcoded scaling factor of 260.0 was obtained via
	trial and error
    */
    TAxis *xaxis = v_dQdx[i]->GetXaxis();
    double integral_bins = v_dQdx[i]->Integral(xaxis->FindBin(low_conv),xaxis->FindBin(high_conv));
    double scale = integral_bins/260.; // 260 -> trial and error

    // fit for the delta ray gaussian
    double dQdx_pars2[3] = { 100.0,gaus_center,50.0};
    // fit using langaufun
    TF1 * f_langaufun = new TF1("f_langaufun",langaufun,low_conv,high_conv,4);
    // parameters before trying to grab them blindly (28/11/23)
    //double p_langaufun[4] = { 80.0, 1350., scale*15000, 180.};
    double p_langaufun[4] = { 80.0, conv_center, scale*15000, 180.};
    f_langaufun->SetParameters(p_langaufun);
    v_dQdx[i]->Fit("f_langaufun","L","",low_conv,high_conv);
    
    //f_LconvG->Draw("same");
    f_langaufun->SetLineColor(kGreen+1);
    f_langaufun->SetLineWidth(3);
    f_langaufun->Draw("same");

    // delta ray gaussian fitting
    g->SetParameters(dQdx_pars2);
    g->SetLineColor(kBlue);
    g->SetLineWidth(3);
    v_dQdx[i]->Fit("g","L Q","",low_gaus,high_gaus);
    // last function to be fitted doesn't need a ->Draw()
    double p0,p1,p2,p3;
    double fp0,fp1,fp2,fp3;
    p0 = f_langaufun->GetParameter(0);
    p1 = f_langaufun->GetParameter(1);
    p2 = f_langaufun->GetParameter(2);
    p3 = f_langaufun->GetParameter(3);
    
    fp0 = f_langaufun->Eval(p0);
    fp1 = f_langaufun->Eval(p1);
    fp2 = f_langaufun->Eval(p2);
    fp3 = f_langaufun->Eval(p3);
    //TLine *line_mpv = new TLine(p1+p2,0,p1+p2,fp1);
    TLine *line_mpv = new TLine(p1,0,p1,fp1);
    line_mpv->SetLineWidth(2);
    line_mpv->SetLineStyle(9);
    line_mpv->Draw("same");
    
    myleg->Clear();
    myleg->AddEntry(f_langaufun,"L(x) * G_{1}(x)");
    myleg->AddEntry(g,"G_{2}(x)");
    myleg->AddEntry(line_mpv,"Fit L_{MPV}");
    myleg->Draw("same");

    // draw info latex text
    l_xbin->DrawLatexNDC(0.5,0.9,Form("#it{x}-bin %i (%.2f < #it{x} < %.2f)",
				      i+1,-202.05+i*bin_width_x, -202.05+(i+1)*bin_width_x));
    l_ntracks->DrawLatexNDC(0.5,0.83,Form("N tracks: %i, N hits: %i",v_tracks_per_bin[i],(int)v_dQdx[i]->GetEntries()));

    l_par->DrawLatexNDC(0.5,0.68,Form("#splitline{#mu_{L} = %.2f, #mu_{G_{2}} = %.2f}{#chi^{2}_{L} = %.2f, #chi^{2}_{G_{2}} = %.2f}",
				      f_langaufun->GetParameter(1),
				      g->GetParameter(1),
				      f_langaufun->GetChisquare(),
				      g->GetChisquare()
				      )
			);

    // use derivative to grab mpv of convolution instead of deconvolved
    double minderiv = 10000;
    double minx = low_conv-1.01;
    for ( double x = low_conv; x < high_conv; x+=1.0 ) {
      double deriv = f_langaufun->Derivative(x);
      if ( deriv < minderiv ) {
	minderiv = deriv;
	minx = x;
      }
    }
    
    double lpg_mpv = f_langaufun->GetParameter(1);
    double error_lpg_mpv = f_langaufun->GetParError(1);
    std::cout << "MPV AND ERROR: " << lpg_mpv << " " << error_lpg_mpv << std::endl;
    //double lpg_mpv = minx;
    
    c1->SaveAs(Form("img/elifetime/%s/charge/dQdx_binning%i_bin%i.pdf",modifier.Data(),binning,i+1));
    c1->SaveAs(Form("img/elifetime/%s/charge/dQdx_binning%i_bin%i.png",modifier.Data(),binning,i+1));
    
    // draw histograms with "raw" landau and gauss for convolution
    TF1 * f_landau_raw = new TF1("f_landau_raw","landau",low_conv,high_conv);
    double p_landau_raw[3] = {
      f_langaufun->GetParameter(2)*0.008,
      f_langaufun->GetParameter(1),
      f_langaufun->GetParameter(0)
    };
    f_landau_raw->SetParameters(p_landau_raw);
    f_landau_raw->Draw("same");
    TF1 * f_gauss_raw = new TF1("f_gauss_raw","gaus",0.0,low_conv);
    double p_gauss_raw[3] = {
      f_langaufun->GetParameter(2)*0.0005,
      0.0,
      f_langaufun->GetParameter(3)
    };
    f_gauss_raw->SetParameters(p_gauss_raw);
    f_gauss_raw->SetLineColor(kMagenta);
    f_gauss_raw->Draw("same");

    c1->SaveAs(Form("img/elifetime/%s/charge/raw/dQdx_binning%i_bin%i_plusraw.pdf",modifier.Data(),binning,i+1));
    c1->SaveAs(Form("img/elifetime/%s/charge/raw/dQdx_binning%i_bin%i_plusraw.png",modifier.Data(),binning,i+1));
    
    if ( lpg_mpv  > 10000 ) continue;
    if ( tpc == 0 ) {
      drift_time = i*bin_width_t;
      h_elifetime_tpc0->SetBinContent(i+1,lpg_mpv);
      h_elifetime_tpc0->SetBinError(i+1,error_lpg_mpv);
    }
    if ( tpc == 1 ) {
      drift_time = (nbins*bin_width_t)-(nbins-i)*bin_width_t;
      h_elifetime_tpc1->SetBinContent(16-(i-14),lpg_mpv);
      h_elifetime_tpc1->SetBinError(16-(i-14),error_lpg_mpv);
    }

  }
  gStyle->SetOptStat(0);

  TLegend *myleg2 = new TLegend(0.5,0.7,0.9,0.95);
  double bin_end =  bin_width_t*nbins/2;
  std::cout << "bin end: " << bin_end << std::endl;

  double min_elifetime_histo = std::min(h_elifetime_tpc0->GetMinimum(),h_elifetime_tpc1->GetMinimum());
  double max_elifetime_histo = std::max(h_elifetime_tpc0->GetMaximum(),h_elifetime_tpc1->GetMaximum());
    
  h_elifetime_tpc0->SetMinimum(std::max(min_elifetime_histo-0.02*min_elifetime_histo,0.0));
  h_elifetime_tpc0->SetMaximum(max_elifetime_histo+0.1*max_elifetime_histo);

  h_elifetime_tpc0->Draw("E1");
  TF1 * e0 = new TF1("e0","[0]*exp(-x/[1])",0.0,1.3);
  double p_e0[2] = { 1350, 10 };
  TF1 * e0_default = new TF1("e0_default","[0]*exp(-x/[1])",0.0,1.3);
  e0_default->SetParameters(p_e0);
  e0_default->SetLineColor(kOrange);
  e0->SetParameters(p_e0);
  //e0->FixParameter(1,td);
  //h_elifetime_tpc0->Fit(e0,"L N","",0.0,bin_end);
  h_elifetime_tpc0->Fit(e0,"L N Q","",0.0,1.0);
  e0->SetLineColor(kBlue);
  e0->SetLineStyle(2);
  e0->SetLineWidth(2);
  e0->Draw("same");
  //e0->Draw();
  //e0_default->Draw("same");
  
  c1->SaveAs(Form("img/elifetime/%s/elifetime_tpc0.pdf",modifier.Data()));
  c1->SaveAs(Form("img/elifetime/%s/elifetime_tpc0.png",modifier.Data()));
  //h_elifetime_tpc1->SetMinimum(1100);
  h_elifetime_tpc1->SetLineColor(kRed);
  h_elifetime_tpc1->Draw("same E1");
  TF1 * e1 = new TF1("e1","[0]*exp(-x/[1])",0.0,1.3);
  double p_e1[2] = { 1350, 10 };
  TF1 * e1_default = new TF1("e1_default","[0]*exp(-x/[1])",0.0,1.3);
  e1_default->SetParameters(p_e1);
  e1_default->SetLineColor(kOrange);
  e1->SetParameters(p_e1);
  //e1->FixParameter(1,td);
  //h_elifetime_tpc1->Fit(e1,"L N","",0.0,bin_end);
  h_elifetime_tpc1->Fit(e1,"L N Q","",0.0,1.0);
  e1->SetLineColor(kRed);
  e1->SetLineStyle(2);
  e1->SetLineWidth(2);
  e1->Draw("same");
  myleg2->AddEntry(h_elifetime_tpc0,"TPC 0");
  myleg2->AddEntry(e0,Form("Fit TPC 0 (#tau_{0} = %.2f #pm %.2f ms)",e0->GetParameter(1),e0->GetParError(1)));
  myleg2->AddEntry(h_elifetime_tpc1,"TPC 1");
  myleg2->AddEntry(e1,Form("Fit TPC 1 (#tau_{1} = %.2f #pm %.2f ms)",e1->GetParameter(1),e1->GetParError(1)));
  myleg2->Draw("same");

  TLatex * l_expo = new TLatex();
  l_expo->SetTextFont(132);
  l_expo->DrawLatexNDC(0.2,0.9,Form("Fit = [0] exp(-#it{t}/[1])"));
  TLatex * l_fudge = new TLatex();
  l_fudge->SetTextFont(132);
  l_fudge->SetTextColor(kRed);
  if ( drift_vel != 160.0 ) {
    l_fudge->DrawLatexNDC(0.2,0.3,Form("Fudged data: drift vel. = %.2f cm/#mu s",drift_vel/1000.));
  }
  c1->SaveAs(Form("img/elifetime/%s/elifetime_tpc1.pdf",modifier.Data()));
  c1->SaveAs(Form("img/elifetime/%s/elifetime_tpc1.png",modifier.Data()));
  
}
