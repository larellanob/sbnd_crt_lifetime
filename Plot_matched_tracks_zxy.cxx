#include "macro/Draw_sbnd_tpc_zxy.cxx"

void draw_aux_lines(double x, double y, double z, TString tpccrt = "crt")
{
  /**
     this function draws dotted lines on the walls of tpc/crt touched
     by the crossing track.

     currently only front-back faces of crt/tpc are implemented
   */
  gROOT->SetStyle("uboone_sty");

  // tpolyline3ds will hold tpc/crt tracks to be drawn
  TPolyLine3D * line_x = new TPolyLine3D(2);
  TPolyLine3D * line_y = new TPolyLine3D(2);
  TPolyLine3D * line_z = new TPolyLine3D(2);
  line_x->SetLineColor(kBlue);
  line_y->SetLineColor(kBlue);
  line_z->SetLineColor(kBlue);
  line_x->SetLineStyle(2);
  line_y->SetLineStyle(2);
  line_z->SetLineStyle(2);
  // SetPoint(z,x,y)
  // hardcoding with up-downstream crts for now
  // up-downstream
  if ( std::abs(z+179.05) < 0.1  || std::abs(z-770.95) < 0.1 ) {
    line_x->SetPoint(0,z,-360,y);
    line_x->SetPoint(1,z,360,y);
    line_y->SetPoint(0,z,x,360);
    line_y->SetPoint(1,z,x,-360);
  }

  if ( tpccrt == "tpc") {
    line_x->SetLineColor(kRed);
    line_y->SetLineColor(kRed);
    line_z->SetLineColor(kRed);
    if ( std::abs(z-5.25) < 0.1  || std::abs(z-504.15) < 0.1 ) {
      line_x->SetPoint(0,z,-200,y);
      line_x->SetPoint(1,z,200,y);
      line_y->SetPoint(0,z,x,200);
      line_y->SetPoint(1,z,x,-200);
    }
  }

  line_x->Draw("same");
  line_y->Draw("same");
  line_z->Draw("same");
}

void Draw_coordinate_system(double theta = 30, double phi = 60)
{
  /**
     Draws the coordinate system on the canvas. call this function
     with the same theta and phi angles set in the gpad
   */
  // convert to rad
  theta = theta*TMath::DegToRad();
  phi = phi*TMath::DegToRad();
  TLatex * coord_syst_x = new TLatex();
  TLatex * coord_syst_y = new TLatex();
  TLatex * coord_syst_z = new TLatex();
  coord_syst_x->SetTextFont(12);
  coord_syst_y->SetTextFont(12);
  coord_syst_z->SetTextFont(12);

  double arr_cx =  0.7; // starting pos in x
  double arr_cy = -0.8; // starting pos in y
  double arr_l = 0.2;  // arrow "stick" length
  // it took me hours to get these angles right, enjoy!
  TArrow *arrx = new TArrow(arr_cx,
			    arr_cy,
			    arr_cx-arr_l*sin(phi),
			    arr_cy+arr_l*sin(theta)*cos(phi),
			    0.01, "|>");
  // y axis "behaves line" z axis
  TArrow *arry = new TArrow(arr_cx,
			    arr_cy,
			    arr_cx,
			    arr_cy+arr_l*cos(theta),
			    0.01,
			    "|>");
  TArrow *arrz = new TArrow(arr_cx,
			    arr_cy,
			    arr_cx+arr_l*cos(phi),
			    arr_cy+arr_l*sin(theta)*sin(phi),
			    0.01,
			    "|>");
  arrx->SetLineWidth(2);
  arry->SetLineWidth(2);
  arrz->SetLineWidth(2);

  // draw axis labels near the pointy end of tarrows
  coord_syst_x->DrawLatex(arr_cx-arr_l*sin(phi)-0.07,
			  arr_cy+arr_l*sin(theta)*cos(phi)-0.05,
			  "x");
  coord_syst_y->DrawLatex(arr_cx-0.02,
			  arr_cy+arr_l*cos(theta)+0.03,
			  "y");
  coord_syst_z->DrawLatex(arr_cx+arr_l*cos(phi)+0.04,
			  arr_cy+arr_l*sin(theta)*sin(phi)+0.01,
			  "z");
  arrx->Draw();
  arry->Draw();
  arrz->Draw();
}
  

void Plot_matched_tracks_zxy()
{
  gROOT->GetStyle("uboone_sty_colz");
  uboone_sty_colz->cd();
  uboone_sty_colz->SetPadTopMargin(0.1);
  uboone_sty_colz->SetPadRightMargin(0.01);
  uboone_sty_colz->SetPadLeftMargin(0.01);

  // open file
  TString filename {"CORRECTED_muon_hitdumper_NS.root"};
  TFile *f = new TFile(filename);

  // read file
  TTreeReader reader("hitdumper/hitdumpertree_corr",f);
  // muon tracks
  TTreeReaderArray<double> muontrk_x1(reader,"muontrk_x1");
  TTreeReaderArray<double> muontrk_y1(reader,"muontrk_y1");
  TTreeReaderArray<double> muontrk_z1(reader,"muontrk_z1");
  TTreeReaderArray<double> muontrk_x2(reader,"muontrk_x2");
  TTreeReaderArray<double> muontrk_y2(reader,"muontrk_y2");
  TTreeReaderArray<double> muontrk_z2(reader,"muontrk_z2");
  // uncorrected muon tracks
  TTreeReaderArray<double> muontrk_uncorrected_x1(reader,"muontrk_uncorrected_x1");
  TTreeReaderArray<double> muontrk_uncorrected_y1(reader,"muontrk_uncorrected_y1");
  TTreeReaderArray<double> muontrk_uncorrected_z1(reader,"muontrk_uncorrected_z1");
  TTreeReaderArray<double> muontrk_uncorrected_x2(reader,"muontrk_uncorrected_x2");
  TTreeReaderArray<double> muontrk_uncorrected_y2(reader,"muontrk_uncorrected_y2");
  TTreeReaderArray<double> muontrk_uncorrected_z2(reader,"muontrk_uncorrected_z2");
  // crt tracks
  TTreeReaderArray<double> ct_x1(reader,"ct_x1");
  TTreeReaderArray<double> ct_y1(reader,"ct_y1");
  TTreeReaderArray<double> ct_z1(reader,"ct_z1");
  TTreeReaderArray<double> ct_x2(reader,"ct_x2");
  TTreeReaderArray<double> ct_y2(reader,"ct_y2");
  TTreeReaderArray<double> ct_z2(reader,"ct_z2");
  // angles
  TTreeReaderArray<double> crt_xz(reader,"crt_thetaxz");
  TTreeReaderArray<double> crt_yz(reader,"crt_thetayz");
  TTreeReaderArray<double> tpc_xz(reader,"tpc_thetaxz");
  TTreeReaderArray<double> tpc_yz(reader,"tpc_thetayz");
  // misc
  // ev = entry in comm trees
  TTreeReaderValue<int> ev(reader,"event");
  TTreeReaderValue<int> n_matched(reader,"n_matched");
  // evt = event number assigned in simulation
  TTreeReaderValue<int> evt(reader,"evt");
  TTreeReaderValue<int> run(reader,"run");
  TTreeReaderValue<int> subrun(reader,"subrun");

  // canvases
  TCanvas *c1 = new TCanvas("c1","mycanvas",1000,1000);
  TCanvas *cxz = new TCanvas("cxz","mycanvas",1000,1000);
  TCanvas *cxz_ang = new TCanvas("cxz_ang","mycanvas",1800,1000);
  TCanvas *cyz = new TCanvas("cyz","mycanvas",1000,1000);
  TCanvas *cangles = new TCanvas("cangles","mycanvas",500,1000);

  // some manually placed labels
  TLatex coord_syst_x;
  TLatex coord_syst_y;
  TLatex coord_syst_z;
  coord_syst_x.SetTextFont(12);
  coord_syst_y.SetTextFont(12);
  coord_syst_z.SetTextFont(12);
  // orientation of the drawing
  double axis_theta = 30; // default 30
  double axis_phi   = 60; // default 60
  double axis_thrad = axis_theta*TMath::DegToRad();
  double axis_phrad = axis_phi*TMath::DegToRad();
  
  gPad->SetTheta(axis_theta);
  gPad->SetPhi(axis_phi);
  TLatex latex_anglestpc;
  latex_anglestpc.SetTextFont(132);
  TLatex latex_anglescrt;
  latex_anglescrt.SetTextFont(132);
  // z x y
  double arr_cx =  0.7; // starting pos in x
  double arr_cy = -0.8; // starting pos in y
  double arr_l = 0.2;  // arrow "stick" length
  // it took me hours to get these angles right, enjoy!
  TArrow *arrx = new TArrow(arr_cx,
			    arr_cy,
			    arr_cx-arr_l*sin(axis_phrad),
			    arr_cy+arr_l*sin(axis_thrad)*cos(axis_phrad),
			    0.01, "|>");
  // y axis "behaves line" z axis
  TArrow *arry = new TArrow(arr_cx,
			    arr_cy,
			    arr_cx,
			    arr_cy+arr_l*cos(axis_thrad),
			    0.01,
			    "|>");
  TArrow *arrz = new TArrow(arr_cx,
			    arr_cy,
			    arr_cx+arr_l*cos(axis_phrad),
			    arr_cy+arr_l*sin(axis_thrad)*sin(axis_phrad),
			    0.01,
			    "|>");
  arrx->SetLineWidth(2);
  arry->SetLineWidth(2);
  arrz->SetLineWidth(2);

  // loop through file
  int events = -1;
  while ( reader.Next() ) {
    events++;
    //if ( events > 26  && *n_matched < 2) continue;
    //if ( *ev != 989 ) continue;
    if ( events > 26  && *n_matched < 2) continue;
    //if ( events > 2) continue;
    
    c1->Clear();
    cxz->Clear();
    cxz_ang->Clear();
    cyz->Clear();
    cangles->Clear();
    cxz_ang->Clear();
    cxz_ang->cd();
    double separation = 0.6;
    TPad *p1 = new TPad("p1","p1",0.0,0.0,separation,1.0);
    TPad *p2 = new TPad("p2","p2",separation,0.0,1.0,1.0);
    p2->SetRightMargin(0.0);
    p1->Draw();
    p2->Draw();
    std::vector<TCanvas *> v_canvases {c1,cxz,cxz_ang,cyz,cangles};
    

    // draw background
    for ( auto c: v_canvases ) {
      c->cd();
      if ( c == cangles ) {
	continue;
      }
      if ( c == c1 ) {
	gPad->SetTheta(axis_theta);
	gPad->SetPhi(axis_phi);
	Draw_coordinate_system(axis_theta,axis_phi);
      } else if ( c == cyz ) {
	gPad->SetTheta(0);
	gPad->SetPhi(0);
	Draw_coordinate_system(0,0);
      } else if ( c == cxz ) {
	gPad->SetTheta(90);
	gPad->SetPhi(0);
	Draw_coordinate_system(90,0);
      } else if ( c == cxz_ang ) {
	//cxz_ang->cd(0);
	p1->cd();
	gPad->SetTheta(90);
	gPad->SetPhi(0);
	Draw_coordinate_system(90,0);
      }
      Draw_invisible_box();
      Draw_sbnd_tpc();
      Draw_sbnd_crt();
    }
    // loop over muon tracks
    int track_line_width = 5;
    for ( int i = 0; i < muontrk_x1.GetSize(); i++ ) {

      TPolyLine3D *pl_tpc = new TPolyLine3D(2);
      pl_tpc->SetLineStyle(9);
      pl_tpc->SetLineWidth(track_line_width);
      pl_tpc->SetLineColor(kRed+i);
      pl_tpc->SetPoint(0,muontrk_z1[i],muontrk_x1[i],muontrk_y1[i]);
      pl_tpc->SetPoint(1,muontrk_z2[i],muontrk_x2[i],muontrk_y2[i]);
      std::cout << muontrk_z2[i]<< " " <<muontrk_x2[i]<< " " <<muontrk_y2[i] << std::endl;
      TPolyLine3D *pl_crt = new TPolyLine3D(2);
      pl_crt->SetLineWidth(track_line_width);
      pl_crt->SetLineColor(kBlue);
      pl_crt->SetPoint(0,ct_z1[i],ct_x1[i],ct_y1[i]);
      pl_crt->SetPoint(1,ct_z2[i],ct_x2[i],ct_y2[i]);
      // uncorrected muontrk
      TPolyLine3D *pl_tpc_uncorr = new TPolyLine3D(2);
      pl_tpc_uncorr->SetLineWidth(track_line_width);
      pl_tpc_uncorr->SetLineColor(kOrange+i);
      pl_tpc_uncorr->SetPoint(0,muontrk_uncorrected_z1[i],muontrk_uncorrected_x1[i],muontrk_uncorrected_y1[i]);
      pl_tpc_uncorr->SetPoint(1,muontrk_uncorrected_z2[i],muontrk_uncorrected_x2[i],muontrk_uncorrected_y2[i]);

      // legend
      TLegend *myleg = new TLegend(0.0,0.8,0.37,1.0);
      myleg->AddEntry(pl_tpc,"TPC track (matched)");
      myleg->AddEntry(pl_crt,"CRT track");
      myleg->AddEntry(pl_tpc_uncorr,"TPC track (uncorrected)");      

      for ( auto c: v_canvases ) {
	c->cd();
	TString filename_modifier;

	if ( c == c1 ) {
	  TLatex l_rse;
	  l_rse.SetTextFont(132);
	  l_rse.DrawLatexNDC(0.01,0.01,Form("Run: %i, Subrun: %i, Event: %i",*run,*subrun,*evt));
	  filename_modifier = "";
	  myleg->Draw("same");
	} else if ( c == cyz ) {
	  filename_modifier = "_yz";
	  
	} else if ( c == cxz ) {
	  filename_modifier = "_xz";
	}
	// need to draw previous angles until exhausted
	int j = 0;
	if ( c == cangles || c == cxz_ang ) {
	  filename_modifier = "_angles";
	  if ( c == cxz_ang ) {
	    filename_modifier = "_xzang";
	    myleg->Draw("same");
	    p2->cd();
	  }
	  while ( j <= i ) {
	    latex_anglescrt.DrawLatex(0.1,0.90-j*0.3,
				      Form("#scale[1.8]{#theta_{XZ}^{CRT} %i: %.2f, #theta_{YZ}^{CRT}: %.2f}",
					   j+1,crt_xz[j],crt_yz[j]));
	    latex_anglestpc.DrawLatex(0.1,0.80-j*0.3,
				      Form("#scale[1.8]{#theta_{XZ}^{TPC} %i: %.2f, #theta_{YZ}^{TPC}: %.2f}",
					   j+1,tpc_xz[j],tpc_yz[j]));
	    j++;
	  }
	  if ( c == cxz_ang ) p1->cd();
	}
	if ( c != cangles ) {
	  // draw stuff
	  // auxiliary lines for crt, front and back
	  draw_aux_lines(ct_x1[i],ct_y1[i],ct_z1[i]);
	  draw_aux_lines(ct_x2[i],ct_y2[i],ct_z2[i]);
	  // auxiliary lines for tpc, front and back
	  draw_aux_lines(muontrk_x1[i],muontrk_y1[i],muontrk_z1[i],"tpc");
	  draw_aux_lines(muontrk_x2[i],muontrk_y2[i],muontrk_z2[i],"tpc");
	  // legend
	  // polylines (tracks) for crt, tpc uncorrected and tpc
	  bool redraw_prev_tpc = false;
	  if ( i > 0 ) {
	    // if more than one track per event, check that CRT track is different
	    // if it's the same, have to redraw the previous muon track so it doesn't get blocked
	    if ( ct_x1[i] == ct_x1[0] && ct_y1[i] == ct_y1[0] && ct_z1[i] == ct_z1[0] &&
		 ct_x2[i] == ct_x2[0] && ct_y2[i] == ct_y2[0] && ct_z2[i] == ct_z2[0] ) {
	      redraw_prev_tpc = true;
	    }
	  }
	  pl_crt->Draw("same");
	  if ( redraw_prev_tpc) {
	    TPolyLine3D *pl_tpc_prev = new TPolyLine3D(2);
	    pl_tpc_prev->SetLineStyle(9);
	    pl_tpc_prev->SetLineWidth(track_line_width);
	    pl_tpc_prev->SetLineColor(kRed+i-1);
	    pl_tpc_prev->SetPoint(0,muontrk_z1[i-1],muontrk_x1[i-1],muontrk_y1[i-1]);
	    pl_tpc_prev->SetPoint(1,muontrk_z2[i-1],muontrk_x2[i-1],muontrk_y2[i-1]);
	    pl_tpc_prev->Draw("same");
	  }
	  pl_tpc_uncorr->Draw("same");
	  pl_tpc->Draw("same");
	}

	// save in the apropriate canvas (viewing angle)
	std::cout << "saving " << filename_modifier << std::endl;
	c->SaveAs(Form("img/matched/poly_ev%i%s_match%i.pdf",*ev,filename_modifier.Data(),i));
	c->SaveAs(Form("img/matched/poly_ev%i%s_match%i.png",*ev,filename_modifier.Data(),i));

      }
    } // end of loop over muon tracks
  } // end of loop though events in the file
}
