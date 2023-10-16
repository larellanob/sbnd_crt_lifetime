#include "macro/Draw_sbnd_tpc_zxy.cxx"

void draw_aux_lines(double x, double y, double z, TString tpccrt = "crt")
{
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

void Plot_matched_tracks_zxy()
{
  gROOT->GetStyle("uboone_sty_colz");
  uboone_sty_colz->cd();
  uboone_sty_colz->SetPadTopMargin(0.1);

  int max_ntracks = 150;

  // open file
  //TString filename {"CORRECTED_muon_hitdumper_NS_231010_1240.root"};
  TString filename {"CORRECTED_muon_hitdumper_NS.root"};
  TFile *f = new TFile(filename);

  // read file
  TTreeReader reader("hitdumper/hitdumpertree_corr",f);
  // evt run subrun
  //TTreeReaderValue<int> run (reader,"run");
  //TTreeReaderValue<int> subrun (reader,"subrun");
  TTreeReaderValue<int> evt (reader,"event"); 
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
  TTreeReaderValue<int> ev(reader,"event");

  // canvases
  TCanvas *c1 = new TCanvas("c1","mycanvas",1000,1000);

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
  TLatex latex_anglescrt;
  latex_anglescrt.SetTextFont(132);
  TLatex latex_anglestpc;
  latex_anglestpc.SetTextFont(132);
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
    if ( events > 20 ) break;
    
    c1->Clear();

    // draw background
    Draw_invisible_box();
    Draw_sbnd_tpc();
    Draw_sbnd_crt();
    // see edges of tarrows
    coord_syst_x.DrawLatex(arr_cx-arr_l*sin(axis_phrad)-0.07,
		     arr_cy+arr_l*sin(axis_thrad)*cos(axis_phrad)-0.05,
		     "x");
    coord_syst_y.DrawLatex(arr_cx-0.02,
			   arr_cy+arr_l*cos(axis_thrad)+0.03,
			   "y");
    coord_syst_z.DrawLatex(arr_cx+arr_l*cos(axis_phrad)+0.04,
		    arr_cy+arr_l*sin(axis_thrad)*sin(axis_phrad)+0.01,
		    "z");
    arrx->Draw();
    arry->Draw();
    arrz->Draw();

    for ( int i = 0; i < muontrk_x1.GetSize(); i++ ) {
      latex_anglescrt.DrawLatex(0.1,0.9,
				Form("#theta_{XZ}^{CRT}: %.2f, #theta_{YZ}^{CRT}: %.2f",
				     crt_xz[i],crt_yz[i]));
      latex_anglestpc.DrawLatex(0.1,0.7,
				Form("#theta_{XZ}^{TPC}: %.2f, #theta_{YZ}^{TPC}: %.2f",
				     tpc_xz[i],tpc_yz[i]));
      TPolyLine3D *pl_tpc = new TPolyLine3D(2);
      pl_tpc->SetLineStyle(9);
      pl_tpc->SetLineWidth(3);
      pl_tpc->SetLineColor(kRed);
      pl_tpc->SetPoint(0,muontrk_z1[i],muontrk_x1[i],muontrk_y1[i]);
      pl_tpc->SetPoint(1,muontrk_z2[i],muontrk_x2[i],muontrk_y2[i]);
      std::cout << muontrk_z2[i]<< " " <<muontrk_x2[i]<< " " <<muontrk_y2[i] << std::endl;
      TPolyLine3D *pl_crt = new TPolyLine3D(2);
      pl_crt->SetLineWidth(3);
      pl_crt->SetLineColor(kBlue);
      pl_crt->SetPoint(0,ct_z1[i],ct_x1[i],ct_y1[i]);
      pl_crt->SetPoint(1,ct_z2[i],ct_x2[i],ct_y2[i]);
      // uncorrected muontrk
      TPolyLine3D *pl_tpc_uncorr = new TPolyLine3D(2);
      pl_tpc_uncorr->SetLineWidth(2);
      pl_tpc_uncorr->SetLineColor(kOrange);
      pl_tpc_uncorr->SetPoint(0,muontrk_uncorrected_z1[i],muontrk_uncorrected_x1[i],muontrk_uncorrected_y1[i]);
      pl_tpc_uncorr->SetPoint(1,muontrk_uncorrected_z2[i],muontrk_uncorrected_x2[i],muontrk_uncorrected_y2[i]);
      
      draw_aux_lines(ct_x1[i],ct_y1[i],ct_z1[i]);
      draw_aux_lines(ct_x2[i],ct_y2[i],ct_z2[i]);
      draw_aux_lines(muontrk_x1[i],muontrk_y1[i],muontrk_z1[i],"tpc");
      draw_aux_lines(muontrk_x2[i],muontrk_y2[i],muontrk_z2[i],"tpc");

      // legend
      TLegend *myleg = new TLegend(0.0,0.8,0.37,1.0);
      myleg->AddEntry(pl_tpc,"TPC track (matched)");
      myleg->AddEntry(pl_crt,"CRT track");
      myleg->AddEntry(pl_tpc_uncorr,"TPC track (uncorrected)");      

      myleg->Draw("same");
      pl_crt->Draw("same");
      pl_tpc_uncorr->Draw("same");
      pl_tpc->Draw("same");
      c1->SaveAs(Form("img/matched/poly_ev%i.pdf",*ev));
      c1->SaveAs(Form("img/matched/poly_ev%i.png",*ev));
    }
  }
}
