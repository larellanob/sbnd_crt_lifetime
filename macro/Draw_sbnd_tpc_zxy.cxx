void Draw_sbnd_tpc()
{
  TPolyLine3D *sbnd_top     = new TPolyLine3D(5);
  TPolyLine3D *sbnd_bottom  = new TPolyLine3D(5);
  TPolyLine3D *sbnd_front   = new TPolyLine3D(5);
  TPolyLine3D *sbnd_back    = new TPolyLine3D(5);
  TPolyLine3D *sbnd_cathode = new TPolyLine3D(5);
  sbnd_bottom->SetPoint(0,  5.25,-202.05,-200);
  sbnd_bottom->SetPoint(1,  5.25, 202.05,-200);
  sbnd_bottom->SetPoint(2,504.15, 202.05,-200);
  sbnd_bottom->SetPoint(3,504.15,-202.05,-200);
  sbnd_bottom->SetPoint(4,  5.25,-202.05,-200);
  
  sbnd_top->SetPoint(0,  5.25,-202.05,200);
  sbnd_top->SetPoint(1,  5.25, 202.05,200);
  sbnd_top->SetPoint(2,504.15, 202.05,200);
  sbnd_top->SetPoint(3,504.15,-202.05,200);
  sbnd_top->SetPoint(4,  5.25,-202.05,200);

  sbnd_front->SetPoint(0,5.25,-202.05,-200);
  sbnd_front->SetPoint(1,5.25, 202.05,-200);
  sbnd_front->SetPoint(2,5.25, 202.05, 200);
  sbnd_front->SetPoint(3,5.25,-202.05, 200);
  sbnd_front->SetPoint(4,5.25,-202.05,-200);

  sbnd_back->SetPoint(0,504.15,-202.05,-200);
  sbnd_back->SetPoint(1,504.15, 202.05,-200);
  sbnd_back->SetPoint(2,504.15, 202.05, 200);
  sbnd_back->SetPoint(3,504.15,-202.05, 200);
  sbnd_back->SetPoint(4,504.15,-202.05,-200);

  sbnd_cathode->SetPoint(0,  5.25, 0,-200);
  sbnd_cathode->SetPoint(1,  5.25, 0, 200);
  sbnd_cathode->SetPoint(2,504.15, 0, 200);
  sbnd_cathode->SetPoint(3,504.15, 0,-200);
  sbnd_cathode->SetPoint(4,  5.25, 0,-200);

  std::vector<TPolyLine3D*> v_sbnd_tpc { sbnd_bottom,
					 sbnd_top,
					 sbnd_front,
					 sbnd_back,
					 sbnd_cathode};
  

  for ( int i = 0; i < v_sbnd_tpc.size(); i++ ) {
    v_sbnd_tpc[i]->Draw("same");
    if ( i != 4 ) { // keep cathode line thinner
      v_sbnd_tpc[i]->SetLineWidth(2);
    }
  }
  

}

void Draw_sbnd_crt()
{
  TPolyLine3D * crt_top    = new TPolyLine3D(5);
  TPolyLine3D * crt_bottom = new TPolyLine3D(5);
  TPolyLine3D * crt_front  = new TPolyLine3D(5);
  TPolyLine3D * crt_back   = new TPolyLine3D(5);
  TPolyLine3D * crt_xp     = new TPolyLine3D(5);
  TPolyLine3D * crt_xm     = new TPolyLine3D(5);
  // z x y
  crt_front->SetPoint(0,-179.05,-360,-360);
  crt_front->SetPoint(1,-179.05, 360,-360);
  crt_front->SetPoint(2,-179.05, 360, 360);
  crt_front->SetPoint(3,-179.05,-360, 360);
  crt_front->SetPoint(4,-179.05,-360,-360);

  crt_back->SetPoint(0,770.95,-360,-360);
  crt_back->SetPoint(1,770.95, 360,-360);
  crt_back->SetPoint(2,770.95, 360, 360);
  crt_back->SetPoint(3,770.95,-360, 360);
  crt_back->SetPoint(4,770.95,-360,-360);

  std::vector<TPolyLine3D * > v_crt{
    crt_front,
    crt_back
  };
  for ( int i = 0; i < v_crt.size(); i++ ) {
    v_crt[i]->Draw("same");
    //v_crt[i]->SetLineWidth(2);
  }
}


void Draw_invisible_box()
{
  // draw invisible axes to increase the size of the empty canvas
  // relative to the sbnd box
  TPolyLine3D *inv_x = new TPolyLine3D(2);
  TPolyLine3D *inv_y = new TPolyLine3D(2);
  TPolyLine3D *inv_z = new TPolyLine3D(2);

  double xhi =  0.0;
  double xlo =  0.0;
  double yhi =  500.0;
  double ylo =  0.0;
  double zlo =  0.0;
  double zhi =  0.0;
  
  
  //std::cout << "Drawing invisible" << std::endl;
  // z x y
  inv_x->SetPoint(0, 0.0, xlo, 0.0);
  inv_x->SetPoint(1, 0.0, xhi, 0.0);
  inv_y->SetPoint(0, 0.0, 0.0, ylo);
  inv_y->SetPoint(1, 0.0, 0.0, yhi);
  inv_z->SetPoint(0, zlo, 0.0, 0.0);
  inv_z->SetPoint(1, zhi, 0.0, 0.0);

  inv_x->SetLineColor(kWhite);
  inv_y->SetLineColor(kWhite);
  inv_z->SetLineColor(kWhite);
  inv_x->Draw("same");
  inv_y->Draw("same");
  inv_z->Draw("same");

}
