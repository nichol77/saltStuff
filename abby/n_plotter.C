{
  
  gROOT->Reset();
  
 
  float christian_depth[8]={10,20,30,40,50,60,70,90};
  float long_depth[3]={20,50,75};
  float baby_depth[3]={10,50,90};
  float christian_n[8]={2.49422469098865,2.48195475233835,2.48195475233835,
			2.46968481368804,2.45741487503773,2.44514493638743,
			2.43287499773712,2.38379524313589};
  float long_n[3]={2.55670935356534,2.46468481368804,2.43400996706227};
  float baby_n[3]={2.50540567258375,2.46859585663283,2.42565107135675};
  float christian_up[8]={2.52231386579657,2.51354202528368,2.51677259522001,
			 2.50733011901294,2.49764925675412,2.4877723971505,
			 2.47773022826549,2.43226980009899};
  float christian_down[8]={2.46632578618417,2.45067658644869,2.44758405860269,
			   2.43263019550368,2.41791848216439,2.40340460254789,
			   2.389056772044,2.33664159731293};
  float long_up[3]={2.58892816372483,2.50501334554841,2.48004000181354};
  float long_down[3]={2.52481215328435,2.42509691905278,2.38909581173286};
  float baby_up[3]={2.53356677348256,2.5089750628488,2.4748549774253};
  float baby_down[3]={2.47743606354905,2.42895871219895,2.37779469013755};
  
  float christian_dist[8]={554.877628937229,554.469465671923,554.061302406616,
			   553.244975876004,552.428649345392,551.61232281478,
			   550.795996284168,549.979669753555};
  float long_dist[3]={553.320486080086,550.871506488249,549.238853427025};
  float baby_dist[3]={555.621506488249,554.397016692331,551.539873835188};
  float c_dist_true[8]={545.039184243403,547.330764034421,546.927855802892,
			548.835294035794,550.761780049626,552.707600594387,
			554.673048204848,565.254165506109};
  float l_dist_true[3]={530.226554303395,547.589364530826,552.84703395045};
  float b_dist_true[3]={543.334241553135,550.220761023597,557.076286384586};

  //for putting systematic error bars on graph
  float errs_low[4]={0.041,0.015,0.018,0.047}; //wander_max, distance,timing, quadrature 
  float errs_high[4]={0.042,0.015,0.018,0.048};  
  float fake_depth_err[4]={20,40,60,80};
  float fake_point_err[4]={2.7,2.7,2.7,2.7};

  TGraphAsymmErrors *glong, *gchrist, *gbaby, *gerror, *gdistchrist, *gdistlong, *gdistbaby;
  TGraphAsymmErrors *gdtruec, *gdtruel, *gdtrueb;
  TCanvas *c1=new TCanvas("c1","c1",1000,800);
  //c1->Divide(3,1);
  TLegend *leg1;
  
  for (int i=0;i<8;i++){
    christian_up[i]=christian_up[i]-christian_n[i];
    christian_down[i]=christian_n[i]-christian_down[i];
    if (i<3){
      long_up[i]=long_up[i]-long_n[i];
      long_down[i]=long_n[i]-long_down[i];
      baby_up[i]=baby_up[i]-baby_n[i];
      baby_down[i]=baby_n[i]-baby_down[i];
    }
  }
  

  gchrist=new TGraphAsymmErrors(8,christian_depth,christian_n, 0,0,
				0,0);
  //christian_down, christian_up);
  glong=new TGraphAsymmErrors(3,long_depth,long_n, 0,0,
			      0,0);
  //long_down, long_up);
  gbaby=new TGraphAsymmErrors(3,baby_depth,baby_n, 0,0,
			      0,0);
  //baby_down, baby_up);
  gerror=new TGraphAsymmErrors(4,fake_depth_err,fake_point_err,0,0,
 		       errs_low, errs_high);
  gdistchrist=new TGraphAsymmErrors(8,christian_depth,christian_dist,0,0,
				   0,0);
  gdistlong=new TGraphAsymmErrors(3,long_depth,long_dist,0,0,
				   0,0);
  gdistbaby=new TGraphAsymmErrors(3,baby_depth,baby_dist,0,0,
				   0,0);
  gdtruec=new TGraphAsymmErrors(8,christian_depth,c_dist_true,0,0,
				   0,0);
  gdtruel=new TGraphAsymmErrors(3,long_depth,l_dist_true,0,0,
				   0,0);
  gdtrueb=new TGraphAsymmErrors(3,baby_depth,b_dist_true,0,0,
				   0,0);
  
  TH2F *h1=new TH2F("Index of Refraction versus Depth (ft)",
		    "Index of Refraction versus Depth (ft)"
		    ,100,0,100,100,2.30,2.76);
  h1->SetXTitle("Depth (ft)");
  h1->SetYTitle("Index of Refraction");
  TH2F *h2=new TH2F("Distance between 1 & 3 versus Depth (ft)",
		    "Distance between 1 & 3 versus Depth (ft)"
		    ,100,0,100,100,525,570);
  h2->SetXTitle("Depth (ft)");
  h2->SetYTitle("Distance between 1 & 3");
  

  //add a legend
  leg1= new TLegend(0.65,0.50,.85,0.65);
  leg1->AddEntry(glong,"LF Antennas","l");
  leg1->AddEntry(gchrist,"Midband Antennas","l");
  leg1->AddEntry(gbaby,"HF Antennas","l");
  
  c1->cd(0);
  gStyle->SetOptStat(kFALSE);
  h1->Draw();
  gchrist->SetMarkerSize(1.1);
  gchrist->SetMarkerStyle(21);
  gchrist->SetLineColor(1);//green
  gchrist->SetMarkerColor(1);//
  gchrist->Draw("P");
  glong->SetMarkerSize(1.1);
  glong->SetMarkerStyle(22);
  glong->SetLineColor(1);//red
  glong->SetMarkerColor(1);//red
  glong->Draw("P");
  gbaby->SetMarkerSize(1.1);
  gbaby->SetMarkerStyle(20);
  gbaby->SetLineColor(1);//blue
  gbaby->SetMarkerColor(1);//blue
  gbaby->Draw("P");
  gerror->SetMarkerSize(0);
  gerror->SetLineColor(1);
  gerror->Draw("P");
  leg1->Draw();


  /*c1->cd(2);
  h2->Draw();
  gdistchrist->SetMarkerSize(1);
  gdistchrist->SetMarkerStyle(21);
  gdistchrist->SetLineColor(3);//green
  gdistchrist->SetMarkerColor(3);//
  gdistchrist->Draw("P");
  gdistlong->SetMarkerSize(1);
  gdistlong->SetMarkerStyle(21);
  gdistlong->SetLineColor(2);//red
  gdistlong->SetMarkerColor(2);//red
  gdistlong->Draw("P");
  gdistbaby->SetMarkerSize(1);
  gdistbaby->SetMarkerStyle(21);
  gdistbaby->SetLineColor(4);//blue
  gdistbaby->SetMarkerColor(4);//blue
  gdistbaby->Draw("P");
  leg1->Draw();
  c1->cd(3);
  h2->Draw();
  gdtruec->SetMarkerSize(1);
  gdtruec->SetMarkerStyle(21);
  gdtruec->SetLineColor(3);//green
  gdtruec->SetMarkerColor(3);//
  gdtruec->Draw("P");
  gdtruel->SetMarkerSize(1);
  gdtruel->SetMarkerStyle(21);
  gdtruel->SetLineColor(2);//red
  gdtruel->SetMarkerColor(2);//red
  gdtruel->Draw("P");
  gdtrueb->SetMarkerSize(1);
  gdtrueb->SetMarkerStyle(21);
  gdtrueb->SetLineColor(4);//blue
  gdtrueb->SetMarkerColor(4);//blue
  gdtrueb->Draw("P");
  leg1->Draw();
  */
}

