void overlay_histos()
{

  //TCut *cuts = new TCut("")

  SNT->Draw("theta_rq*180/TMath::Pi()>>(19,0,190)", "Weight*Normfac/1e6*(abs(Em)<=0.040&&abs(h_delta)<=10.)&&e_delta>-10&&e_delta<22.");
  SNT->Draw("Pm>>(15,0,1.2)", "Weight*Normfac/1e6*(abs(Em)<=0.040&&abs(h_delta)<=10.&&e_delta>-10&&e_delta<22.&&theta_rq*180/TMath::Pi()<40.)");
  

    
}
