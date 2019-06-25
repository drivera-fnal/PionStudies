def set_style(hist, xtitle, ytitle):
  hist.SetTitle(";" + xtitle + ";" + ytitle )

  hist.GetXaxis().SetTitleSize(.05)
  hist.GetXaxis().SetTitleOffset(.85)
  hist.GetXaxis().SetLabelSize(.04)

  hist.GetYaxis().SetTitleSize(.05)
  hist.GetYaxis().SetTitleOffset(.85)
  hist.GetYaxis().SetLabelSize(.04)

def set_eff_style(hist, xtitle, ytitle):
  set_style(hist, xtitle, ytitle)

  hist.SetMarkerStyle(20)
  hist.SetLineColor(1)
  hist.Rebin(2)
  hist.Scale(.5)
  hist.SetMaximum(1.)
  hist.SetMinimum(0.)
