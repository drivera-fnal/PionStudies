def set_style(hist, xtitle, ytitle, form=1 ):
  hist.SetTitle(";" + xtitle + ";" + ytitle )
  
  if form: 
    hist.GetXaxis().SetTitleSize(.05)
    hist.GetXaxis().SetTitleOffset(.85)
    hist.GetXaxis().SetLabelSize(.04)

    hist.GetYaxis().SetTitleSize(.05)
    hist.GetYaxis().SetTitleOffset(.85)
    hist.GetYaxis().SetLabelSize(.04)

def set_eff_style(hist, xtitle, ytitle, rebin=True):
  set_style(hist, xtitle, ytitle)

  hist.SetMarkerStyle(20)
  hist.SetLineColor(1)
  if( rebin ):
    hist.Rebin(2)
    hist.Scale(.5)
  hist.SetMaximum(1.)
  hist.SetMinimum(0.)


def markers(hist, color=1):
  hist.SetMarkerStyle(20)
  hist.SetMarkerSize(.5)
  hist.SetMarkerColor(color)
  hist.SetLineColor(color)

