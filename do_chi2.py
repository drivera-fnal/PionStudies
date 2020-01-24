def do_chi2(e, i, profile):
  dEdXs = [j for j in e.reco_daughter_allTrack_calibrated_dEdX_SCE[i]]
  #dEdXs = [j for j in e.reco_daughter_allTrack_dEdX[i]]
  resRanges = [j for j in e.reco_daughter_allTrack_resRange[i]]


  if( len(dEdXs) < 1 or len(resRanges) < 1 ):
    return [9999., -1.]

  npt = 0
  pid_chi2 = 0.

  for j in range(1, len(dEdXs)-1): 
    if dEdXs[j] > 1000.: continue

    theBin = profile.FindBin(resRanges[j])
    if not ( theBin >= 1 and theBin <= profile.GetNbinsX() ):
      continue


    template_dEdX = profile.GetBinContent(theBin)
    if template_dEdX < 1.e-6:
      template_dEdX = ( profile.GetBinContent( theBin - 1 ) + profile.GetBinContent( theBin + 1 ) ) / 2.

    template_dEdX_err = profile.GetBinError( theBin )
    if( template_dEdX_err < 1.e-6 ):
      template_dEdX_err = ( profile.GetBinError( theBin - 1 ) + profile.GetBinError( theBin + 1 ) ) / 2.            


    dedx_res = 0.04231 + 0.0001783 * dEdXs[j] * dEdXs[j]
    dedx_res *= dEdXs[j]

    pid_chi2 += (  (dEdXs[j] - template_dEdX)**2  / ( template_dEdX_err**2  + dedx_res**2 ) )

    npt += 1
 
  

  if npt == 0:
    return [9999., -1.]

  return [pid_chi2, npt]
