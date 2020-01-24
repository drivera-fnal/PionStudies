def add_stack(s):
  h = s.GetHists().At(0).Clone("")

  for i in range(1, s.GetNhists()):
    h.Add( s.GetHists().At(i) )
  
  return h.Integral()
