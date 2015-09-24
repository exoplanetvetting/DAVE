How to get K2 data from MAST in DAVE

import mastio
path = "/media/external_drive/K2"  #For example
ar = mastio.K2Archive(path)

k2id = 123456789
campaign = 3
fits = ar.getLongCadence(k2id, campaign)
tpf = ar.getLongTpf(k2id, campaign)

Optional Arguments to pyfits.getdata are also optional arguments
to ar's get methods. For example
tpf, hdr = ar.getLongTpf(k2id, campaign, header=True)
