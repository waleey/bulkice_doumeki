#!/bin/bash

#python3 merger/merge.py nakazato-LS220-BH-z0.004-s30.0.fits /home/waly/sntools_output/ibd.kin 10 mdom ibd 88 output 0

python3 merger/merge.py --progenitorModel nakazato-shen-z0.02-t_rev200ms-s50.0.fits --inputFormat SNEWPY-Nakazato_2013 --outfileS nakazato-shen-z0.02-t_rev200ms-s50.0.kin --distance 10 --omModel mdom --simType ibd --depthIndex 88 --outputFolderG output -t 0 -T 500 --runID 0 

#--transformation AdiabaticMSW_NMO
