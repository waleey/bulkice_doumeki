#!/bin/bash

#python3 merge.py nakazato-LS220-BH-z0.004-s30.0.fits /home/waly/sntools_output/radioactivity.kin 10 mdom radioactivity run_pos.mac

#for i in {1..60}
#do
python3 merger/merge.py nakazato-LS220-BH-z0.004-s30.0.fits /home/waly/sntools_output/radioactivity.kin 10 degg radioactivity 88 output 0
#done
