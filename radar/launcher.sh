#!/bin/bash

#Export conda path (on idroclima the standard python is 2.7)
export PATH=/home/silvestro/fp_virtualenv_python3_hyde/bin:$PATH
source activate wradlib

python3 load_vols_in_dir.py
