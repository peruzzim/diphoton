#!/bin/bash
cd /Users/peruzzi/smp13001_withjets/diphoton || exit 1
root -q -b -l template_studies_2d_variablebinning.C+O\(\"$1\",\"$2\",$3,\"$4\"\)



