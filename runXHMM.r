#https://software.broadinstitute.org/gatk/documentation/tooldocs/current/
# XHMM workflow and install
# http://atgu.mgh.harvard.edu/xhmm/download.shtml


wget https://bitbucket.org/statgen/xhmm/get/cc14e528d909.zip
unzip cc14e528d909.zip
cd statgen-xhmm-cc14e528d909
make
7/15/18
######################################################################################################################################
/usr/bin/ld: cannot find -llapack
collect2: error: ld returned 1 exit status
Makefile:142: recipe for target 'build/execs/xhmm' failed
make[1]: *** [build/execs/xhmm] Error 1
make[1]: Leaving directory '/home/jwaldr/statgen-xhmm-cc14e528d909'
Makefile:121: recipe for target 'all' failed
make: *** [all] Error 2
#######################################################################################################################################














# Rterm.exe --vanilla
# https://cran.r-project.org/web/packages/xhmmScripts/xhmmScripts.pdf
# install.packages("xhmmScripts")
library("xhmmScripts")
