REM
REM General batch for Mdelta, Vdelta 
REM Call with 3 arguments:
REM %1 = networks path for space
REM %2 = output file full path
REM %3 = s (float)
REM %4 = Number of spaces (normally 1)
REM
REM Naming convention for output:
REM results\<space>.<Nalpha>.delta.csv
REM
REM e.g.  for a panmictic 128 and Nalpha=8 
REM         results\pan_128.8.delta.csv
REM
REM Example call:
REM delta.bat networks/pan64/ results/pan64.delta.csv 0.02 1
REM
setlocal enableDelayedExpansion
@echo on
echo "Starting Time:" !TIME!
REM p,q,Nialpha,Nibeta, Nc runs
delta 0.05 %3 1000000 %1 %4 1 > %2
echo "Starting N=64, Time:" !TIME!
FOR %%p IN (0.1 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.9 0.95) DO (
echo  %%p
delta %%p %3 1000000 %1 %4 0 >> %2
) 
echo "Finished Time:" !TIME!