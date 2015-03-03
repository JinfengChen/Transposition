echo "Correlation of calls from Relocate, temp and Relocatei. Correlation of calls and read depth, insertion size"
cat mPing_Call_Correlation.R | R --slave
echo "Summary of how many additional calls were in Relocatei compared to Relocate"
cat AdditionRelocaTEi.R | R --slave

