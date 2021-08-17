BEGIN {
  mt = 0; 
  m0 = 0; 
  m1 = 0; 
  m2 = 0
} 
($5 == 1) {
  m1 = m1 + 1
} 
($5 == 2) {
  m2 = m2 + 1
} 
{
  m0 = m0 + 1; 
  mt = mt + $5
} 
END { 
  printf "TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair\n"; 
  printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt, m0, m1, m2, m0/mt, m1/m0, m1/m2 
} 