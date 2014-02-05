#!/bin/bash

for iconf in {0..699}
do
src=$(printf "OmgOmgwave.000900-%06d.RC16x32_B1830Kud013760Ks013710C17610.it07" $iconf)
target=$(printf "OmgOmgwave.000700-%06d.RC16x32_B1830Kud013760Ks013710C17610.it07" $iconf)
#echo "$src ==> $target"
mv $src $target 
done
