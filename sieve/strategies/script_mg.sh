#!/usr/bin/env bash
#[fbb_lpb]
fbb=$1
lpb=$2


cat "Data_ECM-B12["$fbb"_"$lpb"]_["* > "Data_ECM-B12["$fbb"_"$lpb"]_all"
./gfm -fch -fch_in "Data_ECM-B12["$fbb"_"$lpb"]_all" -fch_out "Data_ECM-B12["$fbb"_"$lpb"]"
mkdir "Data_ECM-B12_"$fbb"_"$lpb
mv "Data_ECM-B12["$fbb"_"$lpb""* "Data_ECM-B12_"$fbb"_"$lpb

cat "Data_ECM-M16["$fbb"_"$lpb"]_["* > "Data_ECM-M16["$fbb"_"$lpb"]_all"
./gfm -fch -fch_in "Data_ECM-M16["$fbb"_"$lpb"]_all" -fch_out "Data_ECM-M16["$fbb"_"$lpb"]"
mkdir "Data_ECM-M16_"$fbb"_"$lpb
mv "Data_ECM-M16["$fbb"_"$lpb""* "Data_ECM-M16_"$fbb"_"$lpb

cat "Data_ECM-M12["$fbb"_"$lpb"]_["* > "Data_ECM-M12["$fbb"_"$lpb"]_all"
./gfm -fch -fch_in "Data_ECM-M12["$fbb"_"$lpb"]_all" -fch_out "Data_ECM-M12["$fbb"_"$lpb"]"
mkdir "Data_ECM-M12_"$fbb"_"$lpb
mv "Data_ECM-M12["$fbb"_"$lpb""* "Data_ECM-M12_"$fbb"_"$lpb

cat "Data_PM1["$fbb"_"$lpb"]_["* > "Data_PM1["$fbb"_"$lpb"]_all"
./gfm -fch -fch_in "Data_PM1["$fbb"_"$lpb"]_all" -fch_out "Data_PM1["$fbb"_"$lpb"]"
mkdir "Data_PM1_"$fbb"_"$lpb
mv "Data_PM1["$fbb"_"$lpb""* "Data_PM1_"$fbb"_"$lpb

cat "Data_PP1-27["$fbb"_"$lpb"]_["* > "Data_PP1-27["$fbb"_"$lpb"]_all"
./gfm -fch -fch_in "Data_PP1-27["$fbb"_"$lpb"]_all" -fch_out "Data_PP1-27["$fbb"_"$lpb"]"
mkdir "Data_PP1-27_"$fbb"_"$lpb
mv "Data_PP1-27["$fbb"_"$lpb""* "Data_PP1-27_"$fbb"_"$lpb

cat "Data_PP1-65["$fbb"_"$lpb"]_["* > "Data_PP1-65["$fbb"_"$lpb"]_all"
./gfm -fch -fch_in "Data_PP1-65["$fbb"_"$lpb"]_all" -fch_out "Data_PP1-65["$fbb"_"$lpb"]"
mkdir "Data_PP1-65_"$fbb"_"$lpb
mv "Data_PP1-65["$fbb"_"$lpb""* "Data_PP1-65_"$fbb"_"$lpb
