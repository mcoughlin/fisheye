#!/bin/bash
# to run sqlite3 -separator "," all_sky_sqlite.db < injest_allsky_hp.sql
# just turn this into a bash script?

rm all_sky_sqlite.db
sqlite3  all_sky_sqlite.db "CREATE TABLE medskybrightness (hpindex INT, R REAL, G REAL, B REAL, airmass REAL, mjd REAL);"

for i in $( ls -d output/ut* ); do 
    echo $i
    sqlite3 -separator ","  all_sky_sqlite.db ".import $i/healmaps.dat medskybrightness"
done

sqlite3  all_sky_sqlite.db "CREATE INDEX hpindex_index on medskybrightness (hpindex);"
sqlite3  all_sky_sqlite.db "CREATE INDEX mjd_index on medskybrightness (mjd);"
sqlite3  all_sky_sqlite.db "CREATE INDEX airmass_index on medskybrightness (airmass);"

