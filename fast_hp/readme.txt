
commands for generating sqlite database of the output


CREATE TABLE medskybrightness (hpindex INT, R REAL, G REAL, B REAL,
 airmass REAL, mjd REAL);

.separator ","
.import output/ut012916/healmaps.dat medskybrightness



CREATE INDEX hpindex_index on medskybrightness (hpindex);
CREATE INDEX mjd_index on medskybrightness (mjd);
CREATE INDEX airmass_index on medskybrightness (airmass);
