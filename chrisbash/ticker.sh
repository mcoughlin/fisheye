#!/bin/bash
while [ 1 ] ; do
echo M; ls M/*.phot | wc
echo B; ls B/*.phot | wc
echo G; ls G/*.phot | wc
echo R; ls  R/*.phot | wc
sleep 30
done
