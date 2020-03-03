#!/bin/bash
echo "procesando"
python pospro2.py Input-Al-2d* && 
python pospro4.py Input-Pi-2d* && 
python pospro8.py S-2d* && 
python pospro2.py Z0* && 
python pospro4.py ChiUp-2d* && 
python pospro4.py sigma* && 
python pospro8.py ChiS-2d* && 
python pospro4.py ZS-2d* && 
python pospro2.py Alp-2d* && 
python pospro.py Energia* && 
python pospro4.py Ue* && 
python pospro.py Density* && 

echo "Proceso finalizado"



