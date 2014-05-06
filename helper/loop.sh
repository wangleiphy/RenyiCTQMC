for Temp in 0.1 0.2 0.3 0.4 0.5
do
    python spinlesstV.py -Temp $Temp  > ../data/L8LA4T${Temp}_exact.dat 
done
