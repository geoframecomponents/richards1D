set table "fortran_doc.pgf-plot.table"; set format "%.5f"
set format "%.7e";; set samples 100; set dummy x,y; plot [x=0:8] 0.5*(1+(erf((x-4)/(1*sqrt(2)))));
