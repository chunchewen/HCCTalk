* Orabanche Beta-Binomial.sas
* Fit bb model to orbanche data
* Infile: I:\\Brian\\LDA\\Lectures\\Lecture 16 Quasi-Likelihood and OD\\orabanche.txt
* March 28, 2024
********************;
data one;
 infile "I:\\Brian\\LDA\\Lectures\\Lecture 16 Quasi-Likelihood and OD\\orabanche.txt" firstobs=2;
 input y n cuke oa75 int;
run;
proc print data=one; run;

proc fmm data=one;
 model y/n= cuke oa75 int/dist=bb ;
run;
