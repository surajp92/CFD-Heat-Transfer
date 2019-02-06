clear
ls
rm *.o *,mod *.exe
ls
gfortran -c *.f95
ls
gfortran *.o -o inviscid_burgers.exe
ls
./inviscid_burgers.exe
