# Clearing
binDir="bin/" ; comando="cd $binDir"; eval $comando
comando="rm *.o"; eval $comando
comando="rm *.mod"; eval $comando

# Compiling modules
comando="gfortran -c ../src/utilities.f90"; eval $comando
comando="gfortran -c ../src/solver.f90"; eval $comando
comando="gfortran -c ../src/objective.f90"; eval $comando
comando="gfortran -c ../src/optimize.f90"; eval $comando

# Linking and compiling driver
comando="gfortran ../src/main.f90 utilities.o solver.o objective.o optimize.o -o ../main.exe"; eval $comando
