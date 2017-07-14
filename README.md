# ComputationalPhysicsFortran
Fortran code that accompanies the book "Computational Physics: An Undergraduates Guide"

This replaces the code links that were found [here](https://compphysintro.wordpress.com/)
 
If you're following the book you'll be using the Cygwin terminal on Windows, with the Emacs text editor and be compiling with the gnu Fortran compiler. If your using Unix/Linux/Mac then you'll be using whatever shell - editor combo you prefer. 

To compile a stand alone program you type at the command prompt:

gfortran programName.f90 -o outputName

Obviously 'progName' and 'outputName' are just placeholders for the actual filename and given executable name respectively. If you've aliased the compiler as per the instructions in the book replace gfortran with gfc. Make sure to invoke the compiler from within the same directory as your saved FORTRAN file (or provide the correct path).

For a program that uses a module you type:

gfortran moduleName.f90 programName.f90 -o outputName

For programs that require external libraries you add to the end of the line -llibraryName1 -llibraryName2 and so on. To be clear that is the library name prefixed with '-l'. For example, to use LAPACK subroutines that also uses BLAS you'd add -llapack -lblas to the end of the command above.

To run the programs you have compiled type at the command prompt:

./outputName

Enjoy!

NB: There is a "typo" in the Runge-Kutta-Fehlberg algorithms. One of the coefficients is wrong in the code (this was supposed to be an exercise in the book but was omitted by mistake). The coefficients in the book are correct. To help here is the wikipedia article on the [Runge-Kutta-Fehleberg algorithm](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) - pay attention to the Butcher Tableau. 

