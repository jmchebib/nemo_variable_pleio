## nemo_evolving_pleio
individual-based, genetically explicit, forward-in-time, stochastic, population modelling software that allows for evolving pleiotropic connections between alleles and traits as well as evolving correlational mutations on those traits (based on Guillaume and Rougement 2006) 

#see home page and nemo2 manual for instructions on use. http://nemo2.sourceforge.net/

# Step 0
#requires gsl which can be installed on the linux command line from http://www.gnu.org/software/gsl/

sudo apt-get install gsl

#but depends on your system

sudo apt-get libgsl0-dev

# Step 1
#download / clone source files from this repository into src/ directory

# Step 2
#create VERSION file outside src/ directory

echo 2.3.x > VERSION

# Step 3 
#create bin/directory next to src/ directory

mkdir bin

# Step 4
#move Makefile outside of src/ and bin/ with the following (Makefile may require editing depending on system and location of gsl library):

mv src/Makefile ./

# Step 4b
#move this README.md outside of src/ if it is in there

mv src/README.md ./

# Step 5
#make nemo program file in bin directory by using 'make' command outside of bin/ and src/ directories

make

# Step 6
#use newly created nemo program to run simulations along with init files (see nemo2 manual)

./bin/nemo2.3.x-locked init_file.ini



