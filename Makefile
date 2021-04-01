#------------------------------------------------------

OPTIONS = -O3 -Wall -DDOUBLE=1 -DVERBOSE=1 -Iinclude
COMPILER  = g++ $(OPTIONS)
PGCOMP  = gcc $(OPTIONS) -I$(PGPLOT_DIR) -L$(PGPLOT_DIR) -L/usr/X11R6/lib64/

LIBRARIES = -lm -Llib -lFiEstAS
#PGLIBS = $(LIBRARIES) -lX11 -lcpgplot -lpgplot -lgcc -lstdc++ -lgfortran -lpng
PGLIBS = $(LIBRARIES) -lcpgplot -lpgplot -lgcc -lstdc++ -lgfortran


#------------------------------------------------------
default: FiEstAS random
#------------------------------------------------------
	#-------------------------------------------#
	# FiEstAS v2.0                              #
	#                                           #
	#         Yago Ascasibar (AIP, Summer 2007) #
	#-------------------------------------------#


#------------------------------------------------------
libs: yTree ySnapshot
#------------------------------------------------------

yTree:
#------------------------------------------------------
	$(COMPILER) -c src/lib/yTree/*.cpp -Iinclude
	ar rs lib/libFiEstAS.a *.o
	rm -f *.o

ySnapshot:
#------------------------------------------------------
	$(COMPILER) -c src/lib/ySnapshot/*.cpp -Iinclude
	ar rs lib/libFiEstAS.a *.o
	rm -f *.o

#------------------------------------------------------
FiEstAS: libs ascii gadget
#------------------------------------------------------

ascii:
#------------------------------------------------------
	$(COMPILER) -o bin/FiEstAS_ASCII src/FiEstAS/FiEstAS_ASCII.cpp $(LIBRARIES)

gadget:
#------------------------------------------------------
	$(COMPILER) $(OPTIONS) -o bin/FiEstAS_Gadget src/FiEstAS/FiEstAS_Gadget.cpp $(LIBRARIES)


random:
#------------------------------------------------------
	$(COMPILER) -o bin/randomData src/FiEstAS/RandomData.cpp


#------------------------------------------------------
plots:
#------------------------------------------------------

f1:
#------------------------------------------------------
	$(PGCOMP) -o bin/fig1 src/plots/paper/Fig1.cpp $(PGLIBS)
# 	./bin/fig1 ../data/Ring/1e2/RandomData.txt d
	./bin/fig1 -print ../data/Ring/1e2/RandomData.txt a
	./bin/fig1 -print ../data/Ring/1e2/RandomData.txt b
	./bin/fig1 -print ../data/Ring/1e2/RandomData.txt c
	./bin/fig1 -print ../data/Ring/1e2/RandomData.txt d

f2:
#------------------------------------------------------
	$(PGCOMP) -o bin/fig2 src/plots/paper/Fig2.cpp $(PGLIBS)
	./bin/fig2 -print

f3:
#------------------------------------------------------
	$(PGCOMP) -o bin/fig3 src/plots/paper/Fig3.cpp $(PGLIBS)
	./bin/fig3 -print


#------------------------------------------------------
clean:
#------------------------------------------------------
	rm -f bin/* lib/*.* *~ */*~ */*/*~ */*/*/*~ */*/*/*/*~ */*/*/*/*/*~


#------------------------------------------------------
backup:
#------------------------------------------------------
	echo -n "mkdir Backup/" > xx
	date "+%y%m%d_%H%M" >> xx
	echo -n "cp -r include/ src/  Backup/" >> xx
	date "+%y%m%d_%H%M" >> xx
	sh xx
	rm xx

#------------------------------------------------------
#                              ... Paranoy@ Rulz! ;^D
#------------------------------------------------------
