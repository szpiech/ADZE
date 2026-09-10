# Makefile for ADZE: Allelic Diversity Analyzer

#G++FLAG = -g
G++FLAG = -O2 -std=c++17
I_PATH = /opt/local/include
L_PATH = /opt/local/lib

#For static linking
#LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a

#For dynamic linking
LINK_OPTS = -L$(L_PATH) -lgsl -lgslcblas

adze-1.0 : ADZE_main.o ADZE_pfile.o ADZE_pop.o ADZE_stats.o ADZE_binom.o\
	ADZE_pbar.o ADZE_main_tools.o
	g++ -o adze-1.0 ADZE_main.o ADZE_pfile.o\
		ADZE_stats.o ADZE_pop.o ADZE_binom.o ADZE_pbar.o\
		ADZE_main_tools.o -I$(I_PATH) $(LINK_OPTS) -lm

ADZE_main.o : ADZE_main.cpp ADZE_main_tools.cpp
	g++ $(G++FLAG) -c ADZE_main.cpp -I$(I_PATH)

ADZE_pfile.o : ADZE_pfile.cpp
	g++ $(G++FLAG) -c ADZE_pfile.cpp

ADZE_stats.o : ADZE_stats.cpp
	g++ $(G++FLAG) -c ADZE_stats.cpp

ADZE_pop.o : ADZE_pop.cpp
	g++ $(G++FLAG) -c ADZE_pop.cpp

ADZE_binom.o : ADZE_binom.cpp
	g++ $(G++FLAG) -c ADZE_binom.cpp

ADZE_pbar.o : ADZE_pbar.cpp
	g++ $(G++FLAG) -c ADZE_pbar.cpp

ADZE_main_tools.o : ADZE_main_tools.cpp
	g++ $(G++FLAG) -c ADZE_main_tools.cpp -I$(I_PATH)

clean :
	rm ADZE_main.o
	rm ADZE_pfile.o
	rm ADZE_stats.o
	rm ADZE_pop.o
	rm ADZE_binom.o
	rm ADZE_main_tools.o
	rm ADZE_pbar.o
