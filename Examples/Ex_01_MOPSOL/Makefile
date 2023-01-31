# Makefile MOPSO_Light

CC=gfortran
obj = Precision_defaults_MOPSO_mod.o ObjectiveFunction_MOPSO_mod.o ParetoFront_mod.o File_Manager_ParetoFront_mod.o ParetoFront_EDSD_mod.o SobolSequence_mod.o Swarm_MOPSO_mod.o MOPSO_Light_mod.o


MOPSO_Light_make:
	$(CC) -c ./../../libs/Precision_defaults_MOPSO_mod.f03
	$(CC) -c ./../../libs/ObjectiveFunction_MOPSO_mod.f03
	$(CC) -c ./../../libs/ParetoFront_mod.f03
	$(CC) -c ./../../libs/File_Manager_ParetoFront_mod.f03
	$(CC) -c ./../../libs/ParetoFront_EDSD_mod.f03
	$(CC) -c ./../../libs/SobolSequence_mod.f03
	$(CC) -c ./../../libs/Swarm_MOPSO_mod.f03
	$(CC) -c ./../../libs/MOPSO_Light_mod.f03
	$(CC) -c OF_test1_mod.f03
	$(CC) -c MOPSO_Light_main.f03
	$(CC) -o test1.a $(obj) OF_test1_mod.o MOPSO_Light_main.o

