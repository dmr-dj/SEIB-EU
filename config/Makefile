############################################################################
# Makefile for DA
############################################################################

#Define compiler
#  F_COMP=ifort
#   F_COMP=/opt/intel/oneapi/2024.0/bin/mpiifort
   F_COMP=mpifort
   
#Set compile options (Remove needless options after developing phase for )
   #-fast      :������
   #-openmp    :OpenMP�w���s��L����
   #-traceback :�ُ�I�����Ƀ\�[�X�R�[�h�̍s�ԍ���\��
   #-warn all  :�R���p�C�����ɑS�Ă̌x�����b�Z�[�W��L��
   #-assume byterecl : Direct�A�N�Z�X�̎��A"recl"��byte�P�ʂɂ���
   #
   #OPT_COMPILES= -fast -warn all
   OPT_COMPILES= -assume byterecl -O3 -ipo -g -CB -traceback
   
#Link option 1
   #= -openmp :Enable OpenMP expression
   OPT_LINK1= 
   
#Link option 2
   #= -pg :Enable profiler
   OPT_LINK2=
   
#Define execution file name 
   EXE= go.out
   
#Define Source Files
   F_SRC=  modules.f90 \
           start_mpi.f90 \
           main.f90 \
           initialize.f90 \
           metabolic.f90 \
           output.f90 \
           physics.f90 \
           population_regu.f90 \
           spatial_calc.f90 \
           SFLXALL_SRC_VER_2.7.1.f90 \
           etc.f90
   
$(EXE): $(F_SRC)
	$(F_COMP) -c $(OPT_COMPILES) $(?)
	$(F_COMP) $(OPT_LINK1) -o $(EXE)  *.o $(OPT_LINK2)
   
# Definition for make clean
clean:
	rm -f $(EXE) *.o *.L gmon.out *.mod core
