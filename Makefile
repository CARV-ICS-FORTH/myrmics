
STR_MNT = "[1m[ [36mMNT [0;1m][0m"
STR_CPY = "[1m[ [31mCPY [0;1m][0m"
STR_UMT = "[1m[ [36mUMT [0;1m][0m"


all: myrmics mpi
clean: myrmics_clean mpi_clean


# ===========================================================================
# Myrmics
# ===========================================================================

myrmics:
	@echo
	@echo "        Building Myrmics for MicroBlaze"
	@echo "        -------------------------------"
	@$(MAKE) -f Makefile.myrmics.mb --no-print-directory
	@echo
	@echo "        Building Myrmics for ARM"
	@echo "        ------------------------"
	@$(MAKE) -f Makefile.myrmics.arm --no-print-directory
	@echo

myrmics_clean:
	@echo
	@$(MAKE) -f Makefile.myrmics.mb clean --no-print-directory
	@$(MAKE) -f Makefile.myrmics.arm clean --no-print-directory
	@rm -rf dep/myrmics obj/myrmics


# ===========================================================================
# MPI
# ===========================================================================

mpi:
	@echo
	@echo "        Building MPI for MicroBlaze"
	@echo "        ---------------------------"
	@$(MAKE) -f Makefile.mpi.mb --no-print-directory
	@echo
	@echo "        Building MPI for ARM"
	@echo "        --------------------"
	@$(MAKE) -f Makefile.mpi.arm --no-print-directory
	@echo

mpi_clean:
	@$(MAKE) -f Makefile.mpi.mb clean --no-print-directory
	@$(MAKE) -f Makefile.mpi.arm clean --no-print-directory
	@rm -rf dep/mpi obj/mpi


