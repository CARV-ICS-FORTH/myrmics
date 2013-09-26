ENTRY(ar_vectors)

SECTIONS {

  /* 
     Read-only code. No read-write variables can be placed here.
     These must be explicitly handled with heap allocation, since
     the combination of GCC and ARM board firmware creates a very big
     mess to manage code-data separation on 1MB pages.

     Specifically, during the early periods, this was managed using a 
     big elf with zero-filled space, because the ARM board firmware 
     complains about multiple ELF code sections and it does NOT fully 
     copy data sections which are farer than a few hundred kilobytes.

     Start with exceptions at offset 0, follow with early boot code
     and put everything else afterwards.
  */
  code 0x0 : {

    __linker_start_code = ABSOLUTE(.);

    obj/mpi/arch/arm/exceptions.arm.o   (*)
    obj/mpi/arch/arm/boot.arm.o         (*)
    *                                   (*)
                     
    __linker_end_code = ABSOLUTE(.);


    /* 
       Microblaze ELF file included here, after it's converted to
       a binary blob. Align everything to 64-B so that it can be safely
       DMAed to all Formic boards.
    */

    . = ALIGN(64);

    __mb_start_include = ABSOLUTE(.);

    INCLUDE "mpi.mb.elf.ld";
    
    . = ALIGN(64);

    __mb_end_include = ABSOLUTE(.);

  }
    

}
