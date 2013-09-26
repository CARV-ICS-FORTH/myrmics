
ENTRY(ar_vectors)

SECTIONS {

  code 0x0 : {

    __linker_start_code = ABSOLUTE(.);

    obj/myrmics/arch/mb/vectors.mb.o   (*)

    . = 0x50;
    
    obj/myrmics/arch/mb/boot.mb.o      (.text)
    *                                  (.text)

    *                                  (.rodata.str1.4)
    *                                  (.rodata.cst4)
    *                                  (.rodata)
    
    __linker_end_code = ABSOLUTE(.);
  }


  data : {
    __linker_start_data = ABSOLUTE(.);

    * (.data)
    * (COMMON)

    __linker_end_data = ABSOLUTE(.);
  }


  /*
    If BSS section is used, it must be zeroed during startup: GCC assumes it
    and will depend on it.
  */
    bss : {
      __linker_start_bss = ABSOLUTE(.);

      * (.bss)

      __linker_end_bss = ABSOLUTE(.);
    }


  /*
    The following are enabled with the -mxl-gp-opt gcc (and linker) flag.
    If used, startup scripts must zero-fill .sbss and .sbss2

      * (.sdata)
      * (.sdata2)
      * (.sbss2)
      * (.sbss)
  */

  }
