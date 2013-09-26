/*-----------------------------------------------------------------------*/
/* Program: Stream                                                       */
/* Revision: $Id: stream_bare.c,v 1.1 2012/05/16 13:17:10 lyberis-spree Exp $ */
/* Original code developed by John D. McCalpin                           */
/* Programmers: John D. McCalpin                                         */
/*              Joe R. Zagar                                             */
/*                                                                       */
/* This program measures memory transfer rates in MB/s for simple        */
/* computational kernels coded in C.                                     */
/*-----------------------------------------------------------------------*/
/* Copyright 1991-2005: John D. McCalpin                                 */
/*-----------------------------------------------------------------------*/
/* License:                                                              */
/*  1. You are free to use this program and/or to redistribute           */
/*     this program.                                                     */
/*  2. You are free to modify this program for your own use,             */
/*     including commercial use, subject to the publication              */
/*     restrictions in item 3.                                           */
/*  3. You are free to publish results obtained from running this        */
/*     program, or from works that you derive from this program,         */
/*     with the following limitations:                                   */
/*     3a. In order to be referred to as "STREAM benchmark results",     */
/*         published results must be in conformance to the STREAM        */
/*         Run Rules, (briefly reviewed below) published at              */
/*         http://www.cs.virginia.edu/stream/ref.html                    */
/*         and incorporated herein by reference.                         */
/*         As the copyright holder, John McCalpin retains the            */
/*         right to determine conformity with the Run Rules.             */
/*     3b. Results based on modified source code or on runs not in       */
/*         accordance with the STREAM Run Rules must be clearly          */
/*         labelled whenever they are published.  Examples of            */
/*         proper labelling include:                                     */
/*         "tuned STREAM benchmark results"                              */
/*         "based on a variant of the STREAM benchmark code"             */
/*         Other comparable, clear and reasonable labelling is           */
/*         acceptable.                                                   */
/*     3c. Submission of results to the STREAM benchmark web site        */
/*         is encouraged, but not required.                              */
/*  4. Use of this program or creation of derived works based on this    */
/*     program constitutes acceptance of these licensing restrictions.   */
/*  5. Absolutely no warranty is expressed or implied.                   */
/*-----------------------------------------------------------------------*/
#include <arch.h>
#include <kernel_toolset.h>

/* INSTRUCTIONS:
 *
 *      1) Stream requires a good bit of memory to run.  Adjust the
 *          value of 'N' (below) to give a 'timing calibration' of 
 *          at least 20 clock-ticks.  This will provide rate estimates
 *          that should be good to about 5% precision.
 */

#ifndef NTIMES
#   define NTIMES       10
#endif

/*
 *      3) Compile the code with full optimization.  Many compilers
 *         generate unreasonably bad code before the optimizer tightens
 *         things up.  If the results are unreasonably good, on the
 *         other hand, the optimizer might be too smart for me!
 *
 *         Try compiling with:
 *               cc -O stream_omp.c -o stream_omp
 *
 *         This is known to work on Cray, SGI, IBM, and Sun machines.
 *
 *
 *      4) Mail the results to mccalpin@cs.virginia.edu
 *         Be sure to include:
 *              a) computer hardware model number and software revision
 *              b) the compiler flags
 *              c) all of the output from the test case.
 * Thanks!
 *
 */

# define HLINE "-------------------------------------------------------------\n"

# ifndef MIN
# define MIN(x,y) ((x)<(y)?(x):(y))
# endif
# ifndef MAX
# define MAX(x,y) ((x)>(y)?(x):(y))
# endif



// ===========================================================================
// ===========================================================================
void stream_check_results(float *a, float *b, float *c, int n) {

        float aj,bj,cj,scalar;
        float asum,bsum,csum;
        float epsilon;
        int     j,k;

    /* reproduce initialization */
        aj = 1.0;
        bj = 2.0;
        cj = 0.0;
    /* a[] is modified during timing check */
        aj = 2.0E0 * aj;
    /* now execute timing loop */
        scalar = 3.0;
        for (k=0; k<NTIMES; k++)
        {
            cj = aj;
            bj = scalar*cj;
            cj = aj+bj;
            aj = bj+scalar*cj;
        }
        aj = aj * (float) (n);
        bj = bj * (float) (n);
        cj = cj * (float) (n);

        asum = 0.0;
        bsum = 0.0;
        csum = 0.0;
        for (j=0; j<n; j++) {
                asum += a[j];
                bsum += b[j];
                csum += c[j];
        }

        //kt_printf ("Results Comparison: \n");
        //kt_printf ("        Expected  : %f %f %f \n",aj,bj,cj);
        //kt_printf ("        Observed  : %f %f %f \n",asum,bsum,csum);


#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif
        epsilon = 1.e-3;

        if (abs(aj-asum)/asum > epsilon) {
                kt_printf ("Failed Validation on array a[]\n");
                kt_printf ("        Expected  : %f \n",aj);
                kt_printf ("        Observed  : %f \n",asum);
        }
        else if (abs(bj-bsum)/bsum > epsilon) {
                kt_printf ("Failed Validation on array b[]\n");
                kt_printf ("        Expected  : %f \n",bj);
                kt_printf ("        Observed  : %f \n",bsum);
        }
        else if (abs(cj-csum)/csum > epsilon) {
                kt_printf ("Failed Validation on array c[]\n");
                kt_printf ("        Expected  : %f \n",cj);
                kt_printf ("        Observed  : %f \n",csum);
        }
        else {
                kt_printf ("Solution Validates\n");
        }
}

// ===========================================================================
// ===========================================================================
int stream_bare(int n) {
    float               *a, *b, *c;

    int                 avgtime[4] = {0}, maxtime[4] = {0},
                        mintime[4] = {0x3FFFFFFF, 0x3FFFFFFF,
                                      0x3FFFFFFF, 0x3FFFFFFF};

    char                *label[4] = {"Copy:      ", "Scale:     ",
                                     "Add:       ", "Triad:     "};

    float               bytes[4] = {
                          2 * sizeof(float) * n,
                          2 * sizeof(float) * n,
                          3 * sizeof(float) * n,
                          3 * sizeof(float) * n
                        };

    int                 BytesPerWord;
    register int        j, k;
    float               scalar;
    int                 t, times[4][NTIMES];

    
    /* --- SETUP --- determine precision and check timing --- */
    kt_printf(HLINE);
    kt_printf("STREAM version $Revision: 1.1 $\n");
    kt_printf(HLINE);
    
    // Assign arrays into kernel heap space (too big to fit in stack)
    a = (float *) mm_va_kernel_base(ar_get_core_id());
    b = a + n;
    c = b + n;
    ar_assert((unsigned int) ((c + n) - a) < MM_KERNEL_SIZE);

    BytesPerWord = sizeof(float);
    //kt_printf("This system uses %d bytes per SINGLE PRECISION word.\n",
    //    BytesPerWord);

    //kt_printf(HLINE);
    kt_printf("Array size = %d, Offset = 0, total memory req. = %.2f KB.\n",
        n, (3.0 * BytesPerWord * (float) n) / 1024.0);
    //kt_printf("Each test is run %d times, but only\n", NTIMES);
    //kt_printf("the *best* time for each is used.\n");


    //kt_printf(HLINE);

    /* Get initial value for system clock. */
    for (j=0; j<n; j++) {
        a[j] = 1.0;
        b[j] = 2.0;
        c[j] = 0.0;
        }

    //if  ( (quantum = checktick()) >= 1) 
    //    kt_printf("Your clock granularity/precision appears to be "
    //        "%d microseconds.\n", quantum);
    //else {
    //    kt_printf("Your clock granularity appears to be "
    //        "less than one microsecond.\n");
    //    quantum = 1;
    //}

    ar_timer_reset();
    for (j = 0; j < n; j++)
        a[j] = 2.0E0 * a[j];
    t = ar_timer_get_cycles();

    //kt_printf("Each test below will take on the order"
    //    " of %d clock cycles.\n", t  );
    //kt_printf("Increase the size of the arrays if this shows that\n");
    //kt_printf("you are not getting at least 20 clock cycles per test.\n");

    kt_printf(HLINE);

    /*  --- MAIN LOOP --- repeat test cases NTIMES times --- */

    scalar = 3.0;
    for (k=0; k<NTIMES; k++)
        {
        ar_timer_reset();

        for (j=0; j<n; j++)
            c[j] = a[j];

        times[0][k] = ar_timer_get_cycles();
        
        ar_timer_reset();

        for (j=0; j<n; j++)
            b[j] = scalar*c[j];

        times[1][k] = ar_timer_get_cycles();
        
        ar_timer_reset();

        for (j=0; j<n; j++)
            c[j] = a[j]+b[j];

        times[2][k] = ar_timer_get_cycles();
        
        ar_timer_reset();

        for (j=0; j<n; j++)
            a[j] = b[j]+scalar*c[j];

        times[3][k] = ar_timer_get_cycles();
        }

    /*  --- SUMMARY --- */

    for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
        {
        for (j=0; j<4; j++)
            {
            avgtime[j] = avgtime[j] + times[j][k];
            mintime[j] = MIN(mintime[j], times[j][k]);
            maxtime[j] = MAX(maxtime[j], times[j][k]);
            }
        }
    
    kt_printf("Function   Rate (B/cc)   Avg time     Min time     Max time\n");
    for (j=0; j<4; j++) {
        avgtime[j] = avgtime[j]/(float)(NTIMES-1);

        kt_printf("%s%11.4f  %13d  %11d  %11d\n", label[j],
               (float) bytes[j] / (float) mintime[j],
               avgtime[j],
               mintime[j],
               maxtime[j]);
    }
    kt_printf(HLINE);

    /* --- Check Results --- */
    stream_check_results(a, b, c, n);
    //kt_printf(HLINE);

    return 0;
}
