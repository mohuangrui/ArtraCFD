/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 * - Current file only includes structure type declarations and function    *
 *   prototypes that are shared and accessed by most modular files.         *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 * - Header files should be once-only headers. A standard way to prevent 
 *   include a header file more than once is to enclose the entire area
 *   contents of the file in a conditional.
 * - Do not start the guard symbol with an underscore. Leading underscore 
 *   names are reserved for internal use by the C implementation, breaking
 *   this rule can cause unnecessary and very puzzling errors.
 * - To guarantee uniqueness, the format of the name should be
 *   <PROJECT>_<PATH>_<FILE>_H_
 ****************************************************************************/
#ifndef ARTRACFD_COMMONS_H_ /* if this is the first definition */
#define ARTRACFD_COMMONS_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 * - Every header file should explicitly #include every other header file 
 *   that current header file requires to compile correctly, but no more,
 *   for instance, do not include header files that only the .c file needs. 
 * - The #include directive works by directing the preprocessor to scan the 
 *   specified file as input before continuing the current source file.
 * - #include <***.h> for system header files. It searches for a named file 
 *   in a standard list of system directories.
 * - #include "***.h" for header files of your own program. It searches for 
 *   a named file in the directory containing the current file, such as 
 *   #include "include/foo.h".
 * - Names and order of includes: use standard order for readability and 
 *   to avoid hidden dependencies: 
 *   Related header "foo.h" of the "foo.c" file
 *   C library
 *   C++ library
 *   Other libraries' .h
 *   Your project's .h.
 ****************************************************************************/
#include <stdio.h> /* standard library for input and output */
/****************************************************************************
 *
 *                         Rules for Programming
 *
 * - Always immediately initialize variables in their declaration. Use 0 for
 *   integers, 0.0 for reals, NULL for pointers, '\0' for chars, {'\0'} for
 *   string arrays.
 * - Put one variable declaration per line, and comment them.
 * - Avoid global variables. Use const declarations instead of #define 
 *   wherever possible.
 * - Be aware of pointers. Always keep in mind that assigning a valid memory
 *   location to the pointer before dereferencing the pointer. Or there will
 *   be a segmentation error.
 * - Be aware of string manipulations, most string functions will have
 *   undefined behavior if memory locations of input objects overlap.
 * - Assignment operators always have spaces around them. x = 0;
 *   Other binary operators usually have spaces around them, v = w * x;
 *   Parentheses should have no internal padding. v = w * (x + z);
 *   No spaces separating unary operators and their arguments. x = -5.
 * - Minimize use of vertical whitespace. It's more a principle than a rule:
 *   don't use blank lines when you don't have to. Blank lines at the 
 *   beginning or end of a function very rarely help readability.
 * - When defining a function, parameter order is: inputs, then outputs.
 * - Always return or exit a function with the 0 or 1 code when all of the 
 *   commands executed successfully or unsuccessfully.
 * - Avoid side effects, such as do not use assignment statements in if 
 *   condition, should use ++ and -- on lines by themselves.
 * - Always use the prefix version of ++ and -- (++x, --x) instead of the
 *   postfix version (x++, x--).
 * - Remember that there are only two precedence levels to remember in 
 *   C: multiplication and division come before addition and subtraction.
 *   Everything else should be in parentheses. 
 * - Avoid complex logic like multiply nested ifs. Consider splitting your 
 *   code into multiple procedures, to decrease the level of complexity.
 * - Variables should be declared as locally as possible:
 *    * declare non-constant variables that are used through-out the function
 *      at the top.
 *    * declare constant variables when their values can be determined.
 *    * declare variables that are used in only a local scope of the function
 *      at the point where they are needed, e.g., loop counts.
 * - Use of const whenever possible, using const as much as possible is
 *   compiler-enforced protection from unintended writes to data
 *   that should be read-only.
 *    * Declare variables that should not be changed after initialization.
 *          const double pi = 3.1415; const data;
 *    * Two ways of declaring a const pointer: 
 *      a) the target address which the pointer points to is fixed, but 
 *      the content in the address can be changed:
 *      double * const pointer; (pointer itself is const, can't be changed)
 *      Use const pointer to a large compound variable type is useful for
 *      storage that can be changed in value but not moved in memory, since
 *      address of the storage is fixed, then the pointer can be const.
 *      b) the pointer itself can be changed but the data it points to 
 *      can not be changed.
 *      const double *pointer; (const data but non-const pointer)
 *      Large compound user-defined variable types ('structures' in C and
 *      'classes' in C++) should always be passed by reference or pointer
 *      instead of as a copy. Use a pointer point to const value can pass
 *      the variable without value copying and also prevent the value
 *      being altered.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                          Rules for Naming
 *
 * - The most important consistency rules are those that govern naming. 
 *   The style of a name immediately informs us what sort of thing the
 *   entity is: a type, a variable, a function, a constant, a macro
 *   without requiring us to search for the declaration of that entity.
 *   The pattern-matching engine in our brains relies a great deal on 
 *   these naming rules. 
 * - Function names, variable names, and filenames should be descriptive;
 *   Eschew abbreviation.
 * - Filenames should be all lowercase and can include underscores.
 * - The names of all types - classes, structs, typedefs, and enums have 
 *   the same naming convention: Type names should start with a capital
 *   letter and have a capital letter for each new word. No underscores.
 *   For example, MyExcitingStruct, MyExcitingEnum.
 * - The names of variables and data members are all begin with lowercase
 *   have a capital letter for each new word, with or without underscores.
 * - The names of vectors and tensors start with a capital letter and only 
 *   have a single descriptive word.
 * - Functions should start with a capital letter and have a capital letter
 *   for each new word. No underscores.
 * - Accessors and mutators (get and set functions) should match the name
 *   of the variable they are getting and setting.
 *   my_exciting_member_variable(), set_my_exciting_member_variable(). 
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                         Tips for C Programming I
 *
 *                            Information Flow
 *
 * - Elements of a Program: In computer programming, you need two things: 
 *   data (variables) and instructions (code or functions). Variables are
 *   the basic building blocks of a program. Instructions tell the computer
 *   what to do with the variables.
 * - Program should read like an essay. It should be as clear and easy to
 *   understand as possible. Always comment your programs, which helps you 
 *   organize your thoughts while making you work an art rather than a junk.
 * - A program should be as simple as possible. A programmer should avoid
 *   clever tricks. For example, there are 15 operator precedence rules in C.
 *   These can be simplified into two rules:
 *   1. Multiply and divide come before add and subtract.
 *   2. Put parentheses around everything else.
 * - Structured Programming:
 *   Top-down programming -- start at the top (main) and work your way down.
 *                           The main function acts as the topmost outline.
 *   Bottom-up programming -- write the lowest-level function first, testing
 *                            it and then building on that working set. 
 *                            bottom-up techniques is very useful when
 *                            working with a new and unfamiliar function.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                         Tips for C Programming II
 *
 *                            Modular Programming
 *
 * Several schools of code design exist. In structured programming, you divide
 * the code into modules, then divide the modules into submodules, then 
 * divide the sub-modules into subsubmodules, and so on. 
 *
 * "Modular Programming" is the act of designing and writing programs as
 * interactions among functions that each perform a single well-defined
 * function, and which have minimal side-effect interaction between them. Put
 * differently, the content of each function is cohesive, and there is low
 * coupling between functions. 
 *
 * "Modular Programming" tends to encourage splitting of functionality into two
 * types: "Manager" functions control program flow and primarily contain calls
 * to "Worker" functions that handle low-level details, like moving data
 * between structures. 
 *
 * All have the same basic principle at heart: "Arrange the program's
 * information in the clearest and simplest way possible, and then try to turn
 * it into C code."
 *
 * Information is a key part of any program. The key to any program is 
 * deciding what information is being used and what processing you want to
 * perform on it. Information flow should be analyzed before the design begins.
 *
 * Information hiding is key to good programming. A module should make public
 * only the minimum number of functions and data needed to do the job. The 
 * smaller the interface, the simpler the interface. The simpler the interface
 * the easier it is to use. Also, a simple interface is less risky and less error
 * prone than a complex one.
 *
 * Small, simple interfaces are also easier to design, test, and maintain. Data
 * hiding and good interface design are key to making good modules
 *
 * General Module Design Guidelines
 * - The number of public functions in a module should be small.
 * - The information passed between modules should be limited.
 * - All the functions in a module should perform related jobs.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *               More Details about Modular Programming
 *
 * - Structure code by modules.
 * - A module is a collection of data and functions that perform related tasks.
 * - Modules should be designed to minimize the amount of information that has
 *   to pass between them.
 * - Modules are divided into two parts: public and private.
 * - The public part tells the user how to call the functions in the module.
 *   These public information that is shared between modules are put in a 
 *   header file. So, what goes in a header file is precisely that information
 *   that needs to be communicated. Do not include definitions of structures or
 *   prototypes of functions that are only used in current module source file;
 *   those private information of the module should go within the source file.
 * - As a conclusion, only the declarations, definitions and prototypes needed
 *   to use a module should appear in its header file, this file is always 
 *   used to access the module, it's true even for the source of the module. 
 * - The header should contain all the public information, such as: 
 *   * A comment section describing clearly what the module does and what
 *     is available to the user.
 *   * Common structure definitions.
 *   * Prototypes of all the public functions.
 *   * extern declarations for public variables.
 * - Anything that is internal to the module is private, such as all the
 *   implementations (functions and data) of the module. 
 * - Everything that is not directly usable by the outside world should be 
 *   kept private,
 * - Private information should be put in the source file of the module.
 * - Private functions that will not be called from outside the module should
 *   be declared static. Variables declared outside of a function that are 
 *   not used outside the module should be declared static.
 * - Obviously, the prototypes for static functions should not be put in the
 *   module's header file.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                            C Puzzles
 *
 * - C uses void for two purposes:
 *   In a function declaration, void indicates that the function returns
 *   no value or takes no arguments.
 *   In a pointer declaration, void defines a generic pointer.
 * - static has two meanings:
 *   For function or global variable, static means "private to this file."
 *   If declare a function or global variable as static, it becomes internal.
 *   That is, they are global to only the C file they exist in, not outside.
 *   You cannot access the function or variable through the extern keyword 
 *   from other files in your project. Note: extern keyword is default for
 *   any functions declared without the keyword "static".
 *   For data defined inside a function, it means "variable is allocated from
 *   static memory (instead of temporary stack)." When you declare a local
 *   variable as static, it is created just like any other variable. However,
 *   when the variable goes out of scope (i.e. the block it was local to is
 *   finished) the variable stays in memory, retaining its value. The variable
 *   stays in memory until the program ends. While this behaviour resembles
 *   that of global variables, static variables still obey scope rules and
 *   therefore cannot be accessed outside of their scope. Hence, you need to
 *   pass its address out if access is needed outside its local scope.
 * - In C single quotes identify a single character, while double quotes 
 *   create a string literal. 'a' is a single a character literal, while "a" 
 *   is a string literal containing an 'a' and a null terminator (that is a
 *   2 char array). Note that in C, the type of a character literal is int,
 *   and not char, that is sizeof 'a' is 4 in an architecture where ints 
 *   are 32bit (and CHAR_BIT is 8), while sizeof(char) is 1 everywhere.
 * - Keywords such as const are very useful for optimizing compilers to
 *   produce efficient code. This is particularly true when used in combination
 *   with the restrict keyword from the C99 standard. The const keyword
 *   means that the indicated variable cannot be assigned to. This is 
 *   particularly important for pointers and arrays: to pass an array so that
 *   it cannot be changed by a routine. However, the const keyword cannot
 *   prevent changes due to aliasing. That is, the same memory location can
 *   be referred to through different pointers. This has an effect on
 *   optimization of the code: since the compiler cannot be sure that whether
 *   the const variable will be changed after aliasing, it has to assume that
 *   it can be changed, which prevents the desired optimizations. The restrict
 *   keyword in the C99 standard for C is a way of telling a compiler
 *   to assume that aliasing does not occur for a particular variable. That is,
 *   if a variable is declared restrict, the compiler assumes that changes to
 *   other variables cannot affect that variable. Whether aliasing does or 
 *   does not occur is the responsibility of the programmer.
 * - Definition VS Declaration: definition occurs in only one specifies the 
 *   type of an object; reserves storage for it; is used to create place
 *   new objects, example: int my_array[100]; declaration can occur multiple
 *   describes the type of an object; is used to refer to objects defined.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                  Prototypes of Involved Standard Functions
 *
 * int printf(const char *format, ...);
 * char *strncpy(char *restrict s1, const char *restrict s2, size_t n);
 * int strncmp(const char *s1, const char *s2, size_t n);
 * int strcmp(const char *s1, const char *s2);
 * char *fgets(char *s, int n, FILE *stream);
 * int sscanf(const char *s, const char *format, ...);
 * int snprintf(char *str, size_t size, const char * restrict format, ...)
 * size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
 * size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
 * void free(void *ptr);
 * void *malloc(size_t size);
 *
 ****************************************************************************/
/****************************************************************************
 * 
 *                           Compile with make
 *
 * - http://www.gnu.org/software/make/manual/make.html
 * - Use make program to compile and link programs.
 * - Turn on all the warning flags, then make your program warning free.
 * - When recompiles the program, each changed C source file must be recompiled.
 * - If a header file has changed, each C source file that includes the header
 *   file must be recompiled to be safe.
 * - Each compilation produces an object file corresponding to the source file.
 * - Finally, if any source file has been recompiled, all the object files
 *   whether newly made or saved from previous compilations, must be linked 
 *   together to produce the new executable program.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                     Code Performance and Optimization
 *
 * - Rule of thumb, let the compiler do the job. Code should compiled with 
 *   compiler optimizations. ICC or GCC are the best.
 * - High performance coding requires understanding modern computer hardware. 
 *   The most crucial concept in all this is that of a memory hierarchy.
 *   Optimizing compilers know more about the target architecture(s) than
 *   most programmers, so we should write code that the compiler can make 
 *   best use of.
 * - keep your code general instead of obsessively optimizing your code with 
 *   awkward code structure. Hand optimizations can create odd looking code
 *   that is harder for a compiler to match up to an optimization template.
 *   It is often better to resist the temptation to optimize the code. 
 * - Code performance is not only about algorithms, but also is related to
 *   CPU time and memory reading. To get close to CPU peak, codes should
 *   designed to make best use of hardware, especially memory caches.
 * - Note that it is usually not important to make every routine optimal. 
 *   Often only a small fraction of the code in a large system has a 
 *   significant impact on efficiency. This is sometimes expressed by the
 *   slogan that 95% of the time is spent in 5% of the code. If there is a
 *   bottleneck in your code, have a close look at it. Ask yourself:
 *   could better algorithms be used? Could it be implemented more 
 *   efficiently? Should a higher level of compiler optimization be used?
 *   Should the data structures be re-designed? Repeat the process until 
 *   the system is performing as expected.
 * - Google gperftools can do program performance checking including 
 *   heap-checker, heap-profiler and cpu-profiler.
 * - Valgrind can be used for memory access check and debugging.
 * - Premature optimization is the root of all evil. By using modern compilers,
 *   you do NOT need to concern about:
 *   * Register allocation. Assign commonly used variables to registers for
 *     rapid access.
 *   * Common sub-expression elimination. If an expression appears several 
 *     times, evaluate it once and store the result.
 *   * Loop transformations. Re-order loops to avoid inefficiencies.
 * - But, There are things that you should optimize at a low level:
 *   * Temporal locality: Nearby memory accesses in time should be to nearby
 *     locations in memory. Accessing far-apart memory locations means that 
 *     each time a new memory location is accessed, memory within the CPU 
 *     has to be filled with values at and around that memory location.
 *     C stores its arrays in row-major order. That is, for array a[j][i],
 *     consecutive memory locations hold a[0][0], a[0][1], a[0][2], . . . 
 *     To keep our memory references close together, we should make: the later
 *     position the index in the array, the inner position the index in the
 *     loop. That is, i is the inner loop for a[j][i].
 *   * Memory usage: Try to re-use dynamically allocated memory. This is not
 *     only helpful for avoiding memory leaks, but also avoids time allocating
 *     and freeing memory.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                   Issues Related to  Numerical Computing
 *
 * - Rise the concern about numerical accuracy and reliability whenever conduct
 *   an algorithm or even an operation. Be aware of catastrophic cancellation
 *   in operations, numerical stability of the algorithms.
 * - Recommend double rather than float type for floating point variables. 
 * - Don't test for exact equality between floating point numbers. Don't do
 *   this even if one was assigned to the other: y = x;...if ( x == y )...
 * - Don't subtract nearly equal quantities and then divide by something small.
 *   This often results in catastrophic cancellation and all digits of accuracy
 *   are lost. In general, if you subtract numbers where the first k digits are
 *   equal, you lose k digits of accuracy.
 * - Avoid using floating point numbers as loop counters if exact loop
 *   behaviors are required. Round off errors are unreliable.
 * - Priorities in writing scientific software should be
 *   * correctness,
 *   * numerical stability,
 *   * accurate discretization (including estimating accuracy),
 *   * flexibility,
 *   * efficiency (speed and memory).
 *
 ****************************************************************************/
/****************************************************************************
 * 
 *                             C Textbooks
 *
 * - Reference Style
 *   The C Programming Language (Second edition) - Brian W. Kernighan and
 *   Dennis M. Ritchie (K&R does not address good program design nor good
 *   programming practice)
 * - Beginner
 *   C wikibook http://en.wikibooks.org/wiki/C_Programming
 *   Practical C Programming, 3rd Edition - Steve Oualline
 *   C Programming: A Modern Approach - K. N. King
 * - Intermediate
 *   Algorithms in C - Robert Sedgewick
 * - Above Intermediate
 *   Expert C Programming: Deep C Secrets - Peter van der Linden
 * - Software engineering
 *   Writing Scientific Software: a guide to good style, by Suely Oliveira
 *   and David Stewart
 *
 ****************************************************************************/
/****************************************************************************
 * Data Structure Declarations
 * - Data structures with heterogeneous elements (i.e., elements of different
 *   types) can be defined as struct in C. Combining data used for a common
 *   purpose into a single data structure can provide some level of 
 *   abstraction which can simplify interfaces and other routines.
 * - Use struct to pass a bunch of data at a time, it's simple and elegant.
 *   However, make sure about these:
 * - Always pass structures by reference. That is, use pointers to structures
 *   as function arguments even when nothing in the struct will be modified in
 *   the function (at this circumstance, const modifier should be used). 
 *   Should never do value passing to avoid copying the complete contents
 *   of the structure onto the stack.
 * - Make sure assigning a valid memory location to the pointer before
 *   dereferencing a pointer!
 ****************************************************************************/
/*
 * Define some universe data type for portability and maintenance
 */
typedef double Real;
/*
 * Field variables of flow
 *
 * Conservative variables are vectors with five elements(rho, rho_u, rho_v,
 * rho_w, rho_eT), while each element is a three dimensional array in 3D space.
 * Thus, the conservative variables need to be presented as a 4 dimensional
 * array in 3D flow.
 *
 * Using high order pointers (arrays) is complicated. The pointers in the 
 * array of arrays of arrays waste space and the malloc( ) calls are 
 * expensive, and it is also very time-consuming for nested loops 
 * because of causing lots of cache misses.
 *
 * Maintaining a multi-dimensional array within a single linear array is a
 * common performance technique. High-performance code instead implements a
 * multi-dimensional array as a single linear array with hand-authored array
 * indexing math to keep track of what values are where:
 *
 * value = data[ k * height * depth + j * depth + i ];
 *
 * Since the array is a single large chunk of memory, sweeping through it from
 * start-to-finish creates a regular access pattern that processor prefetchers
 * easily recognize, which enables them to load caches in the background. The
 * result is fewer cache misses and much better performance.
 *
 * This code use nested loops with linear array and index math instead nested
 * loops with multi-dimensional arrays, because we need to store five
 * conservative variables at each (k,j,i), the most efficient index arrangement
 * should be the following:
 *
 * U[kMax][jMax][iMax][5] (Note: NOT U[5][kMax][jMax][iMax]);
 * int k, j, i;
 * int idx; use long type if needed
 * idx = ((k * jMax + j) * iMax + i) * 5;
 * rho    = U[idx+0];
 * rho_u  = U[idx+1];
 * rho_v  = U[idx+2];
 * rho_w  = U[idx+3];
 * rho_eT = U[idx+4];
 */
typedef struct {
    Real *Un; /* store the "old" field data for intermediate calculation */
    Real *U; /* store updating field data, and updated data after every computation  */
    Real *Uswap; /* an auxiliary storage space */
}Field;
/*
 * Space domain parameters
 */
typedef struct {
    int nx; /* mesh number in x */
    int ny; /* mesh number in y */
    int nz; /* mesh number in z */
    int ng; /* number of layers of ghost cells */
    int iMax; /* total node number in x */
    int jMax; /* total node number in y */
    int kMax; /* total node number in z */
    int nMax; /* total node number */
    Real dx; /* mesh size in x */
    Real dy; /* mesh size in y */
    Real dz; /* mesh size in z */
    Real ddx; /* reciprocal of mesh size in x */
    Real ddy; /* reciprocal of mesh size in y */
    Real ddz; /* reciprocal of mesh size in z */
    Real xMin; /* coordinates define the space domain */
    Real yMin; /* coordinates define the space domain */
    Real zMin; /* coordinates define the space domain */
    Real xMax; /* coordinates define the space domain */
    Real yMax; /* coordinates define the space domain */
    Real zMax; /* coordinates define the space domain */
    int *nodeFlag; /* node type integer flag: normal, ghost, solid, etc. */
}Space;
/*
 * Particle Entities
 */
typedef struct {
    int totalN; /* total number of particles */
    Real *headAddress; /* record the head address */
    Real *x; /* x coordinates of the particle center */
    Real *y; /* y coordinates of the particle center */
    Real *z; /* z coordinates of the particle center */
    Real *r; /* radius of the particle */
    Real *u; /* x component of the velocity of particles */
    Real *v; /* y component of the velocity of particles */
    Real *w; /* z component of the velocity of particles */
}Particle;
/*
 * Time domain parameters
 */
typedef struct {
    int restart; /* restart flag */
    Real totalTime; /* total evolution time */
    Real currentTime; /* current time recorder */
    Real dt; /* time step size */
    Real numCFL; /* CFL number */
    int totalStep; /* total number of steps */
    int stepCount; /* step number count */
    int totalOutputTimes; /* total times of exporting computed data */
    int outputCount; /* exporting data count */
}Time;
/*
 * Flow properties and physics parameters
 */
typedef struct {
    Real refMa; /* reference Mach number */
    Real refMu; /* reference dynamic viscosity for Sutherland's law */
    Real refPr; /* reference Prandtl number */
    Real gamma; /* heat capacity ratio */
    Real gasR; /* the gas constant */
    Real cv; /* specific heat capacity at constant volume */
    Real delta; /* numerical dissipation */
    Real refLength; /* characteristic length */
    Real refDensity; /* characteristic density */
    Real refVelocity;  /*characteristic velocity */
    Real refTemperature; /* characteristic temperature */
    int probe[12]; /* store various information of probes */
    Real probePos[11][6]; /* store position information of probes */
}Flow;
/*
 * Domain partition structure
 */
typedef struct {
    int totalN; /* total number of domain partitions */
    int subN; /* inner partitions for each partition */
    int kSub[13]; /* inner decomposition control for each partition */
    int kSup[13];
    int jSub[13];
    int jSup[13];
    int iSub[13];
    int iSup[13];
    int normalZ[7]; /* outer surface normal vector of each inner part */
    int normalY[7];
    int normalX[7];
    int typeBC[7]; /* BC type of each inner part */
    Real valueBC[7][6]; /* BC values of each inner part */
    char name[13][15]; /* store names of each inner part */
    int typeIC[11]; /* list structure for recording regional initial conditions */
    Real valueIC[11][15]; /* queue data structure for storing regional initial values */
}Partition;
/*
 * Program command line arguments and overall control
 */
typedef struct {
    char runMode; /* mode: [s] serial, [i] interact, [t] threaded, [m] mpi, [g] gpu */
    int processorN; /* number of processors */
}Control;
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Command line processor
 *
 * Parameters
 *      lineCommand -- the command to be processed
 * Function
 *      Get rid of end of line, and information after #.
 *      Get rid of before and after tabs, replace between tabs with a space.
 *      Get rid of before and after spaces, retain only one space in words.
 *      If no other information exists, the lineCommand turns to a NULL string.
 * Returns
 *      0 -- successful
 */
extern int CommandLineProcessor(char *lineCommand);
/*
 * Fatal error control
 *
 * Parameters
 *      statement -- the information to show
 * Function
 *      Print information and then exit. Once the process exits, the operating
 *      system is able to free all dynamically allocated memory associated with
 *      the process.
 */
extern void FatalError(const char *statement);
/*
 * Show information to terminal
 *
 * Parameters
 *      statement -- the information to show. If statement is "Session End",
 *      it prints a line asterisks.
 * Function
 *      Print information to standard out.
 */
extern int ShowInformation(const char *statement);
/*
 * Assign Storage
 *
 * Parameters
 *      idxMax -- the maximum number of elements
 *      dataType -- the data type of elements, can be "int", "double", "char",
 *      "Real", "float"
 * Function
 *      Use malloc to assign a linear array of storage with specified data type.
 * Returns
 *      The head address of the assigned storage.
 * Notice
 *      malloc does not initialize the storage; this means that the assigned
 *      memory may contain random or unexpected values!
 * Returns
 *      0 -- successful
 */
extern void *AssignStorage(const int idxMax, const char *dataType);
/*
 * Retrieve storage from a pointer.
 *
 * Parameter
 *      pointer -- the pointer that contains the storage address
 * Function
 *      Use free to free the storage space of the pointer.
 * Notice
 *      Don't free pointer of storage that not allocated by dynamic allocation.
 *      The original pointer becomes to be a wild pointer after being freed, be
 *      aware of this situation. It's a better practice to set pointer back to NULL
 *      after calling free.
 * Returns
 *      0 -- successful
 */
extern int RetrieveStorage(void *pointer);
#endif
/* a good practice: end file with a newline */

