/****************************************************************************
 * Header File                                                              *
 * Last-modified: 19 Jan 2015 07:57:14 AM
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
 *   a named file in the directory containing the current file, 
 *   such as #include "include/foo.h".
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
 * - Avoid complex logic like multiply nested ifs. Consider splitting your 
 *   code into multiple procedures, to decrease the level of complexity.
 * - Variables should be declared as locally as possible.
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
 * More Details
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
 *   pass its adress out if access is needed outside its local scope.
 * - In C single quotes identify a single character, while double quotes 
 *   create a string literal. 'a' is a single a character literal, while "a" 
 *   is a string literal containing an 'a' and a null terminator (that is a
 *   2 char array). Note that in C, the type of a character literal is int,
 *   and not char, that is sizeof 'a' is 4 in an architecture where ints 
 *   are 32bit (and CHAR_BIT is 8), while sizeof(char) is 1 everywhere.
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
 *                     Code Performance and Optimization
 *
 * - Rule of thumb, let the compiler do the job. Code should compiled with 
 *   compiler optimizations. ICC or GCC are the best.
 * - keep your code general instead of obsessively optimizing your code with 
 *   awkward code structure. Hand optimizations can create odd looking code
 *   that is harder for a compiler to match up to an optimization template.
 *   It is often better to resist the temptation to optimize the code. 
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
 * - Use struct to pass a bunch of data at a time, it's simple and elegant.
 *   However, make sure about these:
 * - Try to use regular variables and pointers as struct members, don't use 
 *   extremely large arrays which may decrease performance.
 * - Always pass structures by reference. That is, use pointers to structures
 *   as function arguments even when nothing in the struct will be modified in
 *   the function. Should never do value passing to avoid copying the complete
 *   contents of the structure onto the stack.
 * - Make sure assigning a valid memory location to the pointer before
 *   dereferencing a pointer!
 ****************************************************************************/
/*
 * Field variables of flow
 *
 * Conservative variables are vectors with five elements(rho, rho_u, rho_v,
 * rho_w, E), while each element is a three dimensional array in 3D space.
 * Thus, the conservative variables need to be presented as a 4 dimensional
 * array in 3D flow. Besides, primitive variable Uo has one more element 
 * than U, since Uo stands for primitive variables: rho, u, v, w, p, T.
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
 * loops with multi-dimensional arrays, take an three dimensional array as an
 * example:
 *
 * U[width][height][depth] ~ U[kMax][jMax][iMax];
 * int k, j, i;
 * int idx; use long type if needed
 * idx = (k * jMax * iMax) + (j * iMax) + i; (three multipliers)
 * idx = (k * jMax + j) * iMax + i; (two multipliers)
 * value = U[idx];
 */
typedef struct {
    double *U; /* conservative flow variables at n+1*/
    double *Un; /* conservative flow variables at time n */
    double *Um; /* conservative flow variables at time n-1 */
    double *Uo; /* primitive flow variables */
}Field;
/*
 * Flux variables
 */
typedef struct {
    double *Fx; /* non-viscous flux vector at x direction */
    double *Fy; /* non-viscous flux vector at y direction */
    double *Fz; /* non-viscous flux vector at z direction */
    double *Gx; /* viscous flux vector at x direction */
    double *Gy; /* viscous flux vector at y direction */
    double *Gz; /* viscous flux vector at z direction */
    double *eigenValue; /* eigenvalue */
    double *leftMatrix; /* left matrix */
    double *rightMatrix; /* right matrix */
}Flux;
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
    double dx; /* first use as x length, then update to mesh size in x */
    double dy; /* first use as y length, then update to mesh size in y */
    double dz; /* first use as z length, then update to mesh size in z */
    int *ghostFlag; /* node type integer flag: normal, ghost, solid, etc. */
    int *geoID; /* store the ID of geometry object for each ghost cell */
}Space;
/*
 * Particle Entities
 */
typedef struct {
    int totalN; /* total number of particles */
    double *x; /* x coordinates of the particle center */
    double *y; /* y coordinates of the particle center */
    double *z; /* z coordinates of the particle center */
    double *r; /* radius of the particle */
    double *u; /* x component of the velocity of particles */
    double *v; /* y component of the velocity of particles */
    double *w; /* z component of the velocity of particles */
}Particle;
/*
 * Time domain parameters
 */
typedef struct {
    int restart; /* restart flag */
    double totalTime; /* total evolution time */
    double currentTime; /* current time recorder */
    double dt; /* time step size */
    double numCFL; /* CFL number */
    int totalStep; /* total number of steps */
    int stepCount; /* step number count */
    int totalOutputTimes; /* total times of exporting computed data */
    int outputCount; /* exporting data count */
}Time;
/*
 * Fluid properties
 */
typedef struct {
    double density; /* fluid density */
    double nu; /* kinematic viscosity */
    double alpha; /* thermal diffusivity */
}Fluid;
/*
 * Flow properties and physics parameters
 */
typedef struct {
    double mu; /* generalized normalized dynamic viscosity */
    double heatK; /* generalized normalized thermal conductivity */
    double numMa; /* Mach number */
    double numRe; /* Reynolds number */
    double numPr; /* Prandtl number */
    double gamma; /* heat capacity ratio */
    double gasR; /* the gas constant */
    double cv; /* specific heat capacity at constant volume */
}Flow;
/*
 * Characteristic values for normalization
 */
typedef struct {
    double length; /* characteristic length */
    double density; /* characteristic density */
    double velocity;  /*characteristic velocity */
    double temperature; /* characteristic temperature */
}Reference;
/*
 * Domain partition structure
 */
typedef struct {
    int totalN; /* total number of domain partitions */
    int *idxHead; /* store the malloc address of indices */
    int *kSub; /* domain decomposition control */
    int *kSup; /* domain decomposition control */
    int *jSub; /* domain decomposition control */
    int *jSup; /* domain decomposition control */
    int *iSub; /* domain decomposition control */
    int *iSup; /* domain decomposition control */
    char *nameHead; /* store the malloc address of names */
    int nameLength; /* length of names */
}Partition;
/****************************************************************************
 * Function declaration
 ****************************************************************************/
int CommandLineProcessor(char *);
void FatalError(const char *);
int ShowInformation(const char *);
void *AssignStorage(const int, const char *);
int RetrieveStorage(void *);
#endif
/* a good practice: end file with a newline */

