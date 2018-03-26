/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 *
 *                             References
 *
 * - C wikibook http://en.wikibooks.org/wiki/C_Programming
 * - Practical C Programming by Steve Oualline
 * - Expert C Programming: Deep C Secrets by Peter van der Linden
 * - The Practice of Programming by Brian W. Kernighan and Rob Pike
 * - Effective C++ by Scott Meyers
 * - More Effective C++ by Scott Meyers
 * - Writing Scientific Software: A Guide to Good Style, by Suely Oliveira
 * - C Header File Guidelines by David Kieras
 * - Tips for Optimizing C/C++ Code by Clemson
 * - How to loop through multidimensional arrays quickly by Nadeau software
 * - How expensive is an operation on a CPU by Vincent Hindriksen
 * - Performance Tuning with the RESTRICT Keyword by David H Bartley
 *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 * - Header files should be once-only headers. A standard way to prevent 
 *   including a header file more than once is to enclose the entire
 *   contents of the file in a conditional.
 * - Do not start the guard symbol with an underscore. Leading underscore 
 *   names are reserved for internal use by the C implementation.
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
 * - Avoid side effects, such as do not use assignment statements in if 
 *   condition, should use ++ and -- on lines by themselves.
 * - Always use the prefix version of ++ and -- (++x, --x) instead of the
 *   postfix version (x++, x--).
 * - Remember that there are only two precedence levels to remember in 
 *   C: multiplication and division come before addition and subtraction.
 *   Everything else should be in parentheses. 
 * - Avoid complex logic like multiply nested ifs. Consider splitting your 
 *   code into multiple procedures to decrease the level of complexity.
 * - Variables should be declared as locally as possible:
 *    * declare non-constant variables that are used throughout the function
 *      at the top.
 *    * declare constant variables when their values can be determined.
 *    * declare variables that are used in only a local scope of the function
 *      at the point where they are needed, e.g., loop counts.
 * - Use of const whenever possible, using const as much as possible is
 *   compiler-enforced protection from unintended writes to data
 *   that should be read-only.
 *    * Declare variables that should not be changed after initialization.
 *          const double pi = 3.14159265359;
 *    * Two ways of declaring a const pointer: 
 *      a) the target address which the pointer points to is fixed, but 
 *      the content in the address can be changed:
 *      double *const pointer; (pointer itself is const, can't be changed)
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
 * - Two key steps: Algorithm development. Devise a clear step-by-step solution 
 *   strategy for the problem based on formal logic. Algorithm implementation. 
 *   Translate the algorithm into source code. To be a good translator, the
 *   information flow and the structure of the code shall resemble those in the
 *   solution strategy to minimize mistakes and to improve readability.
 *
 * - Elements of a Program: In computer programming, you need two things: 
 *   data (variables) and instructions (code or functions). Variables are
 *   the basic building blocks of a program. Instructions tell the computer
 *   what to do with the variables:
 *
 *                           Instructions
 *             Input data ------------------> Output data
 *
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
 * the program into modules, then divide the modules into submodules, then 
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
 * hiding and good interface design are key to making good modules.
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
 *   be declared static. Private external variables (variables declared outside
 *   of all functions) that are not used outside the module should be static.
 * - Obviously, the prototypes for static functions should not be put in the
 *   module's header file.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                            C Puzzles
 *
 * - Definition VS Declaration: 
 *   Declaration of a variable/function declares that the variable/function 
 *   exists somewhere in the program. No memory will be allocated by a
 *   declaration. Declaration can occur multiple times and describes the 
 *   type of an object; It is used to refer to objects defined, since a 
 *   declaration of a variable/function is always needed to be given before 
 *   anything that wants to access them. Declarations subject to scope rule.
 *   Definition of a variable/function, apart from the role of declaration, 
 *   it also allocates memory for that variable/function. It is used to create
 *   new objects, example: int my_array[100]; A variable/function can only be
 *   defined once within its scope.
 * - extern keyword: (an object is a variable or a function)
 *   In the C language, an external (global) object is an object defined
 *   outside any function block (external to all functions). A local object
 *   is a object defined inside a function block. External objects are
 *   allocated and initialized when the program starts, and the memory is only
 *   released when the program ends. External objects are globally accessible
 *   and remain in existence permanently, therefore, they can be used to 
 *   communicate data globally between functions. A declaration of an object
 *   must be specified before any thing to access it. An external object 
 *   is directly accessible to all the functions in the same module where the
 *   external object is defined, since definition also serves as declaration.
 *   For functions in other module files to access the object, a declaration
 *   is needed to refer to the object, which is done by the extern keyword.
 *   The extern keyword means "declare without defining". It is a way to
 *   explicitly declare an object without a definition.
 *   Since all external objects are globally accessible in default, their 
 *   scope is global and hence they must be defined exactly once in one of
 *   the modules of the program. For modules that do not define the
 *   external object to access it, a declaration is needed to connect the 
 *   occurrences of the object, which is greatly facilitated by header files
 *   to ensure that all the declarations used are consistent with each other
 *   and with the definition. To simplify the declaration of external objects 
 *   for modules, the usual practice is to collect extern declarations of 
 *   objects in a separate file called a header, which is included by 
 *   #include at the front of each source file. The normal methodology is 
 *   for allocation and actual definitions to go into .c files, and mere 
 *   declarations and prototypes that do not allocate but just describe the
 *   types of objects for others to refer to and access should go to .h files.
 *   The reliable way to declare and define global objects is to use a header
 *   file to contain an extern declaration of the object. The header is 
 *   included by the one source file that defines the object to ensure that
 *   the definition and the declaration are consistent and by all the source
 *   files that reference the object. For each program, one and only one 
 *   source file defines the object. One and only one header file should 
 *   declare the object.
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
 *
 ****************************************************************************/
/****************************************************************************
 * - C uses void for two purposes: In a function declaration, void indicates
 *   that the function returns no value or takes no arguments. In a pointer 
 *   declaration, void defines a generic pointer.
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
 *   it can be changed, which prevents the desired optimizations.
 * - The strict aliasing and restrict keyword are introduced to the C99 standard
 *   to address the aliasing problem. They are declarations of intent given by 
 *   the programmer to the compiler. Strict aliasing means that two objects of
 *   incompatible types cannot refer to the same location in memory, which is
 *   enabled by passing -fstrict-aliasing flag to the compiler. Be sure that all
 *   code can safely run with this rule enabled. The restrict keyword says that
 *   for the lifetime of the pointer, only the pointer itself or a value directly
 *   derived from it (such as pointer + 1) will be used to access the object to
 *   which it points. This limits the effects of pointer aliasing, that is,
 *   each memory block pointed by a restrict pointer is only accessed by the
 *   current pointer. Since the strict aliasing rules prohibit aliasing among
 *   incompatible types, and different restrict pointers of compatible types
 *   always point to different locations, updating one pointer will not 
 *   affect the other pointers, aiding better optimizations. Whether aliasing
 *   does or does not occur is the responsibility of the programmer.
 * - Begin using the restrict keyword immediately. Retrofit old code as soon as
 *   possible. Only use restricted leaf pointers. Use of parent pointers may
 *   break the restrict contract.
 * - Keep loads and stores separated from calculations. This results in better
 *   scheduling in compilers, and makes the relationship between the output 
 *   assembly and the original source clearer.
 *
 ****************************************************************************/
/****************************************************************************
 *
 *                 Best Practices for Scientific Programming
 *
 * Use revision control system
 *
 * - Extremely useful for comparing, recovering, maintenance, etc.
 * - Available options: CVS, Subversion, Github.
 * 
 * Compile with make for automatic build procedures
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
 * Use a robust and exhaustive test suite
 *
 * - Verify the functionality of the software whenever modified.
 * - Should be coupled to the build infrastructure with every release.
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
 *   Optimizing compilers know more about the target architecture than
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
 *   slogan that 95% of the time is spent in 5% of the code. The task is
 *   to first identify the 5% code that really matters. If there is a
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
 * - What Every Computer Scientist Should Know About Floating-Point Arithmetic.
 *   assign a value which is beyond the maximum value of that data type
 *   compiler will assign +INF if number is positive and -INF if number is
 *   negative. If assign a value witch is less than minimum value of
 *   that data type then complier will assign a garbage value or zero.
 * - Floating-point exception handling.
 *   When code generates an overflow, underflow, or divide-by-zero error, the 
 *   result will simply be an infinite or not-a-number value. If undefined 
 *   values are used in other operations, new undefined values are generated. 
 *   Then the program gives inf or NaN as a result. One can use the floating
 *   point exception facilities provided by C in fenv.h to determine and track
 *   a floating-point exceptional condition when it first occurred, which can 
 *   provide useful information for debugging.
 *
 ****************************************************************************/
/****************************************************************************
 * 
 *                             C Textbooks
 *
 * - Reference Style
 *   The C Programming Language by Brian W. Kernighan and Dennis M. Ritchie
 * - Beginner
 *   C wikibook http://en.wikibooks.org/wiki/C_Programming
 *   Practical C Programming by Steve Oualline
 *   C Programming: A Modern Approach by K. N. King
 * - Intermediate
 *   Algorithms in C by Robert Sedgewick
 * - Above Intermediate
 *   Expert C Programming: Deep C Secrets by Peter van der Linden
 * - Software engineering
 *   The Practice of Programming by Brian W. Kernighan and Rob Pike
 *   Advanced Programming in the UNIX Environment
 * - Writing Scientific Software: A Guide to Good Style, by Suely Oliveira
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
 *   dereferencing a pointer.
 ****************************************************************************/
/****************************************************************************
 *
 *             Coding Principles for Numerical Solutions of PDEs
 * 
 * - The structure of the code should strictly resemble the structure of 
 *   the mathematical formulation. For example, the data flow in the code
 *   should follow the data flow in the mathematical formulation. 
 *   The function hierarchy in the code should follow the function hierarchy
 *   in the mathematical formulation; and the function definition and content
 *   in the code should be the same as those in the mathematical formulation. 
 * - The code verification should be thorough. Start from simple problems that
 *   numerical methods could solve exactly, such as uniform or linearly
 *   distributed smooth initial data. In these simple problems, the truncation
 *   error would be zero since the high order derivatives are all zero. 
 *   Therefore, numerical methods would solve these problems exactly. 
 * - Convergence test with convergence rate evaluated should be conducted.
 *   The employed grids should be successively refined with grid number
 *   successively doubled. The chosen problem should be smooth enough to
 *   avoid the affection from discontinuities.
 *
 ****************************************************************************/
/*
 * Define some global integer constants for array bounds, identifiers, etc.
 * enum statement is used instead of macros for handling these magic numbers.
 * An enumerator with = defines its enumeration constant as the value of the
 * constant expression. Using enum member as an array size will not result
 * a variable length array since enumeration constants are constant
 * expressions, and if the size is an integer constant expression and the
 * element type has a known constant size, the array type is not a variable
 * length array type. 
 */
typedef enum {
    /* dimensions related to space */
    DIMS = 3, /* space dimension */
    X = 0,
    Y = 1,
    Z = 2,
    COLLAPSEN = 0, /* dimension collapsed tag */
    COLLAPSEX = 1,
    COLLAPSEY = 2,
    COLLAPSEZ = 3,
    COLLAPSEXY = 5,
    COLLAPSEXZ = 7,
    COLLAPSEYZ = 8, 
    COLLAPSEXYZ = 17, 
    /* dimensions related to field variables */
    DIMU = 5, /* conservative vector: rho, rho_u, rho_v, rho_w, rho_eT */
    DIMUo = 6, /* primitive vector: rho, u, v, w, [p, hT, h], [T, c] */
    DIMT = 3, /* number of time levels to store field data */
    TO = 0, /* the time level for current */
    TN = 1, /* the time level for intermediate */
    TM = 2, /* the time level for intermediate */
    /* parameters related to numerical model */
    PATHN = 30, /* neighbour searching path */
    PATHSEP = 4, /* layer separator in neighbour searching path: pathN, l1N, l2N, l3N */
    NONE = -1, /* invalid flag */
    WENOTHREE = 0, /* 3th order weno */
    WENOFIVE = 1, /* 5th order weno */
    /* parameters related to domain partitions */
    NPART = 13, /* inner region, [west, east, south, north, front, back] x [Boundary, Ghost] */
    NPARTWRITE = 1, /* number of partitions to write data out */
    PIN = 0,
    PWB = 1, 
    PEB = 2, 
    PSB = 3,
    PNB = 4,
    PFB = 5,
    PBB = 6,
    PWG = 7,
    PEG = 8,
    PSG = 9,
    PNG = 10,
    PFG = 11,
    PBG = 12,
    LIMIT = 2, /* number of limits */
    MIN = 0,
    MAX = 1,
    /* parameters related to domain boundary conditions */
    NBC = 7, /* Interior, [west, east, south, north, front, back] x [Boundary] */
    INFLOW = 0, /* boundary condition identifier */
    OUTFLOW = 1,
    SLIPWALL = 2,
    NOSLIPWALL = 3,
    PERIODIC = 4,
    ENTRYBC = 6, /* rho, u, v, w, p, T */
    VARBC = 5, /* rho, u, v, w, p */
    /* parameters related to global and regional initialization */
    NIC = 10, /* maximum number of initializer to support */
    ICGLOBAL = 0, /* global initializer */
    ICPLANE = 1, /* plane initializer */
    ICSPHERE = 2, /* sphere initializer */
    ICBOX = 3, /* box initializer */
    ICCYLINDER = 4, /* cylinder initializer */
    ENTRYIC = 12, /* x1, y1, z1, [x2, Nx], [y2, Ny], [z2, Nz], r, primitive variables */
    VARIC = 5, /* primitive variables: rho, u, v, w, p */
    /* parameters related to geometry */
    DIMTK = 2, /* number of time levels to store kinematic data */
    POLYN = 3, /* polygon facet type */
    EVF = 4, /* edge-vertex-face type */
} Constants;
/*
 * Define some universe data type for portability and maintenance.
 */
typedef double Real; /* real data */
typedef char String[400]; /* string data */
typedef int IntVec[DIMS]; /* integer type vector */
typedef Real RealVec[DIMS]; /* real type vector */
/*
 * Define structures for packing compound data
 *
 * To reduce padding required for data alignment of structures, arranging
 * members of a structure in increasing order by size is desirable.
 * In addition, fields shall be sorted by their frequency or by memory 
 * access pattern. Put frequently accessed elements to small offsets,
 * and if two elements are used at the same time, put them closely to reduce
 * cache misses per structure.
 */
/*
 * Field variables of computational node
 *
 * Using high order pointers for multidimensional arrays wastes space and
 * the malloc calls are expensive, and it is also very time-consuming
 * for nested loops because of causing lots of cache misses.
 *
 * Maintaining a multidimensional array within a single linear array is a
 * common performance technique. High-performance code instead implements a
 * multidimensional array as a single linear array with hand-authored array
 * indexing math to keep track of what values are where:
 *
 * data[kMax][jMax][iMax] mapped to --> data[idxMax]
 * data[k][j][i] mapped to --> data[idx]
 * int kMax, jMax, iMax, k, j, i;
 * int idxMax, idx; use long type if needed
 * idxMax = kMax * jMax * iMax;
 * idx = (k * jMax + j) * iMax + i;
 *
 * Since the array is a single large chunk of memory, sweeping through it from
 * start-to-finish creates a regular access pattern that processor prefetchers
 * easily recognize, which enables them to load caches in the background. The
 * result is fewer cache misses and much better performance.
 *
 * This code use nested loops with linear array and index math instead nested
 * loops with multidimensional arrays. Because we need to store a sequence of
 * different field variables at each (k,j,i), a node data structure is defined
 * to pack data for each node and to achieve high locality of data management.
 *
 */
typedef struct {
    int gid; /* geometry identifier */
    int fid; /* closest face identifier */
    int lid; /* interfacial layer identifier */
    int gst; /* ghost layer identifier */
    Real U[DIMT][DIMU]; /* field data at each time level */
} Node;
/*
 * Domain discretization and partition structure
 */
typedef struct {
    IntVec m; /* mesh number of spatial dimensions */
    IntVec n; /* node number of spatial dimensions */
    int ng; /* number of ghost node layers of global domain */
    int gl; /* number of ghost node layers required for numerical scheme */
    int collapse; /* space collapse flag */
    RealVec d; /* mesh size of spatial dimensions */
    RealVec dd; /* reciprocal of mesh sizes */
    Real tinyL; /* smallest length scale established on grid size */
    int ns[NPART][DIMS][LIMIT]; /* decomposition node range for each partition */
    int np[DIMS][DIMS][LIMIT]; /* computational node range with dimension priority */
    int path[PATHN][DIMS]; /* neighbour searching path */
    int pathSep[PATHSEP]; /* layer separator in neighbour searching path */
    int N[NBC][DIMS]; /* outward surface normal of domain boundary */
    int typeBC[NBC]; /* BC type recorder */
    int countIC; /* flow initializer count */
    int typeIC[NIC]; /* flow initializer type recorder */
    Real valueBC[NBC][ENTRYBC]; /* field values of each boundary */
    Real valueIC[NIC][ENTRYIC]; /* field values of each initializer */
    Real domain[DIMS][LIMIT]; /* coordinates define the space domain */
} Partition;
/*
 * Facet structure
 */
typedef struct {
    RealVec N; /* normal vector */
    RealVec v0; /* vertex */
    RealVec v1; /* vertex */
    RealVec v2; /* vertex */
} Facet;
/*
 * Polyhedron structure
 */
typedef struct {
    int faceN; /* number of faces. 0 for analytical sphere */
    int edgeN; /* number of edges */
    int vertN; /* number of vertices */
    int state; /* dynamic motion indicator */
    int mid; /* material type */
    Real r; /* bounding sphere */
    RealVec O; /* centroid */
    Real I[DIMS][DIMS]; /* inertia matrix */
    Real V[DIMTK][DIMS]; /* translational velocity */
    Real W[DIMTK][DIMS]; /* rotational velocity */
    Real at[DIMTK][DIMS]; /* translational acceleration */
    RealVec g; /* gravitational acceleration */
    Real ar[DIMTK][DIMS]; /* rotational acceleration */
    RealVec Fp; /* pressure force */
    RealVec Fv; /* viscous force */
    RealVec Tt; /* total torque */
    Real to; /* time to end power */
    Real rho; /* density */
    Real T; /* wall temperature */
    Real cf; /* roughness */
    Real area; /* area */
    Real volume; /* volume */
    Real box[DIMS][LIMIT]; /* a bounding box of the polyhedron */
    int (*restrict f)[POLYN]; /* face-vertex list */
    Real (*restrict Nf)[DIMS]; /* face normal */
    int (*restrict e)[EVF]; /* edge-vertex-face list */
    Real (*restrict Ne)[DIMS]; /* edge normal */
    Real (*restrict v)[DIMS]; /* vertex list */
    Real (*restrict Nv)[DIMS]; /* vertex normal */
    Facet *facet; /* facet data */
} Polyhedron;
/*
 * Collision list
 */
typedef struct {
    int gid; /* geometry identifier */
    IntVec N; /* line of impact */
} Collision;
/*
 * Geometry Entities
 */
typedef struct {
    int totN; /* total number of geometries */
    int sphN; /* number of analytical spheres */
    int stlN; /* number of triangulated polyhedrons */
    int colN; /* colliding list pointer and count */
    Polyhedron *poly; /* geometry list */
    Collision *col; /* collision list */
} Geometry;
/*
 * Material properties
 */
typedef struct {
    Real eos; /* equation of state */
} Material;
/*
 * Space domain parameters
 */
typedef struct {
    Node *node; /* field data */
    Geometry geo; /* geometry in space */
    Partition part; /* domain discretization and partition information */
} Space;
/*
 * Time domain parameters
 */
typedef struct {
    int restart; /* restart tag */
    int stepN; /* total number of steps */
    int stepC; /* step number count */
    int writeN; /* field data writing frequency */
    int writeC; /* field data writing count */
    int dataStreamer; /* types of data streamer */
    int pointWriteN; /* point probe writing frequency */
    int lineWriteN; /* line probe writing frequency */
    int curveWriteN; /* body-conformal probe writing frequency */
    int forceWriteN; /* surface force writing frequency */
    int pointProbeN; /* total number of point probes */
    int lineProbeN; /* total number of line probes */
    int curveProbeN; /* body-conformal probe */
    int forceProbeN; /* surface force probe */
    Real end; /* termination time */
    Real now; /* current time recorder */
    Real numCFL; /* CFL number */
    Real (*restrict pp)[DIMS]; /* point probes */
    Real (*restrict lp)[7]; /* line probes */
} Time;
/*
 * Model properties and physics parameters
 */
typedef struct {
    int tScheme; /* temporal discretization scheme */
    int sScheme; /* spatial discretization scheme */
    int multidim; /* multidimensional space method */
    int jacobMean; /* average method for local Jacobian linearization */
    int fluxSplit; /* flux vector splitting method */
    int fsi; /* material interaction trigger */
    int ibmLayer; /* number of interfacial layers using flow reconstruction */
    int mid; /* material identifier */
    int gState; /* gravity state */
    int sState; /* source state */
    Real refMa; /* reference Mach number */
    Real refMu; /* reference dynamic viscosity */
    Real gamma; /* heat capacity ratio */
    Real gasR; /* specific gas constant */
    Real cv; /* specific heat capacity at constant volume */
    Real refL; /* characteristic length */
    Real refRho; /* characteristic density */
    Real refV;  /*characteristic velocity */
    Real refT; /* characteristic temperature */
    RealVec g; /* gravity vector */
    Material mat; /* material database */
} Model;
/*
 * Program command line arguments and overall control
 */
typedef struct {
    char runMode; /* mode: [i] interact, [s] serial, [t] threaded, [m] mpi, [g] gpu */
    int procN; /* number of processors */
} Control;
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Command line processor
 *
 * Function
 *      Get rid of end of line, and information after #.
 *      Get rid of before and after tabs, replace tabs with a space.
 *      Get rid of before and after spaces, retain only one space in words.
 *      If no other information exists, the lineCommand turns to a NULL string.
 */
extern int CommandLineProcessor(char *lineCommand);
/*
 * Fatal error control
 *
 * Function
 *      Print information and then exit. Once the process exits, the operating
 *      system is able to free all dynamically allocated memory associated with
 *      the process.
 */
extern void FatalError(const char *statement);
/*
 * Show information to terminal
 *
 * Function
 *      Print information to standard out. Statement is the information to show. 
 *      If statement is "Session End", it prints a line asterisks.
 */
extern int ShowInformation(const char *statement);
/*
 * Assign Storage
 *
 * Function
 *      Use malloc to assign a linear array of storage. Returns the head address
 *      of the assigned storage. Since malloc does not initialize the storage, 
 *      a call of memset is used to initialize the assigned memory to zero.
 */
extern void *AssignStorage(size_t size);
/*
 * Retrieve storage
 *
 * Function
 *      Use free to free the storage space of the pointer.
 * Notice
 *      Don't free pointer of storage that not allocated by dynamic allocation.
 *      The original pointer becomes to be a wild pointer after being freed, be
 *      aware of this situation. It's a better practice to set pointer back to 
 *      NULL after calling free.
 */
extern int RetrieveStorage(void *pointer);
/*
 * Auxiliary Functions for File Reading
 *
 * Function
 *      Read in lines from the current line until a line matches the lineString.
 *      The file pointer points to the next line of the matched line.
 */
extern int ReadInLine(FILE *filePointer, const char *lineString);
/*
 * Auxiliary Functions for File Writing
 *
 * Function
 *      Search down the file from beginning until a line matches the lineString.
 *      The file pointer points to the matched line.
 */
extern int WriteToLine(FILE *filePointer, const char *lineString);
/*
 * Standard Stream Functions with Checked Return Values
 */
extern void Fgets(char *str, int num, FILE *stream);
extern void Fread(void *ptr, size_t size, size_t count, FILE *stream);
extern void VerifyReadConversion(const int num, const int expect);
#endif
/* a good practice: end file with a newline */

