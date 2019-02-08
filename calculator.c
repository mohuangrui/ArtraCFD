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
 * Required Header Files
 ****************************************************************************/
#include "calculator.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Local Data Structure
 ****************************************************************************/
typedef enum {
    DIMOP = 19, /* number of operator types */
    NOPRD = 50, /* number of operands */
    NOPRT = 50, /* number of operators */
    LENOPRT = 5, /* longest length of operator symbol */
} CalcConst;
typedef struct {
    Real *const bom; /* pointer to stack bottom */
    Real *top; /* pointer to stack top */
    const int size; /* stack size */
    Real space[NOPRD]; /* stack space */
} Operand;
typedef struct {
    char *const bom; /* pointer to stack bottom */
    char *top; /* pointer to stack top */
    const int size; /* stack size */
    char space[NOPRT]; /* stack space */
    const char priority[DIMOP][DIMOP]; /* store the priority of operators */
} Operator;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int SolveExpression(CalcVar *, const char *, Operand *, Operator *);
static int TranslateToMath(char *);
static char QueryPriority(const Operator *, const char, const char);
static int QueryIndex(const char);
static int PushOperator(Operator *, const char);
static int PushOperand(Operand *, const Real);
static int PopOperator(Operator *, char *);
static int PopOperand(Operand *, Real *);
static char GetTopOperator(const Operator *);
static Real GetTopOperand(const Operand *);
static int IsOperator(const char);
static int IsPureUnaryOperator(const char);
static int IsPureBinaryOperator(const char);
static int IsDualOperator(const char);
static int IsDualOperatorActAsUnary(const char *);
static int IsLeftParentheses(const char);
static int IsRightParentheses(const char);
static int IsConstant(const char);
static int IsEndTag(const char);
static int IsDigit(const char);
static int IsDot(const char);
static Real ExtractFirstFloat(char **);
static Real ExtractConstant(const CalcVar *, char **);
static int DoUnary(const char, const Real, Real *);
static int DoBinary(const Real, const char, const Real, Real *);
static void ShowCalcManual(void);
static int SetVariable(CalcVar *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int RunCalculator(void)
{
    CalcVar theVar = {
        .t = 0.0,
        .x = 0.0,
        .y = 0.0,
        .z = 0.0,
        .ans = 0.0,
        .pi = PI,
    };
    ShowInfo("Session");
    ShowInfo("*                   Expression Calculator                  *\n");
    ShowInfo("*                     <By Huangrui Mo>                     *\n");
    ShowInfo("*    Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>     *\n");
    ShowInfo("Session");
    ShowInfo("Enter 'help' for more information\n");
    ShowInfo("Session");
    String str = {'\0'}; /* store the input information */
    while (1) {
        ShowInfo("\nCalc << ");
        ParseCommand(fgets(str, sizeof str, stdin));
        ShowInfo("\n");
        if (0 == strncmp(str, "help", sizeof str)) {
            ShowInfo("Options:\n");
            ShowInfo("[help]       show this information\n");
            ShowInfo("[set]        set variable values\n");
            ShowInfo("[math expr]  compute math expression\n");
            ShowInfo("[manual]     show user manual\n");
            ShowInfo("[exit]       return\n");
            continue;
        }
        if (0 == strncmp(str, "set", sizeof str)) {
            SetVariable(&theVar);
            continue;
        }
        if (0 == strncmp(str, "manual", sizeof str)) {
            ShowCalcManual();
            continue;
        }
        if ('\0' == str[0]) {
            continue;
        }
        if (0 == strncmp(str, "exit", sizeof str)) {
            ShowInfo("Session");
            return 0;
        }
        /* if non of above is true, then compute the expression */
        ComputeExpression(&theVar, str);
        ShowInfo("ans = %.6g\n", theVar.ans);
    }
    return 0;
}
Real ComputeExpression(CalcVar *var, const char *str)
{
    Operand theOprd = {
        .bom = theOprd.space,
        .top = theOprd.space,
        .size = NOPRD - 1,
        .space = {0.0}
    };
    Operator theOprt = {
        .bom = theOprt.space,
        .top = theOprt.space,
        .size = NOPRT - 1,
        .space = {'\0'},
        .priority = { /* + - * / ^ exp ln lg abs sin cos tan ( ) [ ] { } \0 */
            {'>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '>', '<', '>', '<', '>', '>'},
            {'<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '=', '<', 'w', '<', 'w', 'w'},
            {'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w'},
            {'<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', 'w', '<', '=', '<', 'w', 'w'},
            {'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w'},
            {'<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', 'w', '<', 'w', '<', '=', 'w'},
            {'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w'},
            {'<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<', 'w', '<', 'w', '<', 'w', '='}
        }
    };
    SolveExpression(var, str, &theOprd, &theOprt);
    return var->ans;
}
/*
 * The flow control of this program is important, thus, every function
 * call which may result an important error will be monitored.
 * These functional functions return 0 means true, 1 means false.
 */
static int SolveExpression(CalcVar *var, const char *str, Operand *oprd, Operator *oprt)
{
    /* use a new space to avoid modify the original expression */
    String space = {'\0'};
    /*
     * Use as a flag to check whether '+' and '-' are unary operator based on
     * the assumption that if "+" "-" appears as unary operator then it must
     * after a parentheses except at the beginning of an expression. By putting
     * a "(" at the front of the expression, then this rule is always true.
     */
    space[0] = '(';
    /* target expression begins at the second position of storage address */
    char *expr = space + 1;
    /*
     * Copy the expression to space, because need to access expr[2] and even
     * later position when search and convert the expr to math expression,
     * limit the available space for copy is a safety choice.
     */
    strncpy(expr, str, (int)(sizeof space) - LENOPRT);
    if (0 != TranslateToMath(expr)) {
        return 1;
    }
    char toprt = '\0'; /* store the top operator of operator stack */
    char coprt = '\0'; /* store the current operator */
    Real oprdx = 0.0; /* first top operand in stack */
    Real oprdy = 0.0; /* second top operand in stack */
    Real coprd = 0.0; /* store the current operand */
    /* always initialize and reset the stack status */
    oprt->top = oprt->bom;
    oprd->top = oprd->bom;
    PushOperator(oprt,'\0'); /* push a end flag to the operator stack */
    /* calculation loop */
    while (('\0' != *expr) || ('\0' != GetTopOperator(oprt))) {
        if (IsDigit(*expr)) { /* find a operand */
            coprd = ExtractFirstFloat(&expr);
            if (0 != PushOperand(oprd, coprd)) {
                return 1;
            }
            continue;
        }
        if (IsConstant(*expr)) { /* find a constant */
            coprd = ExtractConstant(var, &expr);
            if (0 != PushOperand(oprd, coprd)) {
                return 1;
            }
            continue;
        }
        /* treat everything left as an operator */
        if (!IsOperator(*expr)) {
            ShowWarning("undefined operator in expression");
            return 1;
        }
        coprt = *expr;
        switch (QueryPriority(oprt, GetTopOperator(oprt), coprt)) {
            case '<': /* push the new high-priority operator into stack */
                if (IsDualOperatorActAsUnary(expr)) {
                    /*
                     * Dual operator (they can both be unary and binary)
                     * '+' '-' show as a unary operator, then push a 0 to
                     * operand stack to make them become binary operator
                     */
                    if (0 != PushOperand(oprd, 0)) {
                        return 1;
                    }
                }
                if (0 != PushOperator(oprt, coprt)) {
                    return 1;
                }
                ++expr;
                break;
            case '=': /* dump the two same-priority operators as they are control operators */
                if (0 != PopOperator(oprt, &coprt)) {
                    return 1;
                }
                ++expr;
                break;
            case '>': /* hold the new low-priority operator and finish the old high-priority operator */
                if (0 != PopOperator(oprt, &toprt)) {
                    return 1;
                }
                if (IsPureUnaryOperator(toprt)) {
                    if (0 != PopOperand(oprd, &oprdx)) {
                        return 1;
                    }
                    if (0 != DoUnary(toprt, oprdx, &coprd)) {
                        return 1;
                    }
                    if (0 != PushOperand(oprd, coprd)) {
                        return 1;
                    }
                } else {
                    if (0 != PopOperand(oprd, &oprdx)) {
                        return 1;
                    }
                    if (0 != PopOperand(oprd, &oprdy)) {
                        return 1;
                    }
                    if (0 != DoBinary(oprdy, toprt, oprdx, &coprd)) {
                        return 1;
                    }
                    if (0 != PushOperand(oprd, coprd)) {
                        return 1;
                    }
                }
                break;
            default:
                ShowWarning("can not match parenthesis");
                return 1;
        }
    }
    /*
     * If the loop successfully exit, then it means the expression
     * and operator stack are all processed, now need to check the
     * operand stack, if there are more than one element left,
     * it means something wrong happened.
     */
    if (1 != (oprd->top - oprd->bom)) {
        ShowWarning("wrong expression");
        return 1;
    }
    /* save the result to answer */
    var->ans = GetTopOperand(oprd);
    return 0;
}
/*
 * Obtain the priority between two operators from the priority matrix
 */
static char QueryPriority(const Operator *oprt, const char oprtx, const char oprty)
{
    const int i = QueryIndex(oprtx);
    const int j = QueryIndex(oprty);
    return oprt->priority[i][j];
}
/*
 * Query the index of a operator in the priority matrix
 */
static int QueryIndex(const char op)
{
    int i = 0;
    switch(op) {
        case '+':
            i = 0; break;
        case '-':
            i = 1; break;
        case '*':
            i = 2; break;
        case '/':
            i = 3; break;
        case '^':
            i = 4; break;
        case 'e':
            i = 5; break;
        case 'n':
            i = 6; break;
        case 'g':
            i = 7; break;
        case 'a':
            i = 8; break;
        case 's':
            i = 9; break;
        case 'c':
            i = 10; break;
        case 't':
            i = 11; break;
        case '(':
            i = 12; break;
        case ')':
            i = 13; break;
        case '[':
            i = 14; break;
        case ']':
            i = 15; break;
        case '{':
            i = 16; break;
        case '}':
            i = 17; break;
        case '\0':
            i = 18; break;
        default:
            ShowWarning("unidentified operator");
            break;
    }
    return i;
}
static int PushOperand(Operand *oprd, const Real op)
{
    if ((oprd->top - oprd->bom) >= oprd->size) {
        ShowWarning("operand stack is overflowing...");
        return 1;
    }
    *oprd->top = op;
    ++(oprd->top);
    return 0;
}
static int PopOperand(Operand *oprd, Real *pop)
{
    if (oprd->top == oprd->bom) {
        ShowWarning("no sufficient operands in expression...");
        return 1;
    }
    --(oprd->top);
    *pop = *oprd->top;
    return 0;
}
static Real GetTopOperand(const Operand *oprd)
{
    return oprd->top[-1];
}
static int PushOperator(Operator *oprt, const char op)
{
    if ((oprt->top - oprt->bom) >= oprt->size) {
        ShowWarning("operator stack is overflowing...");
        return 1;
    }
    *oprt->top = op;
    ++(oprt->top);
    return 0;
}
static int PopOperator(Operator *oprt, char *pop)
{
    if (oprt->top == oprt->bom) {
        ShowWarning("no sufficient operators in expression...");
        return 1;
    }
    --(oprt->top);
    *pop = *oprt->top;
    return 0;
}
static char GetTopOperator(const Operator *oprt)
{
    return oprt->top[-1];
}
/*
 * Translate the input string to a specific format that the program
 * can recognize every operator and operand correctly.
 */
static int TranslateToMath(char *str)
{
    char *scanner = str;
    char *receiver = str;
    while ('\0' != *scanner) {
        switch (*scanner) {
            case 'e':
                if (('x' == scanner[1]) && ('p' == scanner[2])) { /* exp */
                    *receiver = 'e';
                    ++receiver;
                    scanner += 3;
                } else {
                    ShowWarning("unknown operator: e..");
                    return 1;
                }
                break;
            case 'l':
                if (('n' == scanner[1])) {/* ln */
                    *receiver = 'n';
                    ++receiver;
                    scanner += 2;
                } else {
                    if (('g' == scanner[1])) { /* lg */
                        *receiver = 'g';
                        ++receiver;
                        scanner += 2;
                    } else {
                        ShowWarning("unknown operator: l..");
                        return 1;
                    }
                }
                break;
            case 'a':
                if (('b' == scanner[1]) && ('s' == scanner[2])) { /* abs */
                    *receiver = 'a';
                    ++receiver;
                    scanner += 3;
                } else {
                    if (('n' == scanner[1]) && ('s' == scanner[2])) { /* ans */
                        *receiver = 'q';
                        ++receiver;
                        scanner += 3;
                    } else {
                        ShowWarning("unknown operator: a..");
                        return 1;
                    }
                }
                break;
            case 's':
                if (('i' == scanner[1]) && ('n' == scanner[2])) { /* sin */
                    *receiver = 's';
                    ++receiver;
                    scanner += 3;
                } else {
                    ShowWarning("unknown operator: s..");
                    return 1;
                }
                break;
            case 'c':
                if (('o' == scanner[1]) && ('s' == scanner[2])) { /* cos */
                    *receiver = 'c';
                    ++receiver;
                    scanner += 3;
                } else {
                    ShowWarning("unknown operator: c..");
                    return 1;
                }
                break;
            case 't':
                if (('a' == scanner[1]) && ('n' == scanner[2])) { /* tan */
                    *receiver = 't';
                    ++receiver;
                    scanner += 3;
                } else { /* t */
                    *receiver = 'u';
                    ++receiver;
                    ++scanner;
                }
                break;
            case 'p':
                if (('i' == scanner[1])) { /* pi */
                    *receiver = 'p';
                    ++receiver;
                    scanner += 2;
                } else {
                    ShowWarning("unknown operator: p..");
                    return 1;
                }
                break;
            case 'x': /* x */
                *receiver = 'x';
                ++receiver;
                ++scanner;
                break;
            case 'y': /* y */
                *receiver = 'y';
                ++receiver;
                ++scanner;
                break;
            case 'z': /* z */
                *receiver = 'z';
                ++receiver;
                ++scanner;
                break;
            case '.':
                if (IsDigit(scanner[1]) && IsDigit(scanner[-1])) {
                    *receiver = *scanner;
                    ++receiver;
                    ++scanner;
                } else {
                    ShowWarning(". is preceded or followed by non digits");
                    return 1;
                }
                break;
            case ' ':
                ++scanner;  /* ignore space */
                break;
            default:
                if (IsDigit(*scanner) || IsPureBinaryOperator(*scanner) ||
                        IsDualOperator(*scanner) || IsLeftParentheses(*scanner) ||
                        IsRightParentheses(*scanner)) { /* legal input single character */
                    *receiver = *scanner;
                    ++receiver;
                    ++scanner;
                } else {
                    ShowWarning("unknown operator in expression");
                    return 1;
                }
                break;
        }
    }
    *receiver='\0';
    return 0;
}
static int IsOperator(const char op)
{
    return (IsPureUnaryOperator(op) || IsPureBinaryOperator(op) ||
            IsDualOperator(op) || IsLeftParentheses(op) ||
            IsRightParentheses(op) || IsEndTag(op));
}
static int IsPureUnaryOperator(const char op)
{
    return (('e' == op) || ('n' == op) || ('g' == op) ||
            ('a' == op) || ('s' == op) ||
            ('c' == op) || ('t' == op));
}
static int IsDualOperator(const char op)
{
    return (('+' == op) || ('-' == op));
}
static int IsDualOperatorActAsUnary(const char *pop)
{
    return (IsDualOperator(pop[0]) && IsLeftParentheses(pop[-1]));
}
static int IsPureBinaryOperator(const char op)
{
    return (('*' == op) || ('/' == op) || ('^' == op));
}
static int IsLeftParentheses(const char op)
{
    return (('(' == op) || ('[' == op) || ('{' == op));
}
static int IsRightParentheses(const char op)
{
    return ((')' == op) || (']' == op) || ('}' == op));
}
static int IsConstant(const char op)
{
    return (('u' == op) || ('x' == op) || ('y' == op) ||
            ('z' == op) || ('p' == op) || ('q' == op));
}
static int IsEndTag(const char op)
{
    return ('\0' == op);
}
static int IsDigit(const char op)
{
    return (('0' <= op) && ('9' >= op));
}
static int IsDot(const char op)
{
    return ('.' == op);
}
static Real ExtractFirstFloat(char **pstr)
{
    const char *fmtI = ParseFormat("%lg");
    char *str = *pstr;
    Real oprd = 0.0;
    Sscanf(str, 1, fmtI, &oprd);
    while (IsDigit(*str)) {
        ++str;
    }
    if (IsDot(*str)) {
        ++str;
    }
    while (IsDigit(*str)) {
        ++str;
    }
    *pstr = str;
    return oprd;
}
static Real ExtractConstant(const CalcVar *var, char **pstr)
{
    char *str = *pstr;
    Real oprd = 0.0;
    switch (*str) {
        case 'u':
            oprd = var->t;
            ++str;
            break;
        case 'x':
            oprd = var->x;
            ++str;
            break;
        case 'y':
            oprd = var->y;
            ++str;
            break;
        case 'z':
            oprd = var->z;
            ++str;
            break;
        case 'p':
            oprd = var->pi;
            ++str;
            break;
        case 'q':
            oprd = var->ans;
            ++str;
            break;
        default:
            ShowWarning("undefined constant value");
            break;
    }
    *pstr = str;
    return oprd;
}
static int DoUnary(const char toprt, const Real oprdx, Real *pcoprd)
{
    const Real zero = 0.0;
    switch (toprt) {
        case 'e':
            *pcoprd = exp(oprdx);
            break;
        case 'n':
            if (zero >= oprdx) {
                ShowWarning("non-positive argument of ln(x)");
                *pcoprd = zero;
                return 1;
            }
            *pcoprd = log(oprdx);
            break;
        case 'g':
            if (zero >= oprdx) {
                ShowWarning("non-positive argument of lg(x)");
                *pcoprd = zero;
                return 1;
            }
            *pcoprd = log10(oprdx);
            break;
        case 'a':
            *pcoprd = fabs(oprdx);
            break;
        case 's':
            *pcoprd = sin(oprdx);
            break;
        case 'c':
            *pcoprd = cos(oprdx);
            break;
        case 't':
            if (zero == cos(oprdx)) {
                ShowWarning("illegal argument of tangent");
                *pcoprd = zero;
                return 1;
            }
            *pcoprd = sin(oprdx) / cos(oprdx);
            break;
        default:
            *pcoprd = zero;
            return 1;
    }
    return 0;
}
static int DoBinary(const Real oprdy, const char toprt, const Real oprdx, Real *pcoprd)
{
    const Real zero = 0.0;
    switch(toprt)
    {
        case '+':
            *pcoprd = oprdy + oprdx;
            break;
        case '-':
            *pcoprd = oprdy - oprdx;
            break;
        case '*':
            *pcoprd = oprdy * oprdx;
            break;
        case '/':
            if (zero == oprdx) {
                ShowWarning("divide by zero");
                *pcoprd = zero;
                return 1;
            }
            *pcoprd = oprdy / oprdx;
            break;
        case '^':
            *pcoprd = pow(oprdy, oprdx);
            break;
        default:
            *pcoprd = zero;
            return 1;
    }
    return 0;
}
static void ShowCalcManual(void)
{
    ShowInfo("\n            Calculator User Manual\n");
    ShowInfo("Operators:   +, -, *, /, x^y, exp(x), ln(x), lg(x), abs(x), sin(x), cos(x), tan(x)\n");
    ShowInfo("             x, y are numbers or variables or expressions; angle should be radian;\n");
    ShowInfo("Variables:   t, x, y, z, ans, pi;\n");
    ShowInfo("Parenthesis: (), [], {}\n");
    ShowInfo("Example:     1.5*sin(-pi/6)-[cos(pi/3)]^2+ln{exp[5*lg(abs(-100))]} = 9\n");
}
static int SetVariable(CalcVar *var)
{
    const char *fmtI = ParseFormat("%lg");
    ShowInfo("\n t=");
    Sread(stdin, -1, fmtI, &(var->t));
    ShowInfo("\n x=");
    Sread(stdin, -1, fmtI, &(var->x));
    ShowInfo("\n y=");
    Sread(stdin, -1, fmtI, &(var->y));
    ShowInfo("\n z=");
    Sread(stdin, -1, fmtI, &(var->z));
    return 0;
}
/* a good practice: end file with a newline */

