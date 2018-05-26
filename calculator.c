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
/*
 * Operand
 */
typedef struct {
    double *const base; /* pointer to stack bottom */
    double *top; /* pointer to stack top */
    const int stacksize;
} OperandStack;
/*
 * Operator
 */
typedef struct {
    char *const base; /* pointer to stack bottom */
    char *top; /* pointer to stack top */
    const int stacksize;
    const char priority[19][19]; /* store the priority of operators */
} OperatorStack;
/*
 * Parameters
 */
typedef struct {
    const double pi; /* PI */
    double answer; /* store the answer of the latest calculation */
    int radianMode; /* angle in radian mode is default */
    double angleFactor; /* factor transfer to radian mode */
} Parameter;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ComputeExpression(const char *, OperandStack *, OperatorStack *, Parameter *);
static int TranslateCommandToMathExpression(char *);
static char QueryPriority(const OperatorStack *, const char, const char);
static int QueryIndex(const char);
static int PushOperatorToStack(OperatorStack *, const char);
static int PushOperandToStack(OperandStack *, const double);
static int PopOperatorFromStack(OperatorStack *, char *);
static int PopOperandFromStack(OperandStack *, double *);
static char PeakTopElementOfOperatorStack(const OperatorStack *);
static double PeakTopElementOfOperandStack(const OperandStack *);
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
static int SetAngleMode(Parameter *);
static void HelpCalculator(void);
static double ReadFirstFloat(char **);
static double ReadConstant(const Parameter *, char **);
static int UnaryOperation(const Parameter *, const char, const double, double *);
static int BinaryOperation(const double, const char, const double, double *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ExpressionCalculator(void)
{
    /*
     * Data declaration and initialization
     */
    double operandStackSpace[100] = {0.0};
    OperandStack theOperandStack = {
        .base = operandStackSpace,
        .top = operandStackSpace,
        .stacksize = sizeof operandStackSpace - 1
    };
    char operatorStackSpace[100] = {'\0'};
    OperatorStack theOperatorStack = {
        .base = operatorStackSpace,
        .top = operatorStackSpace,
        .stacksize = sizeof operatorStackSpace - 1,
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
    Parameter theParameter = {
        .pi = acos(-1),
        .answer = 0.0,
        .radianMode = 1,
        .angleFactor = 1.0
    };
    String currentLine = {'\0'}; /* store the input information */
    /* main loop */
    fprintf(stdout, "**********************************************************\n\n");
    fprintf(stdout, "Enter 'help' for a brief user manual of calculator\n");
    fprintf(stdout, "**********************************************************\n\n");
    while (1) {
        fprintf(stdout, "\nArtraCFD Calculator<< ");
        Fgets(currentLine, sizeof currentLine, stdin);
        CommandLineProcessor(currentLine); /* process current command */
        fprintf(stdout, "\n");
        if (0 == strncmp(currentLine, "help", sizeof currentLine)) {
            HelpCalculator();
            continue;
        }
        if (0 == strncmp(currentLine, "set", sizeof currentLine)) {
            SetAngleMode(&theParameter);
            continue;
        }
        if (0 == strncmp(currentLine, "exit", sizeof currentLine)) {
            return 0;
        } 
        if ('\0' == currentLine[0]) { /* no useful information in the command */
            fprintf(stdout, "\n");
            continue;
        }
        /* if non of above is true, then compute the expression */
        ComputeExpression(currentLine, &theOperandStack, &theOperatorStack, &theParameter);
    }
}
static int ComputeExpression(const char *currentLine, OperandStack *operandStack, 
        OperatorStack *operatorStack, Parameter *parameter)
{
    /*
     *  Every call will define a new space to store the command to avoid
     *  unclear interference between each inputted data which has different
     *  content and lengths.
     */
    String expression = {'\0'};
    /*
     * Use as a flag to judge whether '+' and '-' are unary operator based on 
     * the assumption that if "+" "-" appears as unary operator then it must
     * after a paratheses except at the beginning of an expression. By puting
     * a "(" at the initial of the expression, then this rule is always true.
     */
    expression[0] = '(';
    /* target command begins at the second position of storage address */
    char *command = expression + 1; 
    /*
     * Copy the information in current line to expression space, because need
     * to access command[2] and even later data when search and convert the
     * command to math expression, reduce and limit the available space for
     * copy is a safety choice.
     */
    strncpy(command, currentLine, (int)(sizeof expression) - 5);
    /*
     * The flow control of this program is important, thus, every function
     * call which may result an important error will be monitored.
     * These functional functions return 0 means true, 1 means false.
     */
    /*
     * Translate the command to math expression
     */
    if (0 != TranslateCommandToMathExpression(command)) {
        return 1; /* failed */
    }
    char *pointer = command; /* pointer to point the expression */
    char headOperator = '\0'; /* store the top operator of operator stack */
    char currentOperator = '\0'; /* store the current operator in command */
    double currentOperand = 0.0; /* store the current operand in command */
    double operandA = 0.0; /* first top operand in stack */
    double operandB = 0.0; /* second top operand in stack */
    /*
     * Always initialize and reset the stack status
     */
    operandStack->top = operandStack->base;
    operatorStack->top = operatorStack->base;
    PushOperatorToStack(operatorStack,'\0'); /* push a end flag to the stack of operator */
    /*
     * Calculation loop
     */
    while (('\0' != pointer[0]) || ('\0' != PeakTopElementOfOperatorStack(operatorStack))) {
        if (0 == IsDigit(pointer[0])) { /* find a operand */
            /* 
             * Read this float to current operand, note the read function will
             * update the pointer to the first character after the float data.
             */
            currentOperand = ReadFirstFloat(&pointer); 
            if (0 != PushOperandToStack(operandStack, currentOperand)) {
                return 1; /* failed */
            }
            continue;
        }
        if (0 == IsConstant(pointer[0])) { /* find a constant */
            /* 
             * Read this constant to current operand, note the read function will
             * update the pointer to the first character after the constant.
             */
            currentOperand = ReadConstant(parameter, &pointer);
            if (0 != PushOperandToStack(operandStack, currentOperand)) {
                return 1; /* failed */
            }
            continue;
        }
        /*
         * Now, treat everything left as an operator
         */
        if (0 != IsOperator(pointer[0])) {
            ShowInformation("undefined operator in expression");
            return 1;
        }
        currentOperator = pointer[0];
        switch (QueryPriority(operatorStack, PeakTopElementOfOperatorStack(operatorStack), currentOperator)) {
            case '<': 
                /*
                 * The priority of the head operator in the stack is lower
                 * than current operator in command, thus need to push
                 * current operator into stack.
                 */
                if (0 == IsDualOperatorActAsUnary(pointer)) {
                    /*
                     * Dual operator (they can both be unary and binary )
                     * '+' '-' show as a unary operator, then push a 0 to
                     * operand stack to make them become binary operator 
                     */
                    if (0 != PushOperandToStack(operandStack, 0)) {
                        return 1;
                    }
                }
                if (0 != PushOperatorToStack(operatorStack, currentOperator)) {
                    return 1;
                }
                ++pointer;
                break;
            case '=': 
                /*
                 * If the head operator in stack and current operator in
                 * command have the same priority, they both need to be
                 * dumped since they are control operators like parentheses
                 * and end tag.
                 */
                if (0 != PopOperatorFromStack(operatorStack, &currentOperator)) {
                    return 1;
                }
                ++pointer;
                break;
            case '>': 
                /*
                 * If the head operator in stack has higher priority than
                 * the current operator in command, then the current
                 * operator need to wait and can not go into the stack.
                 * At the same time, the head operator need to finish its
                 * calculation.
                 */
                if (0 != PopOperatorFromStack(operatorStack, &headOperator)) {
                    return 1;
                }
                if (0 == IsPureUnaryOperator(headOperator)) { /* unary operator */
                    if (0 != PopOperandFromStack(operandStack, &operandA)) {
                        return 1;
                    }
                    if (0 != UnaryOperation(parameter, headOperator, operandA, &currentOperand)) {
                        return 1;
                    }
                    if (0 != PushOperandToStack(operandStack, currentOperand)) {
                        return 1;
                    }
                } else {/* binary operator */
                    if (0 != PopOperandFromStack(operandStack, &operandA)) {
                        return 1;
                    }
                    if (0 != PopOperandFromStack(operandStack, &operandB)) {
                        return 1;
                    }
                    if (0 != BinaryOperation(operandB, headOperator, operandA, &currentOperand)) {
                        return 1;
                    }
                    if (0 != PushOperandToStack(operandStack, currentOperand)) {
                        return 1;
                    }
                }
                break;
            default: 
                ShowInformation("can not match parenthesis");
                return 1;
        }
    }
    /*
     * If the loop successfully exit, then it means the command,
     * and operator stack are all processed, now need to check the
     * operand stack, if there are more than one element left, 
     * it means something wrong happened.
     */
    if (1 != (operandStack->top - operandStack->base)) {
        ShowInformation("error, wrong expression");
        parameter->answer = 0.0; /* reset answer */
        return 1;
    }
    /*
     * Finally, get the final answer
     */
    /* save the result to answer */
    parameter->answer = PeakTopElementOfOperandStack(operandStack); 
    /* output the results */
    fprintf(stdout, "ans = %.6g\n", parameter->answer);
    return 0;
}
/*
 * Obtain the priority between two operators from the priority matrix
 */
static char QueryPriority(const OperatorStack *operatorStack, const char operatorA, const char operatorB)
{
    const int i = QueryIndex(operatorA);
    const int j = QueryIndex(operatorB);
    return operatorStack->priority[i][j];
}
/*
 * Query the index of a operator in the priority matrix
 */
static int QueryIndex(const char operator)
{
    int i = 0;
    switch(operator) {
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
            ShowInformation("unidentified operator"); 
            break;
    }
    return i;
}
/*
 * Push an element to the operand stack
 */
static int PushOperandToStack(OperandStack *operandStack, const double currentOperand)
{
    if ((operandStack->top - operandStack->base) >= operandStack->stacksize) {
        ShowInformation("operand stack is overflowing...");
        return 1;
    }
    operandStack->top[0] = currentOperand;
    ++(operandStack->top);
    return 0;
}
/*
 * Pop an element from operand stack
 */
static int PopOperandFromStack(OperandStack *operandStack, double *operandAddress)
{
    if (operandStack->top == operandStack->base) {
        ShowInformation("no sufficient operands in expression...");
        return 1;
    }
    --(operandStack->top);
    operandAddress[0] = operandStack->top[0];
    return 0;
}
/*
 * Get the top element from operand stack
 */
static double PeakTopElementOfOperandStack(const OperandStack *operandStack)
{
    return operandStack->top[-1];
}
/*
 * Push an element to the operator stack
 */
static int PushOperatorToStack(OperatorStack *operatorStack, const char currentOperator)
{
    if ((operatorStack->top - operatorStack->base) >= operatorStack->stacksize) {
        ShowInformation("operator stack is overflowing...");
        return 1;
    }
    operatorStack->top[0] = currentOperator;
    ++(operatorStack->top);
    return 0;
}
/*
 * Pop an element from the operator stack
 */
static int PopOperatorFromStack(OperatorStack *operatorStack, char *operatorAddress)
{
    if (operatorStack->top == operatorStack->base) {
        ShowInformation("no sufficient operators in expression...");
        return 1;
    }
    --(operatorStack->top);
    operatorAddress[0] = operatorStack->top[0];
    return 0;
}
/*
 * Get the top element from the operator stack
 */
static char PeakTopElementOfOperatorStack(const OperatorStack *operatorStack)
{
    return operatorStack->top[-1];
}
/*
 * Translate the input command to specific format that program can recognize
 * every operator correctly.
 */
static int TranslateCommandToMathExpression(char *command)
{
    char *scanner = command; /* scanner of the string */
    char *receiver = command; /* receiver to rewrite the string */
    int test = 1; /* test condition, default is false 1 */
    while ('\0' != scanner[0]) {
        switch (scanner[0]) {
            case 'e': 
                if (('x' == scanner[1]) && ('p' == scanner[2])) { /* use "e" stands for "exp" */
                    receiver[0] = 'e';
                    ++receiver; /* update receiver to newest receiving position */
                    scanner = scanner + 3; /* update scanner to newest scanning position */
                } else {
                    ShowInformation("unknown operator: e..");
                    return 1;
                }
                break;
            case 'l':
                if (('n' == scanner[1])) {/* use "n" stands for "ln" */
                    receiver[0] = 'n';
                    ++receiver;
                    scanner = scanner + 2;
                } else {
                    if (('g' == scanner[1])) { /* use "g" stands for "lg" */
                        receiver[0] = 'g';
                        ++receiver;
                        scanner = scanner + 2;
                    } else {
                        ShowInformation("unknown operator: l..");
                        return 1;
                    }
                }
                break;
            case 'a': 
                if (('b' == scanner[1]) && ('s' == scanner[2])) { /* use "a" stands for "abs" */
                    receiver[0] = 'a';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    if (('n' == scanner[1]) && ('s' == scanner[2])) { /* use 'q' for keyword "ans" */
                        receiver[0] = 'q';
                        ++receiver;
                        scanner = scanner + 3;
                    } else {
                        ShowInformation("unknown operator: a..");
                        return 1;
                    }
                }
                break;
            case 's': 
                if (('i' == scanner[1]) && ('n' == scanner[2])) { /* use "s" stands for "sin" */
                    receiver[0] = 's';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    ShowInformation("unknown operator: s..");
                    return 1;
                }
                break;
            case 'c': 
                if (('o' == scanner[1]) && ('s' == scanner[2])) { /* use "c" stands for "cos" */
                    receiver[0] = 'c';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    ShowInformation("unknown operator: c..");
                    return 1;
                }
                break;
            case 't': 
                if (('a' == scanner[1]) && ('n' == scanner[2])) { /* use "t" stands for "tan" */
                    receiver[0] = 't';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    ShowInformation("unknown operator: t..");
                    return 1;
                }
                break;
            case 'p': 
                if (('i' == scanner[1])) { /* use "p" stands for "pi" */
                    receiver[0] = 'p';
                    ++receiver;
                    scanner = scanner + 2;
                } else {
                    ShowInformation("unknown operator: p..");
                    return 1;
                }
                break;
            case '.':
                if (0 != (IsDigit(scanner[1]) + IsDigit(scanner[-1]))) {
                    ShowInformation(". is preceded or followed by non digits");
                    return 1;
                }
                receiver[0] = scanner[0];
                ++receiver;
                ++scanner;
                break;
            case ' ': 
                ++scanner;  /* ignore space */
                break;
            default:
                test = IsDigit(scanner[0]) * IsPureBinaryOperator(scanner[0]) * 
                    IsDualOperator(scanner[0]) * IsLeftParentheses(scanner[0]) * 
                    IsRightParentheses(scanner[0]); /* if any one is true, it's true */
                if (0 == test) { /* legal input single character */
                    receiver[0] = scanner[0];
                    ++receiver;
                    ++scanner;
                } else {
                    ShowInformation("unknown operator in expression");
                    return 1;
                }
                break;
        }
    }
    /* add a '\0' at the end of translated expression */
    receiver[0]='\0';
    return 0;
}
static int IsOperator(const char character)
{
    int testCondition = 1; /* default is false 1 */
    testCondition = IsPureUnaryOperator(character) *
        IsPureBinaryOperator(character) * IsDualOperator(character) *
        IsLeftParentheses(character) * IsRightParentheses(character) *
        IsEndTag(character); /* if any one is true, it's true */
    if (0 == testCondition) {
        return 0;
    }
    return 1;
}
static int IsPureUnaryOperator(const char character)
{
    if (('e' == character) || ('n' == character) || ('g' == character) ||
            ('a' == character) || ('s' == character) ||
            ('c' == character) || ('t' == character)) {
        return 0;
    }
    return 1;
}
static int IsDualOperator(const char character)
{
    if (('+' == character) || ('-' == character)) {
        return 0;
    }
    return 1;
}
static int IsDualOperatorActAsUnary(const char *pointer)
{
    int testCondition = IsDualOperator(pointer[0]) + 
        IsLeftParentheses(pointer[-1]); /* if both are true, it's true */
    if (0 == testCondition) { /* it means is a unary operator */
        return 0;
    }
    return 1;
}
static int IsPureBinaryOperator(const char character)
{
    if (('*' == character) || ('/' == character)
            || ('^' == character)) {
        return 0;
    }
    return 1;
}
static int IsLeftParentheses(const char character)
{
    if (('(' == character) || ('[' == character) || ('{' == character)) {
        return 0;
    }
    return 1;
}
static int IsRightParentheses(const char character)
{
    if ((')' == character) || (']' == character) || ('}' == character)) {
        return 0;
    }
    return 1;
}
static int IsConstant(const char character)
{
    if (('p' == character) || ('q' == character)) {
        return 0;
    }
    return 1;
}
static int IsEndTag(const char character)
{
    if (('\0' == character)) {
        return 0;
    }
    return 1;
}
static int IsDigit(const char character)
{
    if (('0' <= character) && ('9' >= character)) {
        return 0;
    } 
    return 1;
}
static int IsDot(const char character)
{
    if ('.' == character) {
        return 0;
    }
    return 1;
}
static double ReadFirstFloat(char **pointerAddress)
{
    int nscan = 0; /* read conversion count */
    char *string = *pointerAddress; /* copy the command address */
    double operand = 0.0;
    /* first, read a float to operand */
    nscan = sscanf(string, "%lg", &operand);
    VerifyReadConversion(nscan, 1);
    /* then update the pointer to latest position*/
    while (0 == IsDigit(string[0])) {
        ++string;
    }
    if (0 == IsDot(string[0])) {
        ++string;
    }
    while (0 == IsDigit(string[0])) {
        ++string;
    }
    *pointerAddress = string; /* get the updated address */
    return operand; /* return the float value */
}
static double ReadConstant(const Parameter *parameter, char **pointerAddress)
{
    char *string = *pointerAddress; /* copy the command address */
    double operand = 0.0;
    switch (string[0]) {
        case 'p': 
            operand = parameter->pi;
            ++string;
            break;
        case 'q': 
            operand = parameter->answer;
            ++string;
            break;
        default: 
            ShowInformation("undefined constant value");
            break;
    }
    *pointerAddress = string; /* get the updated address */
    return operand; /* return the float value */
}
static int UnaryOperation(const Parameter *parameter, 
        const char headOperator, const double operandA, double *currentOperandAddress)
{
    switch (headOperator) {
        case 'e':
            currentOperandAddress[0] = exp(operandA);
            break;
        case 'n':
            if (0.0 >= operandA) {
                ShowInformation("negative argument of ln(x)");
                currentOperandAddress[0] = 0;
                return 1;
            }
            currentOperandAddress[0] = log(operandA);
            break;
        case 'g':
            if (0.0 >= operandA) {
                ShowInformation("negative argument of lg(x)");
                currentOperandAddress[0] = 0;
                return 1;
            }
            currentOperandAddress[0] = log10(operandA);
            break;
        case 'a':
            currentOperandAddress[0] = fabs(operandA);
            break;
        case 's':
            currentOperandAddress[0] = sin(operandA * parameter->angleFactor);
            break;
        case 'c':
            currentOperandAddress[0] = cos(operandA * parameter->angleFactor);
            break;
        case 't':
            if (0.0 == cos(operandA * parameter->angleFactor)) {
                ShowInformation("negative argument of tangent");
                currentOperandAddress[0] = 0;
                return 1;
            }
            currentOperandAddress[0] = sin(operandA * parameter->angleFactor) / cos(operandA * parameter->angleFactor);
            break;
        default: 
            currentOperandAddress[0] = 0;
            return 1;
    }
    return 0;
}
static int BinaryOperation(const double operandB, 
        const char headOperator, const double operandA, double *currentOperandAddress)
{
    switch(headOperator)
    {
        case '+': 
            currentOperandAddress[0] = operandB + operandA;
            break;
        case '-':
            currentOperandAddress[0] = operandB - operandA;
            break;
        case '*':
            currentOperandAddress[0] = operandB * operandA;
            break;
        case '/':
            if (0.0 == operandA) {
                ShowInformation("negative argument of divide");
                currentOperandAddress[0] = 0;
                return 1;
            }
            currentOperandAddress[0] = operandB / operandA;
            break;
        case '^':
            currentOperandAddress[0] = pow(operandB, operandA);
            break;
        default: 
            currentOperandAddress[0] = 0;
            return 1;
    }
    return 0;
}
static int SetAngleMode(Parameter *parameter)
{
    fprintf(stdout, "Set mode by order number\n");
    fprintf(stdout, "\n 1 Angle in radian\n 2 Angle in degree\n\nSet:");
    String currentLine = {'\0'}; /* store current line */
    int nscan = 0; /* read conversion count */
    Fgets(currentLine, sizeof currentLine, stdin);
    nscan = sscanf(currentLine, "%d", &(parameter->radianMode));
    VerifyReadConversion(nscan, 1);
    fprintf(stdout, "\n");
    if (1 == parameter->radianMode) {
        parameter->angleFactor = 1;
        ShowInformation("*** Set mode: angle in radian ***");
    } else {
        if (2 == parameter->radianMode) {
            parameter->angleFactor = parameter->pi / 180;
            ShowInformation("*** Set mode: angle in degree ***");
        } else {
            ShowInformation("warning, unknown command...");
            parameter->angleFactor = 1;
            ShowInformation("*** Reset to default mode: angle in radian ***");
        }
    }
    return 0;
}
static void HelpCalculator(void)
{
    fprintf(stdout, "Operation options:\n\n");
    fprintf(stdout, "[help]              show this information\n");
    fprintf(stdout, "[set]               set angle mode in radian (default) or degree\n");
    fprintf(stdout, "math expression     calculator the inputted math expression\n");
    fprintf(stdout, "[exit]              return to ArtraCFD\n\n");
    fprintf(stdout, "                   Expression Calculator\n");
    fprintf(stdout, "Notice: Please avoid ambiguity expressions and use parenthesis:\n");
    fprintf(stdout, "        \'()\',\'[]\',\'{}\' to make semantic clear\n");
    fprintf(stdout, "Support: +, -, *, /, x^y, exp(x), ln(x), lg(x), abs(x), sin(a), cos(a), tan(a), pi\n");
    fprintf(stdout, "         where, x,y are numbers or expressions; a is a degree or radian;\n");
    fprintf(stdout, "         keyword \'ans\' is used to access the calculated value of last expression;\n");
    fprintf(stdout, "         space and tab are ignored in expression calculation;\n");
    fprintf(stdout, "Example: 1.5*sin(-pi/6)-[cos(pi/3)]^2+ln{exp[5*lg(abs(-100))]}\n\n");
}
/* a good practice: end file with a newline */

