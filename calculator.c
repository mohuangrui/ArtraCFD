/****************************************************************************
 * Expression Calculator                                                    *
 * Last-modified: 20 Jan 2015 09:53:39 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a calculator for mathematical expressions            *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "calculator.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
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
    double *base; /* pointer to stack bottom */
    double *top; /* pointer to stack top */
    const int stacksize;
} OperandStack;
/*
 * Operator
 */
typedef struct {
    char *base; /* pointer to stack bottom */
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
static char GetTopElementOfOperatorStack(const OperatorStack *);
static double GetTopElementOfOperandStack(const OperandStack *);
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
    double operandStackSpace[100] = {0};
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
        .pi = 3.14159265358979323846264338327950288419716939937510,
        .answer = 0,
        .radianMode = 1,
        .angleFactor = 1.0
    };
    char currentLine[250] = {'\0'}; /* store the input information */
    /* main loop */
    printf("**********************************************************\n\n");
    printf("Enter 'help' for a brief user manual of calculator\n");
    printf("**********************************************************\n\n");
    while (1) {
        printf("\nArtraCFD Calculator<< ");
        fgets(currentLine, sizeof currentLine, stdin); /* read a line */
        CommandLineProcessor(currentLine); /* process current command */
        printf("\n");
        if (strncmp(currentLine, "help", sizeof currentLine) == 0) {
            HelpCalculator();
            continue;
        }
        if (strncmp(currentLine, "set", sizeof currentLine) == 0) {
            SetAngleMode(&theParameter);
            continue;
        }
        if (strncmp(currentLine, "quit", sizeof currentLine) == 0) {
            return 0;
        } 
        if (currentLine[0] == '\0') { /* no useful information in the command */
            printf("\n");
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
     *  unclear interference between each inputted data. 
     */
    char expression[250] = {'\0'};
    /*
     * Use as a flag to judge whether '+' and '-' are unary operator based on 
     * the assumption that if "+" "-" appears as unary operator then it must
     * after a "(" except at the beginning of expression.
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
    strncpy(command, currentLine, sizeof expression - 5);
    /*
     * The flow control of this program is important, thus, every function
     * call which may result an important error will be monitored.
     * These functional functions return 0 means true, 1 means false.
     */
    /*
     * Translate the command to math expression
     */
    if (TranslateCommandToMathExpression(command) != 0) {
        return 1; /* failed */
    }
    char *pointer = command; /* pointer to point the expression */
    char headOperator = '\0'; /* store the top operator of operator stack */
    char currentOperator = '\0'; /* store the current operator in command */
    double currentOperand = 0; /* store the current operand in command */
    double operandA = 0; /* first top operand in stack */
    double operandB = 0; /* second top operand in stack */
    /*
     * Always initialize and reset the stack status
     */
    operandStack->top = operandStack->base;
    operatorStack->top = operatorStack->base;
    PushOperatorToStack(operatorStack,'\0'); /* push a end flag to the stack of operator */
    /*
     * Calculation loop
     */
    while ((*pointer != '\0') || (GetTopElementOfOperatorStack(operatorStack) != '\0')) {
        if (IsDigit(*pointer) == 0) { /* find a operand */
            /* 
             * read this float to current operand, note the read function will
             * update the pointer to the first character after the float data.
             */
            currentOperand = ReadFirstFloat(&pointer); 
            if (PushOperandToStack(operandStack, currentOperand) != 0) {
                return 1; /* failed */
            }
            continue;
        }
        if (IsConstant(*pointer) == 0) { /* find a constant */
            /* 
             * read this constant to current operand, note the read function will
             * update the pointer to the first character after the constant.
             */
            currentOperand = ReadConstant(parameter, &pointer);
            if (PushOperandToStack(operandStack, currentOperand) != 0) {
                return 1; /* failed */
            }
            continue;
        }
        /*
         * now, treat everything left as an operator
         */
        if (IsOperator(*pointer) != 0) {
            ShowInformation("Undefined operator");
            return 1;
        }
        currentOperator = *pointer;
        switch (QueryPriority(operatorStack, GetTopElementOfOperatorStack(operatorStack), currentOperator)) {
            case '<': 
                /*
                 * the priority of the head operator in the stack is lower
                 * than current operator in command, thus need to push
                 * current operator into stack.
                 */
                if (IsDualOperatorActAsUnary(pointer) == 0) {
                    /*
                     * dual operator (they can both be unary and binary )
                     * '+' '-' show as a unary operator, then push a 0 to
                     * operand stack to make them become binary operator 
                     */
                    if (PushOperandToStack(operandStack, 0) != 0) {
                        return 1;
                    }
                }
                if (PushOperatorToStack(operatorStack, currentOperator) != 0) {
                    return 1;
                }
                ++pointer;
                break;
            case '=': 
                /*
                 * if the head operator in stack and current operator in
                 * command have the same priority, they both need to be
                 * dumped since they are control operators like parentheses
                 * and end tag.
                 */
                if (PopOperatorFromStack(operatorStack, &currentOperator) != 0) {
                    return 1;
                }
                ++pointer;
                break;
            case '>': 
                /*
                 * if the head operator in stack has higher priority than
                 * the current operator in command, then the current
                 * operator need to wait and can not go into the stack.
                 * At the same time, the head operator need to finish its
                 * calculation.
                 */
                if (PopOperatorFromStack(operatorStack, &headOperator) != 0) {
                    return 1;
                }
                if (IsPureUnaryOperator(headOperator) == 0) { /* unary operator */
                    if (PopOperandFromStack(operandStack, &operandA) != 0) {
                        return 1;
                    }
                    if (UnaryOperation(parameter, headOperator, operandA, &currentOperand) !=  0) {
                        return 1;
                    }
                    if (PushOperandToStack(operandStack, currentOperand) != 0) {
                        return 1;
                    }
                } else {/* binary operator */
                    if (PopOperandFromStack(operandStack, &operandA) != 0) {
                        return 1;
                    }
                    if (PopOperandFromStack(operandStack, &operandB) != 0) {
                        return 1;
                    }
                    if (BinaryOperation(operandB, headOperator, operandA, &currentOperand) !=  0) {
                        return 1;
                    }
                    if (PushOperandToStack(operandStack, currentOperand) != 0) {
                        return 1;
                    }
                }
                break;
            default: 
                ShowInformation("Can't match parenthesis");
                return 1;
        }
    }
    /*
     * if the loop successfully exit, then it means the command,
     * and operator stack are all processed, now need to check the
     * operand stack, if there are more than one element left, 
     * it means something wrong happened.
     */
    if (operandStack->top - operandStack->base != 1) {
        ShowInformation("Error, wrong expression");
        parameter->answer = 0; /* reset answer */
        return 1;
    }
    /*
     * finally, get the final answer
     */
    /* save the result to answer */
    parameter->answer = GetTopElementOfOperandStack(operandStack); 
    /* output the results */
    printf("ans = %.6lg\n", parameter->answer);
    return 0;
}
/*
 * Obtain the priority between two operators from the priority matrix
 */
static char QueryPriority(const OperatorStack *operatorStack, const char operatorA, const char operatorB)
{
    int i = 0; /* index of a operator */
    int j = 0; /* index of a operator */
    i = QueryIndex(operatorA);
    j = QueryIndex(operatorB);
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
            ShowInformation("Unidentified operator"); 
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
        ShowInformation("Operand stack is overflowing...");
        return 1;
    }
    *operandStack->top = currentOperand;
    ++operandStack->top;
    return 0;
}
/*
 * Pop an element from operand stack
 */
static int PopOperandFromStack(OperandStack *operandStack, double *operandAddress)
{
    if (operandStack->top == operandStack->base) {
        ShowInformation("Encountering an empty operand stack...");
        return 1;
    }
    --operandStack->top;
    *operandAddress = *operandStack->top;
    return 0;
}
/*
 * Get the top element from operand stack
 */
static double GetTopElementOfOperandStack(const OperandStack *operandStack)
{
    return (*(operandStack->top - 1));
}
/*
 * Push an element to the operator stack
 */
static int PushOperatorToStack(OperatorStack *operatorStack, const char currentOperator)
{
    if ((operatorStack->top - operatorStack->base) >= operatorStack->stacksize) {
        ShowInformation("Operator stack is overflowing...");
        return 1;
    }
    *operatorStack->top = currentOperator;
    ++operatorStack->top;
    return 0;
}
/*
 * Pop an element from the operator stack
 */
static int PopOperatorFromStack(OperatorStack *operatorStack, char *operatorAddress)
{
    if (operatorStack->top == operatorStack->base) {
        ShowInformation("Encountering an empty operator stack...");
        return 1;
    }
    --operatorStack->top;
    *operatorAddress = *operatorStack->top;
    return 0;
}
/*
 * Get the top element from the operator stack
 */
static char GetTopElementOfOperatorStack(const OperatorStack *operatorStack)
{
    return (*(operatorStack->top - 1));
}
/*
 * Translate the input command to specific format that program can recognize
 * every operator correctly, such as use:
 * "e" "n" "g" stands for "exp" "ln" "lg" respectively.
 * Note: receiver always points to the newest receiving position.
 */
static int TranslateCommandToMathExpression(char *command)
{
    char *scanner = command; /* scanner of the string */
    char *receiver = command; /* receiver to rewrite the string */
    int test = 1; /* test condition, default is false 1 */
    while (scanner[0] != '\0') {
        switch (scanner[0]) {
            case 'e': 
                if ((scanner[1] == 'x') && (scanner[2] == 'p')) { /* use "e" stands for "exp" */
                    receiver[0] = 'e';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    ShowInformation("unknown operator: e..");
                    return 1;
                }
                break;
            case 'l':
                if ((scanner[1] == 'n')) {/* use "n" stands for "ln" */
                    receiver[0] = 'n';
                    ++receiver;
                    scanner = scanner + 2;
                } else {
                    if ((scanner[1] == 'g')) { /* use "g" stands for "lg" */
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
                if ((scanner[1] == 'b') && (scanner[2] == 's')) { /* use "a" stands for "abs" */
                    receiver[0] = 'a';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    if ((scanner[1] == 'n') && (scanner[2] == 's')) { /* use 'q' for keyword "ans" */
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
                if ((scanner[1] == 'i') && (scanner[2] == 'n')) { /* use "s" stands for "sin" */
                    receiver[0] = 's';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    ShowInformation("unknown operator: s..");
                    return 1;
                }
                break;
            case 'c': 
                if ((scanner[1] == 'o') && (scanner[2] == 's')) { /* use "c" stands for "cos" */
                    receiver[0] = 'c';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    ShowInformation("unknown operator: c..");
                    return 1;
                }
                break;
            case 't': 
                if ((scanner[1] == 'a') && (scanner[2] == 'n')) { /* use "t" stands for "tan" */
                    receiver[0] = 't';
                    ++receiver;
                    scanner = scanner + 3;
                } else {
                    ShowInformation("unknown operator: t..");
                    return 1;
                }
                break;
            case 'p': 
                if ((scanner[1] == 'i')) { /* use "p" stands for "pi" */
                    receiver[0] = 'p';
                    ++receiver;
                    scanner = scanner + 2;
                } else {
                    ShowInformation("unknown operator: p..");
                    return 1;
                }
                break;
            case '.':
                if ((IsDigit(scanner[1]) + IsDigit(scanner[-1])) != 0) {
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
                    IsRightParentheses(scanner[0]);
                if (test == 0) { /* legal input single character */
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
        IsEndTag(character);
    if (testCondition == 0) {
        return 0;
    }
    return 1;
}
static int IsPureUnaryOperator(const char character)
{
    if ((character == 'e') || (character == 'n') || (character == 'g') ||
            (character == 'a') || (character == 's') ||
            (character == 'c') || (character == 't')) {
        return 0;
    }
    return 1;
}
static int IsDualOperator(const char character)
{
    if ((character == '+') || (character == '-')) {
        return 0;
    }
    return 1;
}
static int IsDualOperatorActAsUnary(const char *pointer)
{
    int testCondition = IsDualOperator(*pointer) + 
        IsLeftParentheses(*(pointer - 1));
    if (testCondition == 0) { /* it means is a unary operator */
        return 0;
    }
    return 1;
}
static int IsPureBinaryOperator(const char character)
{
    if ((character == '*') || (character == '/')
            || (character == '^')) {
        return 0;
    }
    return 1;
}
static int IsLeftParentheses(const char character)
{
    if ((character == '(') || (character == '[') || (character == '{')) {
        return 0;
    }
    return 1;
}
static int IsRightParentheses(const char character)
{
    if ((character == ')') || (character == ']') || (character == '}')) {
        return 0;
    }
    return 1;
}
static int IsConstant(const char character)
{
    if ((character == 'p') || (character == 'q')) {
        return 0;
    }
    return 1;
}
static int IsEndTag(const char character)
{
    if ((character == '\0')) {
        return 0;
    }
    return 1;
}
static int IsDigit(const char character)
{
    if ((character >= '0') && (character <= '9')) {
        return 0;
    } 
    return 1;
}
static int IsDot(const char character)
{
    if (character == '.') {
        return 0;
    }
    return 1;
}
static double ReadFirstFloat(char **pointerAddress)
{
    char *string = *pointerAddress; /* copy the command address */
    double operand = 0;
    /* first, use sscanf read a float to operand */
    sscanf(string, "%lg", &operand);
    /* then update the pointer to latest position*/
    while (IsDigit(*string) == 0) {
        ++string;
    }
    if (IsDot(*string) == 0) {
        ++string;
    }
    while (IsDigit(*string) == 0) {
        ++string;
    }
    *pointerAddress = string; /* get the updated address */
    return operand; /* return the float value */
}
static double ReadConstant(const Parameter *parameter, char **pointerAddress)
{
    char *string = *pointerAddress; /* copy the command address */
    double operand = 0;
    switch (*string) {
        case 'p': 
            operand = parameter->pi;
            ++string;
            break;
        case 'q': 
            operand = parameter->answer;
            ++string;
            break;
        default: 
            ShowInformation("Undefined constant value");
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
            *currentOperandAddress = exp(operandA);
            break;
        case 'n':
            if (operandA <= 0) {
                ShowInformation("Negative argument of ln(x)");
                *currentOperandAddress = 0;
                return 1;
            }
            *currentOperandAddress = log(operandA);
            break;
        case 'g':
            if (operandA <= 0) {
                ShowInformation("Negative argument of lg(x)");
                *currentOperandAddress = 0;
                return 1;
            }
            *currentOperandAddress = log10(operandA);
            break;
        case 'a':
            *currentOperandAddress = fabs(operandA);
            break;
        case 's':
            *currentOperandAddress = sin(operandA * parameter->angleFactor);
            break;
        case 'c':
            *currentOperandAddress = cos(operandA * parameter->angleFactor);
            break;
        case 't':
            *currentOperandAddress = sin(operandA * parameter->angleFactor) / cos(operandA * parameter-> angleFactor);
            break;
        default: 
            *currentOperandAddress = 0;
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
            *currentOperandAddress = operandB + operandA;
            break;
        case '-':
            *currentOperandAddress = operandB - operandA;
            break;
        case '*':
            *currentOperandAddress = operandB * operandA;
            break;
        case '/':
            *currentOperandAddress = operandB / operandA;
            break;
        case '^':
            *currentOperandAddress = pow(operandB, operandA);
            break;
        default: 
            *currentOperandAddress = 0;
            return 1;
    }
    return 0;
}
static int SetAngleMode(Parameter *parameter)
{
    printf("Set mode by order number\n");
    printf("\n 1 Angle in radian\n 2 Angle in degree\n\nSet:");
    char currentLine[250] = {'\0'}; /* store current line */
    fgets(currentLine, sizeof currentLine, stdin); /* read a line */
    sscanf(currentLine, "%d", &(parameter->radianMode));
    printf("\n");
    if (parameter->radianMode == 1) {
        parameter->angleFactor = 1;
        ShowInformation("*** Set mode: angle in radian ***");
    } else {
        if (parameter->radianMode == 2) {
            parameter->angleFactor = parameter->pi / 180;
            ShowInformation("*** Set mode: angle in degree ***");
        } else {
            ShowInformation("Warning, unknown command...");
            parameter->angleFactor = 1;
            ShowInformation("*** Reset to default mode: angle in radian ***");
        }
    }
    return 0;
}
static void HelpCalculator(void)
{
    printf("Operation options:\n");
    printf("[help]              show this information\n");
    printf("[set]               set angle mode in radian (default) or degree\n");
    printf("math expression     calculator the inputted math expression\n");
    printf("[quit]              return to ArtraCFD\n\n");
    printf("                   Expression Calculator\n");
    printf("Notice: Please avoid ambiguity expressions and use parenthesis:\n");
    printf("        \'()\',\'[]\',\'{}\' to make semantic clear\n");
    printf("Support: +, -, *, /, x^y, exp(x), ln(x), lg(x), abs(x), sin(a), cos(a), tan(a), pi\n");
    printf("         where, x,y are numbers or expressions; a is a degree or radian;\n");
    printf("         keyword \'ans\' is used to access the calculated value of last expression;\n");
    printf("         space and tab are ignored in expression calculation;\n");
    printf("Example: 1.5*sin(-pi/6)-[cos(pi/3)]^2+ln{exp[5*lg(abs(-100))]}\n\n");
}
/* a good practice: end file with a newline */

