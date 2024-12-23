#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Point{
    float x; 
    float y;
} Point; 

typedef struct Function{
    float alpha;
    float beta; 
} Function;

typedef struct Vector{
    Point* a;
    Point* b;  
    Function* func;     
    int isNull;
} Vector; 

Point* createPoint(float x, float y); 
void printVector(Vector* vector);
Vector* createVector(Point* a,Point* b);
Vector* createNullVector();
Vector* createDerivateVectorFromVector(Vector* vector); 
// void setVectorA(Vector* vector, Point* a); 
// void setVectorB(Vector* vector, Point* b); 
int square(int x); 
int almostEquals(float a, float b); 
int calculateDistance(Vector* vector); 
//void setFunction(Function* func, int alpha, int beta); 
void findFunc(Vector* vector);
float calculateAffineFunction(Function* func, float x);
Function* calculateDerivate(Vector* vector, int x); 
void printFunction(Function* func, float start, float end, float step); 
Vector* createUserVector(); 
Vector* createVectorsForPolynomiousFunction(float x2, float x, float c, float start, float end);

#endif