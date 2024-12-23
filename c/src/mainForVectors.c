#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

int main(){

    Point* A = createPoint(-1,4); // A(2,3)
    Point* B = createPoint(3,2); // B(-1,6)

    printf("*------------------------------------------------------*\n"); 
    printf("Tests for Vectors ! \nStep one, creation.\n");

    Vector* firstVector = createVector(A,B); 
    printVector(firstVector); 

    printf("*------------------------------------------------------*\n"); 
    printf("Step Two : Print Function (not as a graph tho :/)\n");

    printFunction(firstVector->func, -5.00, 5.00, 0.5); 

    printf("*------------------------------------------------------*\n"); 
    printf("Step Three : Create Own Vectors\n");


    Vector* secondVector = createUserVector(); 
    printVector(secondVector);
    printFunction(secondVector->func,-4.00,-1.00,0.1);

    printf("*------------------------------------------------------*\n"); 

/**
 * Okayyyyyyyyyyy now i can create vectors, and derivate them 
 * Firstly, I wanted to estimate the 2nd degree function with affine functions represented by vectors
 * Now, I noticed that if i take two points with a 1 foot step it wouldnt work 
 * so, I thought that i can take ONE point by creating a vector with a step of +- 0.01
 * It's now correct if we take the abs of the derivate (which's returned in float)
 */


    printf("*------------------------------------------------------*\n"); 

    
    free(A);
    free(B);
    free(firstVector);
    free(secondVector);
    return 0; 
}