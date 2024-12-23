#include "functions.h"

/**
 * @brief Function used to create a point
 * @param x The x-axis 
 * @param y The y-axis
 * @return A pointer to the new point. 
 */
Point* createPoint(float x, float y){
    Point* newPoint = (Point*)malloc(sizeof(Point)); 
    newPoint->x = x; 
    newPoint->y = y; 
    return newPoint; 
}

/**
 * @brief Function used to print caracteristic of the vector
 * @param vector A pointer to the wanted vector 
 * @return A pointer to the first char of the string
 */
void printVector(Vector* vector){

    if (vector->isNull){
        printf("This vector's null");
        return; 
    }

    Point* A = vector->a;
    Point* B = vector->b; 
    
    Function* func = vector->func; 
 
    float alpha, beta, x1, x2, y1, y2;

    x1 = A->x; 
    x2 = B->x; 
    y1 = A->y; 
    y2 = B->y; 
    alpha = func->alpha; 
    beta = func->beta; 

    printf("This vector caracteristics are : \n");
    printf("A(%.1f , %.1f)\nB(%.1f , %.1f)\n",x1,y1,x2,y2);
    if (alpha == 1.0f) {
        if (beta == 0.0f) {
            printf("f(x) = x.\n");
        } else {
            printf("f(x) = x %c %.2f.\n", (beta > 0) ? '+' : '-', fabs(beta));
        }
    } else if (alpha == -1) {
        if (beta == 0) {
            printf("f(x) = -x.\n");
        } else {
            printf("f(x) = -x %c %.2f.\n", (beta > 0) ? '+' : '-', fabs(beta));
        }
    } else {
        if (beta == 0) {
            printf("f(x) = %.2fx.\n", alpha);
        } else {
            printf("f(x) = %.2fx %c %.2f.\n", alpha, (beta > 0) ? '+' : '-', fabs(beta));
        }
    }
}

/**
 * @brief Function used to create a new vector 
 * If you want to create a single point vector (using graph origin as a point) use 0 as a point :)
 * If point A and point B are missing, it creates an empty vector
 * @param a The first point of the vector
 * @param b A pointer to a second point.
 * @return A pointer to the new vector
 */
Vector* createVector(Point* a, Point* b){
    if (a == NULL && b == NULL){return (createNullVector());}
    Vector* newVector = (Vector*)malloc(sizeof(Vector)); 
    newVector->isNull = 0; 
    newVector->a = a; 
    newVector->b = b;
    findFunc(newVector); 
    return newVector; 
}

/**
 * @brief Function created when I thought we needed a different function to create null vectors
 * Used in createVector()
 * @return A point to the new vector.
 */
Vector* createNullVector(){
    Vector* newVector = (Vector*)malloc(sizeof(Vector));
    newVector->isNull = 1; 
    newVector->a = NULL;
    newVector->b = NULL; 
    newVector->func = NULL;
    return newVector;
}

/**
 * @brief Function used to create a new vector directly from an another one. 
 * It's the derivate vector 
 * @param vector The "father" vector 
 * @return A pointer to the "son" vector 
 */
Vector* createDerivateVectorFromVector(Vector* vector){

    if (vector->isNull || vector->func == NULL){
        fprintf(stderr,"The Vector or its function can't be null\n");
        exit(EXIT_FAILURE);
    }

    Vector* derivateVector = (Vector*)malloc(sizeof(Vector)); 
    Function* derivateFunction = (Function*)malloc(sizeof(Function));

    derivateFunction = calculateDerivate(vector,vector->a->x);

    derivateVector->isNull = 0; 
    derivateVector->a = vector->a; 
    derivateVector->b = vector->b; 
    derivateVector->func = derivateFunction;

    return derivateVector;
    
}

// void setVectorA(Vector* vector, Point* a){
//     vector->a = a;
// }

// void setVectorB(Vector* vector, Point* b){
//     vector->b =b; 
// }

int square(int x){
    return x*x;
}

/**
 * 
 */
int almostEquals(float a, float b){
    float epsilon = 1e-6; 
    return fabs(a-b) < epsilon; 
}

/**
 * @brief Function used to calculate the distance of the vector 
 * We assume that the vector is positioned at the center of the graph
 * @param vector The vector whose distance is calculated 
 * @return An int equals to the distance of the vector
 */
int calculateDistance(Vector* vector){

    Point* a = vector->a;
    Point* b = vector->b;
    float result; 

    if (a == NULL && b == NULL){
        printf("Distance of a null vector is not doable");
        return 1;  
    }
    if (a == NULL) {
        result = sqrt(square(b->x) + square(b->y));
        return result;
    } else {
        float x = b->x - a->x;
        float y = b->y - a->y;
        result = sqrt(square(x) + square(y));
        return result; 
    }

    fprintf(stderr, "Error : Vector is not recognized\n");
    exit(EXIT_FAILURE); 
}

// /**
//  * @brief Function used to set the parameter of the function correctly (affine function)
//  * @param func A pointer to the function you want to set
//  * @param alpha The director coeficient. Mandatory 
//  * @param beta The other coeficient
//  */
// void setFunction(Function* func, int* alpha, int beta){

//     if (alpha == NULL || alpha == 0){
//         fprintf(stderr,"The director coefficient can't be 0 or null.\n");
//         return; 
//     }

//     func->alpha= alpha;
//     func->beta = beta; 
// }

/**
 * @brief Function used to find the function associated to the vector 
 * We find this function by using the two points associated to the vector 
 * @param vector The vector
 * @return A struct "Function" used to determinate the function associated to the pointer
 */
void findFunc(Vector* vector){

    float alpha,beta; 
    Function* vectorFunction = (Function*)malloc(sizeof(Function));

    float x1 = vector->a->x;
    float x2 = vector->b->x; 
    float y1 = vector->a->y;
    float y2 = vector->b->y;

    if (x2 - x1 != 0.0f){
        alpha = (y2-y1) / (x2-x1);
    } else {
        fprintf(stderr,"Division by zero\n"); 
        exit(EXIT_FAILURE);
    }

    if (y1 != y2){
        beta = y1 - alpha * x1;
    } else {
        beta = 0; 
    }   

    vectorFunction->alpha = alpha; 
    vectorFunction->beta = beta;
    vector->func = vectorFunction; 
}

/**
 * @brief Function used to calculate a value for f(x)
 * @param func Pointer to the function we want to calculate the point.
 * @param x Value of the antecedant 
 * @return The value of f(x) (x image, in int) 
 */
float calculateAffineFunction(Function* func, float x){

    float alpha,beta;

    alpha = func->alpha; 
    beta = func->beta; 

    return (alpha * x + beta); 
}

/**
 * @brief Function used to derivate from a value 
 * As there is a lot of "classic" derivate functions, we'll use the finite difference method
 * f'(x) = (f(x+h) - f(x-h))/2*h where h's a small number beside 1 
 * Tangent's equation'll be used to find the function associated (speed in physics)
 * y = f'(a)(x-a) + f(a)
 * Knowing that y = f'(a) * (x-a) + f(a)
 *              y = f'(a) * x - f'(a) * a + f(a)
 * We also know f'(a), a and f(a) 
 *               <=> f'(a) * x - Cst + Cst
 * Cst + Cst = beta and f'(a) = alpha
 * <=> beta = f'(a) * a + f(a)
 * @param value The value you want to derivate
 * @param x The point to derivate on
 * @return A pointer to the tangent's function freshly derivated 
 */
Function* calculateDerivate(Vector* vector, int x){

    Function* function = vector->func; 

    if (vector->isNull || function == NULL){
        fprintf(stderr,"This vector (or it's function) is null\n");
        return NULL;
    }   
    
    Function* derivateFunction = (Function*)malloc(sizeof(Function));
    
    int epsilon, fx ,fx1, fx2, alpha, beta; 
    epsilon = 1e-3; 

    fx = calculateAffineFunction(function,x);
    fx1 = calculateAffineFunction(function,(x+epsilon)); 
    fx2 = calculateAffineFunction(function,(x-epsilon)); 

    alpha = (fx1 - fx2) / 2 * epsilon;

    beta = alpha * x + fx; 

    derivateFunction->alpha = alpha;
    derivateFunction->beta = beta;
    
    return derivateFunction;
}

/**
 * @brief Function used to print a function board
 * @param func A pointer to the function
 * @param start A float number used to strart the x-values
 * @param end A float number used to end the x-values
 * @param step A float number used as a step for the function
 */
void printFunction(Function* func, float start, float end, float step){
    
    printf("--------------------------------------\n"); 
    printf(" x \t\t f(x) \t\t f'(x)\n");
    printf("--------------------------------------\n"); 

    float i; 

    for (i = start ; i <= end ; i += step){

        float valueForI, primeValueForI; 

        valueForI = calculateAffineFunction(func,i);
        primeValueForI = func->alpha; 

        printf("%.2f \t\t %.2f \t\t %.2f\n", i,valueForI,primeValueForI );
    }
}

/**
 * @brief Function used to allow user to create vectors
 * @return A pointer to the freshly create vector
 */
Vector* createUserVector(){

    float x1,x2,y1,y2; 

    printf("Welcome in the vectors creation pannel! \nPlease enter the x-axis for A : ?\n");
    scanf("%f",&x1);
    printf("Please enter the y-axis for A : ?\n");
    scanf("%f",&y1);
    printf("Please enter the x-axis for B : ?\n");
    scanf("%f",&x2); 
    printf("Please enter the y-axis for B : ?\n");
    scanf("%f",&y2);

    Point* A = createPoint(x1,y1); 
    Point* B = createPoint(x2,y2); 

    Vector* vector = createVector(A,B); 
    findFunc(vector); 
    return vector; 
}

/**
 * @brief Function used to create vectors 
 * Polynomious functions (2nd Degree) : ax² + bx + c
 * @param x2 A float (x² coefficient, a)
 * @param x A float (x coefficient, b)
 * @param c The constant value (c)
 * @param start A float to start the function need to be a round float
 * @param end A float to end the function need to be a round float
 * @return A tab of the vectors created
*/
Vector* createVectorsForPolynomiousFunction(float x2, float x, float c, float start, float end){
    
    float size = end - start;  
    Vector* vectors = (Vector*)malloc(size * sizeof(Vector));

    for (float i = start ; i < size ; i += 1.00){
        float xA = i - 0.01; 
        float xB =  i + 0.01; 
        
    }

    return vectors;
}