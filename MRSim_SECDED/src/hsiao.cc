/**
 * @file: hsiao.cc
 * @author: Jungrae Kim <dale40@gmail.com>
 * CODEC implementation (HSIAO)
 */

#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include "hsiao.hh"
#include <iostream>

//--------------------------------------------------------------------
//// From "A Class of Optimal Minimum Odd-weight-column SEC-DED Codes" / Hsiao / 1970
//// (72,64) code
//// 8C1 (=8) + 8C3 (=56) + 8C5 (=8 e.a.)
//// total weight = 8 + 3*56 + 5*8 = 216
//// average number of 1s in each row = 27
//// average number of 1s in Hb = 26
                    /*           8                          8                        8                     8                          8                        8                        8                        8            --> # of columns = 64(== bitK)  */
uint8_t HSIAO_PT[] = {1, 1, 1, 1, 1, 1, 1, 1,  0, 0, 1, 0, 0, 1, 1, 0,  0, 1, 0, 0, 1, 0, 0, 1,  1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 1,  0, 0, 0, 1, 1 ,1, 0, 0,  1, 1, 1, 0, 0, 0, 0, 0,   
                      1, 1, 1, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1,  0, 0, 1, 0, 0, 1, 1, 0,  0, 1, 0, 0, 1, 0, 0, 1,  1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 1,  0, 0, 0, 1, 1, 1, 0, 0,  //
                      0, 0, 0, 1, 1, 1, 0, 0,  1, 1, 1, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1,  0, 0, 1, 0, 0, 1, 1, 0,  0, 1, 0, 0, 1, 0, 0, 1,  1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 1,  //
                      0, 0, 0, 1, 0, 0, 1, 1,  0, 0, 0, 1, 1, 1, 0, 0,  1, 1, 1, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1,  0, 0, 1, 0, 0, 1, 1, 0,  0, 1, 0, 0, 1, 0, 0, 1,  1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0,  //  # of rows = 8 (== bitR)
                      0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 1,  0, 0, 0, 1, 1, 1, 0, 0,  1, 1, 1, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1,  0, 0, 1, 0, 0, 1, 1, 0,  0, 1, 0, 0, 1, 0, 0, 1,  1, 0, 0, 1, 0, 0, 0, 0,  // 
                      1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 1,  0, 0, 0, 1, 1, 1, 0, 0,  1, 1, 1, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1,  0, 0, 1, 0, 0, 1, 1, 0,  0, 1, 0, 0, 1, 0, 0, 1,  //
                      0, 1, 0, 0, 1, 0, 0, 1,  1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 1,  0, 0, 0, 1, 1, 1, 0, 0,  1, 1, 1, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1,  0, 0, 1, 0, 0, 1, 1, 0,  //
                      0, 0, 1, 0, 0, 1, 1, 0,  0, 1, 0, 0, 1, 0, 0, 1,  1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 1, 0, 0, 1, 1,  0, 0, 0, 1, 1, 1, 0, 0,  1, 1, 1, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 1, 1, 1};

Hsiao::Hsiao(const char *name, int _bitN, int _bitR)
    : BinaryLinearCodec(name, _bitN, _bitR) {
    
    /*
        To do:
        1. make generate matrix
        2. make parity check matrix

        Hint:
        G = [Ik: P], shpae is kxn
        H = [-P.T:I(n-k)], shape is nx(n-k)  
    */
    for(int col=0;col<bitK;col++){
        for(int row=0;row<bitN;row++){
            if(row<bitR){ gMatrix[col*bitN+row] = HSIAO_PT[row*bitK+col]; }
            else{ gMatrix[col*bitN+row] = ((row-bitR)==col)?1:0; }
        }
    }
    for(int col=0;col<bitR;col++){
        for(int row=0;row<bitN;row++){
            if(row<bitR){ hMatrix[col*bitN+row] = (row==col)?1:0;}
            else{ hMatrix[col*bitN+row] = HSIAO_PT[col*bitK+row-bitR];}
        }
    }
    
    //print(stdout);
    verifyMatrix(); // Automatically check whether yout matrix is correct or not.
}

ErrorType Hsiao::decode(ECCWord *msg, ECCWord *decoded, std::set<int>* correctedPos) {
    // step 1: copy the message data
    decoded->clone(msg);
    
    // step 2: generate syndrome
    bool synError = genSyndrome(msg);
    //std::cout << "The synError is : "<<synError << std::endl;


    // error is detectd.
    if(synError){
        // search syndrome in hMatrix using array member "syndrom"
        for(int col=0;col<bitN;col++){
            bool syncomp = true;
            for(int row=0;row<bitR;row++){
                if(syndrom[row] != hMatrix[row*bitN+col]){ syncomp = false; break; }
            }
            if(syncomp){
                decoded->invBit(col);
                if(decoded->isZero()) return CE;
                else return SDC;
            }
        }
        return DUE;
    }else{
        // Note that we assume non-error msg is filtered.
        return SDC;
    }
}

/**
 * 
 *  @return: wheter there is syndrome or not in received codeword.
 *           true: error is detected in msg.
 *           false: error is not detected in msg.
 * 
 */
bool Hsiao::genSyndrome(ECCWord *msg) {
    bool synError = false;
    for(int row=0;row<bitR;row++){
        syndrom[row]=0;
        for(int col=0;col<bitN;col++){
            syndrom[row] ^= hMatrix[row*bitN+col]&msg->getBit(col);
        }
        if(syndrom[row]) synError = true;
    }     
    return synError;
}


/**
 * @brief Determine whether gMatrix and hMatrix is correct or not.
 * 
 * 
 */
void Hsiao::verifyMatrix()
{
    int max_row_weight = 0;
    int min_row_weight = INT_MAX;
    for (int i=0; i<bitR; i++) {
        int row_weight = 0;
        for (int j=0; j<bitN; j++) {
            row_weight += hMatrix[i*bitN+j];
        }
        if (row_weight > max_row_weight) max_row_weight = row_weight;
        if (row_weight < min_row_weight) min_row_weight = row_weight;
    }

    if ((max_row_weight - min_row_weight) > 1) {
        fprintf(stderr, "Invalid Hsiao matrix: Max row weight: %d Min row weight: %d\n", max_row_weight, min_row_weight);
        exit(-1);
    }

    // G x Ht = 0
    for (int i=0; i<bitK; i++) {
        for (int j=0; j<bitR; j++) {
            int temp = 0;
            for (int k=0; k<bitN; k++) {
                temp ^= gMatrix[i*bitN+k] & hMatrix[j*bitN+k];
            }
            if (temp!=0) {
                fprintf(stderr, "BAD H and G matrix %d\n", i);
                exit(-1);
            }
        }
    }

}

