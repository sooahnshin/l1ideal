// This code is built upon the C code 'anominate.c' from the R package 'anominate' which is under GPL-2 license.
// Originally, 'anominate.c' implements univariate slice sampling to draw samples from the posterior distribution of their own model - alpha-nominate.
// Here, the original code has been modified to implement multivariate slice sampling using hyper-rectangle to draw samples from the posterior distribution of L1 norm multidimensional ideal point model.
// Original C code came from: Royce Carroll, Christopher Hare, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal (2017). Alpha-NOMINATE: Ideal Point Estimator. R package version 0.6. URL http://k7moa.com/alphanominate.htm

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#define MAXVOTES 100000
#define SLICE_W 4
#define _USE_MATH_DEFINES
#include <math.h>

typedef struct{
    int row, col, vote;
}cell;

/*  blockData is a generic structure used to pass information back and forth, mainly into the slicer
 'Data' is the raw data, 'ideal/yealoc/nayloc' are locations, 'length' is Data length
 'dim' is number of dimensions to be updated
 'beta' is a signal-to-noise ratio.  'weight' is the weight of the first dimension.
 */
typedef struct{
    cell **Data;
    double *ideal, *yealoc, *nayloc;
    int length, dim, nrow, ncol;
    double beta, weight;
} blockData;

typedef struct{
    cell ***rowData;
    double *ideal, *yealoc, *nayloc;
    int *rowLengths, nrow, ncol, dims;
    double beta, weight;
} betaBlock;

double *slice(double (*fp)(double *, void *), double *init, void *ptr, double w, int dims);
double slice_weight(double (*fp)(double *, void *), double *init, void *ptr, double w, int dims);
double slice_beta(double (*fp)(double *, void *), double *init, void *ptr, double w, int dims);
double l1LogLike(int vote, double *x, double *yea, double *nay, double beta, double weight, int dims);
double legisLoglike(double *par, void *legis_ptr);
double nayLoglike(double *par, void *nay_ptr);
double yeaLoglike(double *par, void *yea_ptr);
void readDataFromVector(int *inputVector, cell **data, int *nvotes, int *nrow, int *ncol);
void formatData(cell *data, int nvotes, int nrow, int ncol, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData);
void readInitsFromVector(double *initIdeal, double *initBill, double **idealpts, double **yealocs, double **naylocs, double **beta, double **weight, int nrow, int ncol, int dims);
void writeDataOutput(double *output, int *op, int nrow, int ncol, double **idealpts, double **yealocs,
                     double **naylocs, double **beta, double **weight,int dims);
void updateYea(cell ****colData,int **colLengths,double **idealpts,double **yealocs,double **naylocs,double **beta, double **weight,int nrow,int ncol,int dims);
void updateNay(cell ****colData,int **colLengths, double **idealpts,double **yealocs,double **naylocs,double **beta, double **weight,int nrow,int ncol,int dims);
void updateLegis(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta, double **weight,int nrow,int ncol,int dims);
void sampleData(int Nsamples,int nrow,int ncol,int **rowLengths,int **colLengths,cell ****rowData,cell ****colData,double **idealpts,double **yealocs,double **naylocs, double **beta, double **weight, int verbose, int thin, int dims, double *output);
void freeData(int nrow, int ncol, cell **data, double **idealpts, double **yealocs, double **naylocs, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData, double **beta, double **weight);

double weightLegisLL(double *par, void *legis_ptr);
double weightLogLike(double *par, void *weight_ptr);
void updateWeight(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta, double **weight,int nrow,int ncol,int dims);
double betaLegisLL(double *par, void *legis_ptr);
double betaLogLike(double *par, void *beta_ptr);
void updateBeta(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta, double **weight,int nrow,int ncol,int dims);
void Cl1ideal(int *dataVector, double *initIdeal, double *initBill, double *output, int *thin, int *ncol, int *nrow, int *Nsamples, int *dims, int *verbose);


double *slice(double (*fp)(double *, void *), double *init, void *ptr, double w, int dims){
    double *x,y,*L,*R;
    int i;
    short int flag;
    blockData *sliceData = (blockData *) ptr;
    
    x = (double *) malloc(dims * sizeof(double));
    L = (double *) malloc(dims * sizeof(double));
    R = (double *) malloc(dims * sizeof(double));
    
    y = fp(init, sliceData) - rexp(1);
    for(i=0;i<dims;i++){
        L[i] = init[i] - w*runif(0,1);
        R[i] = L[i] + w;
    }
    flag=1;
    while(flag==1) {
        for(i=0;i<dims;i++){
            x[i] = L[i] + runif(0,1)*(R[i]-L[i]);
        }
        if(y < fp(x,sliceData)) flag=0;
        for(i=0;i<dims;i++){
            if(x[i]<init[i]) L[i] = x[i];
            else R[i] = x[i];
        }
    }
    free(L);
    free(R);
    return(x);
}

double slice_weight(double (*fp)(double *, void *), double *init, void *ptr, double w, int dims){
    double x,y,L,R;
    short int flag;
    blockData *sliceData = (blockData *) ptr;
    
    y = fp(init, sliceData) - rexp(1);
    L = *init - w*runif(0,1);
    R = L + w;
    if(L<0) L=0;
    if(R>1) R=1;
    flag=1;
    while(flag==1) {
        x = L + runif(0,1)*(R-L);
        if(y < fp(&x,sliceData)) flag=0;
        if(x<*init) L = x;
        else R = x;
    }
    
    return(x);
}

double slice_beta(double (*fp)(double *, void *), double *init, void *ptr, double w, int dims){
    double x,y,L,R;
    short int flag;
    blockData *sliceData = (blockData *) ptr;
    
    y = fp(init, sliceData) - rexp(1);
    L = *init - w*runif(0,1);
    R = L + w;
    flag=1;
    while(flag==1) {
        x = L + runif(0,1)*(R-L);
        if(y < fp(&x,sliceData)) flag=0;
        if(x<*init) L = x;
        else R = x;
    }
    
    return(x);
}

double l1LogLike(int vote, double *x, double *yea, double *nay, double beta, double weight, int dims){
    
    double yeaUtil, nayUtil, diffUtil;
    double DistYea = 0, DistNay=0;
    
    if (dims==1) {
        DistYea = weight*fabs(x[0]-yea[0]);
        DistNay = weight*fabs(x[0]-nay[0]);
    } else {
        DistYea = weight*fabs(x[0]-yea[0])+(1-weight)*fabs(x[1]-yea[1]);
        DistNay = weight*fabs(x[0]-nay[0])+(1-weight)*fabs(x[1]-nay[1]);
    }
    
    yeaUtil = -beta*DistYea;
    nayUtil = -beta*DistNay;
    diffUtil = yeaUtil - nayUtil;
    
    if(vote==1) return(pnorm(diffUtil, 0, 1, 1, 1));
    if(vote==0) return(pnorm(-diffUtil, 0, 1, 1, 1));
    
    return(0);
}

//Takes one legislator (pair of ideal pts), all roll call locations, gives loglike
double legisLoglike(double *par, void *legis_ptr){
    
    int i,j,k;
    double loglike=0.0, *yea, *nay, prior=0.0;
    blockData *legisData = (blockData *)legis_ptr;
    
    yea = (double *) malloc((*legisData).dim * sizeof(double));
    nay = (double *) malloc((*legisData).dim * sizeof(double));
    
    for(i=0; i<(*legisData).length; i++) {
        
        for(j=0;j<(*legisData).dim;j++){
            yea[j] = (*legisData).yealoc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
            nay[j] = (*legisData).nayloc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
        }
        
        loglike += l1LogLike((*(*legisData).Data[i]).vote, par, yea, nay, (*legisData).beta, (*legisData).weight, (*legisData).dim);
    }
    
    for(k=0;k<(*legisData).dim;k++){
        prior += pow(par[k],2);
    }
    
    free(yea);
    free(nay);
    return(loglike - prior);
    
}

//Takes one pair yea/nay locations, all ideal pts, gives loglike
double nayLoglike(double *par, void *nay_ptr){
    
    int i,j,k;
    double loglike=0.0, *ideal, prior=0.0;
    blockData *nayData = (blockData *)nay_ptr;
    
    ideal = (double *) malloc((*nayData).dim * sizeof(double));
    
    for(i=0; i<(*nayData).length; i++){
        
        for(j=0;j<(*nayData).dim;j++){
            ideal[j] = (*nayData).ideal[(*(*nayData).Data[i]).row -1 + j*(*nayData).nrow];
        }
        loglike += l1LogLike((*(*nayData).Data[i]).vote, ideal, (*nayData).yealoc,par, (*nayData).beta, (*nayData).weight, (*nayData).dim);
    }
    
    for(k=0;k<(*nayData).dim;k++){
        prior += pow(par[k],2);
    }

    free(ideal);
    return(loglike - prior);
}


//Takes one pair yea/nay locations, all ideal pts, gives loglike
double yeaLoglike(double *par, void *yea_ptr){
    
    int i,j,k;
    double loglike=0.0, *ideal, prior=0.0;
    blockData *yeaData = (blockData *)yea_ptr;
    
    ideal = (double *) malloc((*yeaData).dim * sizeof(double));
    
    for(i=0; i<(*yeaData).length; i++){
        
        for(j=0; j<(*yeaData).dim; j++){
            ideal[j] = (*yeaData).ideal[(*(*yeaData).Data[i]).row -1 + j*(*yeaData).nrow];
        }
        loglike += l1LogLike((*(*yeaData).Data[i]).vote, ideal, par, (*yeaData).nayloc, (*yeaData).beta, (*yeaData).weight, (*yeaData).dim);
    }
    
    for(k=0;k<(*yeaData).dim;k++){
        prior += pow(par[k],2);
    }
    
    free(ideal);
    return(loglike - prior);
    
}

void readDataFromVector(int *inputVector, cell **data, int *nvotes, int *nrow, int *ncol) {
    int idx,i,j,nRow,nCol,nVotes;
    cell *p = (cell *) malloc((*nrow)*(*ncol)*sizeof(cell));
    
    nRow=*nrow;
    nCol=*ncol;
    Rprintf("Number of roll calls:\t%i\n",nCol);
    Rprintf("Number of legislators:\t%i\n",nRow);
    
    
    idx=0;
    nVotes = 0;
    for (i=0;i<nCol;i++) {
        for (j=0;j<nRow;j++) {
            if (inputVector[idx] != -1) {
                p[nVotes].row = j+1;
                p[nVotes].col = i+1;
                p[nVotes].vote = inputVector[idx];
                nVotes++;
            }
            idx++;
        }
    }
    *data = (cell *) realloc(p, nVotes * sizeof(cell));
    *nvotes = nVotes;
}


void formatData(cell *data, int nvotes, int nrow, int ncol, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData) {
    
    int i;
    int *rLengths, *cLengths, *rIndex, *cIndex;
    cell ***rData, ***cData;
    
    rLengths = calloc(nrow, sizeof(int));
    cLengths = calloc(ncol, sizeof(int));
    rIndex = calloc(nrow, sizeof(int));
    cIndex = calloc(ncol, sizeof(int));
    rData = (cell ***) malloc(nrow * sizeof(cell **));
    cData = (cell ***) malloc(ncol * sizeof(cell **));
    
    for(i=0; i<nvotes; i++) {
        rLengths[data[i].row - 1]++;
        cLengths[data[i].col - 1]++;
    }
    
    for(i=0; i<nrow; i++) rData[i] = (cell **) malloc(rLengths[i] * sizeof(cell *));
    for(i=0; i<ncol; i++) cData[i] = (cell **) malloc(cLengths[i] * sizeof(cell *));
    for(i=0; i<nvotes; i++) {
        rData[data[i].row-1][rIndex[data[i].row-1]] = &data[i];
        rIndex[data[i].row-1]++;
        cData[data[i].col-1][cIndex[data[i].col-1]] = &data[i];
        cIndex[data[i].col-1]++;
        
    }
    
    *rowLengths=rLengths;
    *colLengths=cLengths;
    *rowData=rData;
    *colData=cData;
    free(rIndex);
    free(cIndex);
}

void readInitsFromVector(double *inputIdeal,
                         double *inputBill,
                         double **idealpts, double **yealocs, double **naylocs,
                         double **beta, double **weight, int nrow, int ncol, int dims){
    
    int i;
    double *idealptr,*yeaptr,*nayptr,*betaptr,*weightptr;
    
    idealptr = (double *) malloc(dims * nrow * sizeof(double));
    yeaptr = (double *) malloc(dims * ncol * sizeof(double));
    nayptr = (double *) malloc(dims * ncol * sizeof(double));
    betaptr = (double *) malloc(sizeof(double));
    weightptr = (double *) malloc(sizeof(double));
    
    for(i=0;i<(nrow*dims);i++) idealptr[i]=inputIdeal[i];

    for(i=0;i<(ncol*dims);i++) {
        yeaptr[i] = inputBill[2*i];
        nayptr[i] = inputBill[2*i+1];
    }
    Rprintf("Number of dimensions:\t%i\n",dims);
    
    betaptr[0] = 10.0;
    if(dims==1) weightptr[0] = 1.0;
    else weightptr[0] = 0.5;
    
    
    *idealpts = idealptr;
    *yealocs = yeaptr;
    *naylocs = nayptr;
    *beta = betaptr;
    *weight = weightptr;
}



void writeDataOutput(double *output, int *op, int nrow, int ncol, double **idealpts, double **yealocs,
                     double **naylocs, double **beta, double **weight, int dims) {
    int i;
    int outputpos = *op;
    output[outputpos++] = (*beta)[0];
    for(i=0;i<(nrow*dims);i++) output[outputpos++] = (*idealpts)[i];
    for(i=0;i<(ncol*dims);i++) output[outputpos++] = (*yealocs)[i];
    for(i=0;i<(ncol*dims);i++) output[outputpos++] = (*naylocs)[i];
    output[outputpos++] = (*weight)[0];
    *op = outputpos;
}


void updateYea(cell ****colData,int **colLengths,double **idealpts,double **yealocs,double **naylocs,double **beta, double **weight, int nrow,int ncol,int dims){
    
    int i,j;
    double *yea,*nay,*par;
    blockData block;
    
    yea = (double *) malloc(dims * sizeof(double));
    nay = (double *) malloc(dims * sizeof(double));
    par = (double *) malloc(dims * sizeof(double));
    
    block.beta = (*beta)[0];
    block.weight = (*weight)[0];
    block.ideal = *idealpts;
    block.dim = dims;
    block.nrow = nrow;
    block.ncol = ncol;
    
    for(i=0;i<ncol;i++){
        
        for(j=0;j<dims;j++){
            yea[j] = (*yealocs)[i + j*ncol];
            nay[j] = (*naylocs)[i + j*ncol];
        }
        
        block.yealoc = yea;
        block.nayloc = nay;
        block.Data = (*colData)[i];
        block.length = (*colLengths)[i];
        par = slice(yeaLoglike,yea,&block,SLICE_W,dims);
        for(j=0;j<dims;j++){
            (*yealocs)[i + j*ncol] = par[j];
        }
    }
    
    free(yea);
    free(nay);
    free(par);
    
}


void updateNay(cell ****colData,int **colLengths, double **idealpts,double **yealocs,double **naylocs,double **beta,double **weight,int nrow,int ncol,int dims){
    
    int i,j;
    double *yea,*nay,*par;
    blockData block;
    
    yea = (double *) malloc(dims * sizeof(double));
    nay = (double *) malloc(dims * sizeof(double));
    par = (double *) malloc(dims * sizeof(double));
    
    block.beta = (*beta)[0];
    block.weight = (*weight)[0];
    block.ideal = *idealpts;
    block.dim = dims;
    block.nrow = nrow;
    block.ncol = ncol;
    
    for(i=0;i<ncol;i++){
        
        for(j=0;j<dims;j++){
            yea[j] = (*yealocs)[i + j*ncol];
            nay[j] = (*naylocs)[i + j*ncol];
        }
        
        block.yealoc = yea;
        block.nayloc = nay;
        block.Data = (*colData)[i];
        block.length = (*colLengths)[i];
        par = slice(nayLoglike,nay,&block,SLICE_W,dims);
        for(j=0;j<dims;j++){
            (*naylocs)[i + j*ncol] = par[j];
        }
    }
    
    free(yea);
    free(nay);
    free(par);
}


void updateLegis(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **weight,int nrow,int ncol,int dims){
    
    int i,j;
    double *ideal,*par;
    blockData block;
    
    ideal = (double *) malloc(dims * sizeof(double));
    par = (double *) malloc(dims * sizeof(double));
    
    block.beta = (*beta)[0];
    block.weight = (*weight)[0];
    block.yealoc = *yealocs;
    block.nayloc = *naylocs;
    block.dim = dims;
    block.nrow = nrow;
    block.ncol = ncol;
    
    for(i=0;i<nrow;i++){
        
        for(j=0;j<dims;j++){
            ideal[j] = (*idealpts)[i + j*nrow];
        }
        
        block.ideal = ideal;
        block.Data = (*rowData)[i];
        block.length = (*rowLengths)[i];
        par = slice(legisLoglike,ideal,&block,SLICE_W,dims);
        for(j=0;j<dims;j++){
            (*idealpts)[i + j*nrow] = par[j];
        }
    }
    
    free(ideal);
    free(par);
    
}

void sampleData(int Nsamples,int nrow,int ncol,int **rowLengths,int **colLengths,cell ****rowData,cell ****colData,double **idealpts,double **yealocs,double **naylocs, double **beta, double **weight, int verbose, int thin, int dims, double *output){
    
    int i;
    
    int outputpos = 0; // keep track of current position in the output container
    
    for(i=1;i<=Nsamples;i++){
        updateYea(colData,colLengths,idealpts,yealocs,naylocs,beta,weight,nrow,ncol,dims);
        updateNay(colData,colLengths,idealpts,yealocs,naylocs,beta,weight,nrow,ncol,dims);
        updateLegis(rowData,rowLengths,idealpts,yealocs,naylocs,beta,weight,nrow,ncol,dims);
        if(dims==2) updateWeight(rowData,rowLengths,idealpts,yealocs,naylocs,beta,weight,nrow,ncol,dims);
        updateBeta(rowData,rowLengths,idealpts,yealocs,naylocs,beta,weight,nrow,ncol,dims);
        if(i % thin==0) writeDataOutput(output,&outputpos,nrow,ncol,idealpts,yealocs,naylocs,beta,weight,dims);
        if(i % verbose==0) Rprintf("\tIteration %i of %i ...\n", i, Nsamples);
        R_FlushConsole();
    }
    
    
}

void freeData(int nrow, int ncol, cell **data, double **idealpts, double **yealocs, double **naylocs, int **rowLengths, int **colLengths, cell ****rowData, cell ****colData, double **beta, double **weight){
    
    int i;
    free(*data);
    free(*beta);
    free(*weight);
    free(*idealpts);
    free(*yealocs);
    free(*naylocs);
    for(i=0;i<nrow;i++) free((*rowData)[i]);
    for(i=0;i<ncol;i++) free((*colData)[i]);
    free(*rowData);
    free(*colData);
    free(*rowLengths);
    free(*colLengths);
    
}

double weightLegisLL(double *par, void *legis_ptr){
    int i,j;
    double loglike=0.0, *yea, *nay, *templegis;
    blockData *legisData = (blockData *)legis_ptr;
    
    yea = (double *) malloc((*legisData).dim * sizeof(double));
    nay = (double *) malloc((*legisData).dim * sizeof(double));
    templegis = (double *) malloc((*legisData).dim * sizeof(double));
    
    for(i=0; i<(*legisData).dim; i++) {
        templegis[i] = (*legisData).ideal[i];
    }
    
    for(i=0; i<(*legisData).length; i++) {
        
        for(j=0;j<(*legisData).dim;j++){
            yea[j] = (*legisData).yealoc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
            nay[j] = (*legisData).nayloc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
        }
        
        loglike += l1LogLike((*(*legisData).Data[i]).vote, templegis, yea, nay, (*legisData).beta, *par,(*legisData).dim);
    }
    
    free(yea);
    free(nay);
    free(templegis);
    return(loglike);
}
double weightLogLike(double *par, void *weight_ptr){
    int i, j;
    double loglike=0.0, *ideal,temppar;
    blockData block;
    
    temppar = *par;
    
    betaBlock *weightData = (betaBlock *)weight_ptr;
    block.weight = temppar;
    block.yealoc = (*weightData).yealoc;
    block.nayloc = (*weightData).nayloc;
    block.nrow = (*weightData).nrow;
    block.ncol = (*weightData).ncol;
    block.dim = (*weightData).dims;
    block.beta = (*weightData).beta;
    
    ideal = (double *) malloc((*weightData).dims * sizeof(double));
    
    for(i=0;i<block.nrow;i++){
        
        for(j=0;j<(*weightData).dims;j++){
            ideal[j] = (*weightData).ideal[i + j*block.nrow];
        }
        block.ideal = ideal;
        block.Data = (*weightData).rowData[i];
        block.length = (*weightData).rowLengths[i];
        
        loglike += weightLegisLL(&temppar, &block);
    }
    
    free(ideal);
    return(loglike);
    
}
void updateWeight(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **weight,int nrow,int ncol,int dims){
    
    double temp;
    betaBlock block;
    
    temp = (*weight)[0];
    
    block.beta = (*beta)[0];
    block.yealoc = *yealocs;
    block.nayloc = *naylocs;
    block.dims = dims;
    block.nrow = nrow;
    block.ncol = ncol;
    block.ideal = *idealpts;
    block.rowData = *rowData;
    block.rowLengths = *rowLengths;
    
    (*weight)[0] = slice_weight(weightLogLike,&temp,&block,SLICE_W,dims);
}


//Takes one legislator (pair of ideal pts), all roll call locations, gives log likelihood
//Essentially legisLoglike, but without a prior.  Used in betaLogLike
double betaLegisLL(double *par, void *legis_ptr){
    
    int i,j;
    double loglike=0.0, *yea, *nay, *templegis;
    blockData *legisData = (blockData *)legis_ptr;
    
    yea = (double *) malloc((*legisData).dim * sizeof(double));
    nay = (double *) malloc((*legisData).dim * sizeof(double));
    templegis = (double *) malloc((*legisData).dim * sizeof(double));
    
    for(i=0; i<(*legisData).dim; i++) {
        templegis[i] = (*legisData).ideal[i];
    }
    
    for(i=0; i<(*legisData).length; i++) {
        
        for(j=0;j<(*legisData).dim;j++){
            yea[j] = (*legisData).yealoc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
            nay[j] = (*legisData).nayloc[(*(*legisData).Data[i]).col -1 + j*(*legisData).ncol];
        }
        
        loglike += l1LogLike((*(*legisData).Data[i]).vote, templegis, yea, nay, *par, (*legisData).weight, (*legisData).dim);
    }
    
    free(yea);
    free(nay);
    free(templegis);
    return(loglike);
}


double betaLogLike(double *par, void *beta_ptr){
    
    int i, j;
    double loglike=0.0, *ideal,temppar;
    blockData block;
    
    temppar = *par;
    if(temppar < 0 ) return(-1e307);
    
    betaBlock *betaData = (betaBlock *)beta_ptr;
    block.beta = temppar;
    block.weight = (*betaData).weight;
    block.yealoc = (*betaData).yealoc;
    block.nayloc = (*betaData).nayloc;
    block.nrow = (*betaData).nrow;
    block.ncol = (*betaData).ncol;
    block.dim = (*betaData).dims;
    
    ideal = (double *) malloc((*betaData).dims * sizeof(double));
    
    for(i=0;i<block.nrow;i++){
        
        for(j=0;j<(*betaData).dims;j++){
            ideal[j] = (*betaData).ideal[i + j*block.nrow];
        }
        block.ideal = ideal;
        block.Data = (*betaData).rowData[i];
        block.length = (*betaData).rowLengths[i];
        
        loglike += betaLegisLL(&temppar, &block);
    }
    
    free(ideal);
    return(loglike);
    
}

void updateBeta(cell ****rowData,int **rowLengths,double **idealpts,double **yealocs,double **naylocs,double **beta,double **weight,int nrow,int ncol,int dims){
    
    double temp;
    betaBlock block;
    
    temp = (*beta)[0];
    
    block.weight = (*weight)[0];
    block.yealoc = *yealocs;
    block.nayloc = *naylocs;
    block.dims = dims;
    block.nrow = nrow;
    block.ncol = ncol;
    block.ideal = *idealpts;
    block.rowData = *rowData;
    block.rowLengths = *rowLengths;
    
    (*beta)[0] = slice_beta(betaLogLike,&temp,&block,SLICE_W,dims);
    
}

void Cl1ideal(int *dataVector, double *initIdeal, double *initBill,
              double *output, int *thin, int *ncol, int *nrow,
              int *Nsamples, int *dims, int *verbose) {
    cell *data;
    cell ***rowData, ***colData;
    int nvotes;
    int *rowLengths, *colLengths;
    double *idealpts, *yealocs, *naylocs, *beta, *weight;
    
    readDataFromVector(dataVector,&data,&nvotes,nrow,ncol);
    formatData(data,nvotes,*nrow,*ncol,&rowLengths,&colLengths,&rowData,&colData);
    readInitsFromVector(initIdeal,initBill,&idealpts,&yealocs,&naylocs,&beta,&weight,*nrow,*ncol,*dims);
    sampleData(*Nsamples,*nrow,*ncol,&rowLengths,&colLengths,&rowData,&colData,&idealpts,
               &yealocs,&naylocs,&beta,&weight,*verbose,*thin,*dims,output);
    freeData(*nrow,*ncol,&data,&idealpts,&yealocs,&naylocs,&rowLengths,&colLengths,&rowData,&colData,&beta,&weight);
}
