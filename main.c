#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

/*For minimization*/
#include <nlopt.h>
#include <cubature.h>

#include <c_fornberg.h>//for derivatives.

/*Header for this file*/
#include "magic.h"

/*Optimized weights for operations.*/
//#include "weights.h"

//This must be changed!
#define ZEROTOL 1e-5 //TODO: definitely change this
#define MIN_BOUND ZEROTOL
//TODO: arrays of timestamps and iterations to analyse.

/*Constants that should have no baring on behavior outside of this file.*/
#define NO_ERROR 0
#define TRANSLATION_SHIFT 0
#define GEOMETRIC_SHIFT 1
//scaling codes
#define SC_NONE '0'
#define SC_BOTH '1'
#define SC_ONLY_TRANS '2'
#define SC_ONLY_GEO '3'

#define INITIAL_VAR 1.00 //initial value of variables.

double safe_doubles[100];

//make life easy
struct params_basis_set global_pbs = {
    .init = (af_init_t *)&params_basis_set_init,
    .destroy = (af_free_t *)&params_basis_set_free
};

struct params_variational_master global_pvm = {
    .init = (af_init_t *)&params_variational_master_init,
    .destroy = (af_free_t *)&params_variational_master_free
};

//{{{ optimization
//determine the brent fornberg weights


//}}}
//{{{ printing
void
print_data(struct params_variational_master *pvm)
{
    size_t i;
    double x, V, psi;
    FILE *fp = fopen(pvm->pbs->data_file, "w");
    if(!fp)
    {
        fprintf(stderr, "Could not create file %s.\n", pvm->pbs->data_file);
        return;
    }

    size_t n;
    for(i=0;i<1000;i++)//TODO: make this configurable
    {
        x = -2.0 + i /100.0;
        morse_potential(&x, &pvm->pvmm[0], &V);
        fprintf(fp, "%g %g ", x, V);
        for(n=0;n<=pvm->pbs->maxn;n++)
        {
            gs_known(&x, &pvm->pvmm[n], &psi);
            fprintf(fp, "%g ", pow(psi,2.0) + (pvm->H_val[n] - pvm->H_val[0]));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return;
}

void
print_variational_params(pvmm *p)
{
    size_t j;

    braket(p->pvm->pbs, p->pvm->S[p->n][p->n], &p->pvm->S_val[p->n]);
    braket(p->pvm->pbs, p->pvm->H[p->n][p->n], &p->pvm->H_val[p->n]);
    printf("f(%.2zu)[%.4zu]=%.4E\t", p->n, p->pvm->progress, p->pvm->H_val[p->n]);
    for(j=0;j<var_nvars(p);j++)
        printf("%.4G\t", var_get_blind(p,j));
    printf("\n");

    return;
}

#define QUOTE(A) #A
#define STR(A) QUOTE(A)

void
print_summary(pvmm *p)
{
    /*Give an update on the state of the calculation.*/
    printf("////////Summary///////////////////////////\n");
    printf("f(%.2zu) overlap\t\t=\t%g\n"    , p->n, p->pvm->S_val[p->n]);
    printf("Iterations\t\t=\t%zu\n"         , p->pvm->progress);
    printf("Energy eigenvalue\t=\t%g\n"     , p->pvm->H_val[p->n]);
    printf("Normalization Constant\t=\t%g\n", 1.0/sqrt(fabs(p->pvm->R_val[p->n][p->n])));
    printf("Number of variables\t=\t%zu\n"  , var_nvars(p));
    printf("Minimization algorithm\t=\t%s\n", report_algorithm(p->pvm->pbs));
    printf("Scaling code\t\t=\t%c\n"        , p->pvm->pbs->scaling_code);
    printf("//////////////////////////////////////////\n");

    return;
}


void
print_params(struct params_variational_master *pvm)
{
    size_t i,j;
    FILE *fp = fopen(pvm->pbs->params_file, "w");
    if(!fp)
    {
        fprintf(stderr, "Could not create file %s.\n", pvm->pbs->params_file);
        return;
    }

    for(i=0; i<=pvm->pbs->maxn; i++)
    {//cycle through excitation states
        for(j=0; j<var_nvars(&pvm->pvmm[i]);j++)
            fprintf(fp, "%.20E\n", var_get_blind(&pvm->pvmm[i], j));
//        fprintf(fp, "\n");
    }

    fclose(fp);

    return;
}

//todo: return int
void
read_params(struct params_variational_master *pvm)
{  /* to be run after pbs is set up. Only run if the file exists and we are not calculating!
    */
    size_t i, j;
    double var;
    FILE *fp = fopen(pvm->pbs->params_file, "r");
    if(!fp)
    {
        fprintf(stderr, "Could not open file %s.\n", pvm->pbs->params_file);
        return;
    }


    for(i=0; i<=pvm->pbs->maxn; i++)
    {
        for(j=0;j<var_nvars(&pvm->pvmm[i]);j++)
            if(1 != fscanf(fp, "%lf", &var))
            {
                fprintf(stderr, "The params_file, %s, is not of the correct format.\n",
                    pvm->pbs->params_file);
                exit(1);/*The data file is corrupted. Allow user to delete it.*/
            }
        var_push_blind(&pvm->pvmm[i], j, var);
    }

    /*Restore normalization constants*/
    /*It is easier to do this than actually save them.*/
    for(i=0; i<=pvm->pbs->maxn; i++)
        gs_learn(&pvm->pvmm[i]);

}

//}}}
//{{{ anonymous functions

static inline af_leaf *
af_compile(size_t num, ...)
{
    af_leaf *t;
    af_leaf *r;
    size_t i;
    va_list args;
    if(!num)
        return NULL;
    t = malloc(sizeof(*t));
    r = t;
    va_start(args, num);
    for(i = 0; i < num; i++, t = t->next)
    {
        t->f = va_arg(args, af_t *);
        t->p = va_arg(args, void *);
        if(num == i + 1)
            t->next = NULL;
        else
            t->next = malloc(sizeof(*t));
    }
    va_end(args);
    return r;
}

static inline int
af_eval(const double *x, af_leaf *p, double *ret)
{
    double tmp = 0.0;
    double psum;
    size_t i;

    psum = 1.0;
    for(i=0;p; p = p->next,i++)
    {
        p->f(x, p->p, &tmp);
        psum *= tmp;
    }
    /*Assuming one dimensional function.*/

    *ret = psum;

    return 0;
}

void
af_eval_integration_wrapper(unsigned ndim, const double *x, af_leaf *fdata, unsigned fdim, double *fval)
{   //make cubature happy and keep the functions similar
    (void) ndim;
    (void) fdim;

    af_eval(x, fdata, fval);

    return;
}

void
af_free(af_leaf *tree)
{
    af_leaf *b;
    while(tree)
    {
        b = tree;
        tree = tree->next;
        free(b);
    }
}

//}}}
//{{{ utility functions
/*These functions are meant to clear up the code for wavefunctions.*/


/*We are hardcoding in the skew factors. This way parameters can be more easily read if part of the basis set is changed.*/

size_t
var_nvars(pvmm *p)
{   /*Return mask size to allocate.*/
    return (p->n/2) * 2.0 * p->pvm->pbs->var_growth +
        p->pvm->pbs->var_npol_base + p->pvm->pbs->var_nexp_base + 2;
}

double
var_get(pvmm *p, var_t type, size_t j)
{
    double ret;
    switch(type)
    {
        case EXP:
            ret = p->pvm->masks[p->n][p->n_pol+j+2];
            break;
        case POL:
            ret = p->pvm->masks[p->n][j+2];
            break;
        default:
            ret = p->pvm->masks[p->n][j];
    }
    if(!isnormal(ret)&&ret!=0)
    {
        fprintf(stderr,
            "var_get(p->n=%zu,%s,%zu) is attempting to return '%g.' "
            "This should never happen.\n", p->n, (type==EXP?"EXP":type==POL?"POL":"SPE"), j, ret);
        exit(2);
    }

    return ret;
}

void
var_push(pvmm *p, var_t type, size_t j, double x)
{
    switch(type)
    {
        case EXP:
            p->pvm->masks[p->n][p->n_pol+j+2] = x;
            return;
        case POL:
            p->pvm->masks[p->n][j+2] = x;
            return;
       default:
            p->pvm->masks[p->n][j] = x;
            return;
    }
}

/*Blind types are here, because I don't want to rewrite the file writing part. #YOLO*/
double
var_get_blind(pvmm *p, size_t j)
{
    return p->pvm->masks[p->n][j];
}

void
var_push_blind(pvmm *p, size_t j, double x)
{
    printf("var_push_blind()\n");
    p->pvm->masks[p->n][j] = x;
    return;
}


//}}}
//{{{ data managment functions

void
params_variational_master_init(struct params_basis_set *pbs, struct params_variational_master *pvm)
{
    size_t i, j;
    pvm->pbs       = pbs;

    pvm->progress  = 0;

    pvm->R_val     = malloc(sizeof(*pvm->R_val)*(1+pbs->maxn));
    pvm->masks     = malloc(sizeof(*pvm->masks)*(1+pbs->maxn));
    pvm->S_val     = calloc(sizeof(*pvm->S_val), (1+pbs->maxn));
    pvm->H_val     = calloc(sizeof(*pvm->H_val), (1+pbs->maxn));

    pvm->pvmm      = malloc(sizeof(*pvm->pvmm)*(1+pbs->maxn));
    pvm->H         = malloc(sizeof(*pvm->H)*(1+pbs->maxn));
    pvm->S         = malloc(sizeof(*pvm->S)*(1+pbs->maxn));

    /*initiate the pvm masks and create space for wave function parameters*/
    for(i=0;i<=pbs->maxn;i++)
    {
        pvm->pvmm[i].n     = i;
        pvm->pvmm[i].error = NO_ERROR;
        pvm->pvmm[i].pvm   = pvm;

        pvm->pvmm[i].n_pol = pvm->pbs->var_npol_base + pvm->pbs->var_growth/2;
        pvm->pvmm[i].n_exp = pvm->pbs->var_nexp_base + pvm->pbs->var_growth/2;
        pvm->pvmm[i].n_vars = 2 + pvm->pvmm[i].n_pol + pvm->pvmm[i].n_exp;

        pvm->R_val[i]      = calloc(sizeof(**pvm->R_val), (1+pbs->maxn));
        pvm->masks[i]      = malloc(sizeof(**pvm->masks) * var_nvars(&pvm->pvmm[i]));
        /*Let's initialize these so that we have some shape to our
         * functions to begin with. */
        for(j=2;j<=var_nvars(&pvm->pvmm[i]);j++)
            pvm->masks[i][j] = INITIAL_VAR;
        /*The default skews should have reasonable values.*/
        pvm->masks[i][TRANSLATION_SHIFT] = 0.0;
        pvm->masks[i][GEOMETRIC_SHIFT]   = 1.0;
    }

    /*initiate function pointers*/
    for(i=0;i<=pbs->maxn;i++)
    {
        pvm->H[i] = malloc(sizeof(**pvm->H)*(1+pbs->maxn));
        pvm->S[i] = malloc(sizeof(**pvm->S)*(1+pbs->maxn));
    }

    for(i=0;i<=pbs->maxn;i++)
        for(j=0;j<=pbs->maxn;j++)
        {
            pvm->H[i][j] = af_compile(2, &gs_known, &pvm->pvmm[i], &ket_hamiltonian_wave, &pvm->pvmm[j]);
            pvm->S[i][j] = af_compile(2, &gs_known, &pvm->pvmm[i], &gs_known, &pvm->pvmm[j]);
        };


    return;
}

void
params_variational_master_free(struct params_variational_master *pvm)
{
    size_t i, j;

    /*free function pointers*/
    for(i=0;i<=pvm->pbs->maxn;i++)
    {
        for(j=0;j<=pvm->pbs->maxn;j++)
        {
            af_free(pvm->H[i][j]);
            af_free(pvm->S[i][j]);
        }
        free(pvm->R_val[i]);
        free(pvm->masks[i]);
        free(pvm->H[i]);
        free(pvm->S[i]);
    }

    free(pvm->H);
    free(pvm->S);

    /*free the arrays*/
    free(pvm->R_val);
    free(pvm->masks);
    free(pvm->S_val);
    free(pvm->H_val);
    free(pvm->pvmm );

//    free(pvm);

    return;
}

//}}}
//{{{ wave functions

double
real_deriv_second(af_t *f, const double *x, pvmm *p)
{   /*Calculate a linear, equidistant derivative of order 2. //'M'*/
    size_t i;
    double tmp, ret = 0.0, y[1];

    for(i=0;i<p->pvm->pbs->deriv_accuracy;i++)
    {
        y[0] = x[0] + i * p->pvm->pbs->deriv_step + p->pvm->pbs->deriv_step/2.0
             - p->pvm->pbs->deriv_accuracy * p->pvm->pbs->deriv_step/2.0;
        f(y, p, &tmp);
        ret += p->pvm->pbs->deriv_weights[i] * tmp;
    }

    return ret;
}


//  void
//  wave(const double *x, pvmm *p, double *ret)
//  { /*The basis set implementation.*/
//      /*Based of the basis sets shown in Koch, Schuck and Wacker*/
//      double y = *x;
//      size_t i;
//      double poly, expr = 0.0;
//
//      /*shift axis if necessary*/
//      switch(p->pvm->pbs->scaling_code)
//      {
//          case SC_NONE:
//              break;
//          case SC_BOTH:
//              y -= var_get(p, TRANSLATION_SHIFT);
//              y *= var_get(p, GEOMETRIC_SHIFT);
//              break;
//          case SC_ONLY_TRANS:
//              y -= var_get(p, TRANSLATION_SHIFT);
//              break;
//          case SC_ONLY_GEO:
//              y *= var_get(p, GEOMETRIC_SHIFT);
//              break;
//          default:
//              fprintf(stderr, "The scaling code is invalid. Aborting.\n");
//              exit(1);
//      }
//
//      /*Set polynomial first term*/
//      if(!(p->n%2))
//          poly = 1.0;
//      else
//          poly = y;
//
//
//      /*calculate polynomial term*/
//      for(i=1;i<=p->pvm->pbs->n_pol;i++)
//      {
//          double A = pow(-1.0, i+1.0);
//          double k = var_get(p, i+1);
//          double po = pow(y, 2.0 * i + (p->n%2));
//
//          poly += pow(-1.0, i+1.0) * var_get(p, i+1) * pow(y, 2.0 * i + (p->n%2));
//      }
//      /*calculate exponential term*/
//      for(i=1;i<=p->pvm->pbs->n_exp;i++)
//      {
//          expr -= var_get(p, i+1+p->pvm->pbs->n_pol) * pow(y, 2.0 * i);
//  //        printf("var_get =%g\n", var_get(p, i+1+p->pvm->pbs->n_pol));
//      }
//
//      *ret = poly * exp(expr);
//
//      return;
//  }

void
wave(const double *x, pvmm *p, double *ret)
{
    double y, poly, expr = 0.0;
    size_t i;

     /*Let us shift the axis.*/
     y = *x;
     switch(p->pvm->pbs->scaling_code)
     {
         case SC_NONE:
            break;
         case SC_BOTH:
             y -= var_get(p, SPE, TRANSLATION_SHIFT);
             y *= var_get(p, SPE, GEOMETRIC_SHIFT);
             break;
         case SC_ONLY_TRANS:
             y -= var_get(p, SPE, TRANSLATION_SHIFT);
             break;
           case SC_ONLY_GEO:
            y *= var_get(p, SPE, GEOMETRIC_SHIFT);
            break;
        default:
            fprintf(stderr, "The scaling code is invalid. Aborting.\n");
            exit(1);
    }
    /**/
    /*Next, initiate the polynomial and exponential terms.*/
    if(!(p->n%2))//Is even?
        poly = 1.0;
    else
        poly = y;
    for(i=0; i<p->n_pol; i++) //TODO: clean this up
        poly += pow(-1.0, i+2.0) * var_get(p, POL, i) * pow(y, (2.0 * (i+1.0)+(p->n%2)));
    for(i=0; i<p->n_exp; i++)
        expr -= var_get(p, EXP, i) * pow(y, (i+1.0)*2.0);


    *ret = poly * exp(expr);
    /**/

    return;
}

void
braket(struct params_basis_set *pbs, af_leaf *t, double *ret)
{
    double err;

    pcubature(1, (integrand)&af_eval_integration_wrapper, t, 1, &pbs->xmin,
        &pbs->xmax, 0, 0, pbs->epsrel, ERROR_INDIVIDUAL, ret, &err);

//    printf("braket = %g\n", *ret);
    //TODO: don't ignore 'err'
    return;
}

void
ket_hamiltonian_wave(const double *x, pvmm *p, double *ret)
{
    double V;
    double T;
    double tmp;

    morse_potential(x, p, &V);
    gs_known(x, p, &tmp);
    V *= tmp;
    T = -real_deriv_second((af_t *)&gs_known, x, p)/2.0;
    T *= 1.0/(2.0);//input mass
    *ret = T + V;

    return;
}

//potential energy waves:

void
morse_potential(const double *x, pvmm *p, double *ret)
{   /*The potential energy function. TODO: add variety to this?*/
    *ret = p->pvm->pbs->D * pow(1.0 - exp(- p->pvm->pbs->b * *x), 2.0);

    return;
}


//Grahm-Schmidt orthogonalization functions:
void
gs_learn(pvmm *p)
{  /* For use with a wave function that is freshly adjusted.
    * The function should be normal and we orthogonalize here.
    * Also, the diagonal of the coefficients' matrix is the normalization factors, because otherwise
    * it would be a diaganol of ones.
    * As a note if (p->pvmm->n == 0) we return V_k, this makes life easy and saves us a function.
    *
    * Calculate:
    * U_k = V_k - S [j=1;j<n-1;j++] {<V_k, e_j> e_j}
    * Write:
    * <V_k,e_j> under p->pvm->R_val[k,j]
    * Note only writing a portion of the values to avoid repeats!
    */
    size_t i;
    double s;

    p->pvm->R_val[p->n][p->n] = 1.0; //Set this to 1.0 to give the un-normalized values of wave().
    for(i=0; i<=p->n; i++)
    {
        braket(p->pvm->pbs, p->pvm->S[p->n][i], &s);
//        p->pvm->R_val[p->n][i] = p->pvm->S[p->n][i]();
        if(!isnormal(s))
        { //catch boundary conditions.
            if(!(p->n == i)&&(0 == s))
            {
                fprintf(stderr, "Tried to orthogonalize orthogonal functions, "
                                "this could probably be avoided (optimize me!).\n");
                continue; //it's ok (orthogonal functions)
            }

            /*I realize that the overlap might not be zero, but this 'just werks.'*/
            fprintf(stderr, "The overlap is %g. Returning ZERO_OVERLAP.\n", s);
            p->error = ZERO_OVERLAP;
        }
        p->pvm->R_val[p->n][i] = s;
    }

    return;
}

void
gs_known(const double *x, pvmm *p, double *ret)
{ /*An implementation of the Grahm-Schmidt orthogonalization process.*/
    size_t i;
    double tmp;

    wave(x, p, ret);
    for(i=0;i<p->n;i++)
    {
        gs_known(x, &p->pvm->pvmm[i], &tmp);
        *ret -= p->pvm->R_val[p->n][i] * tmp; //NOTE: we are only writing the upper triangle.
    }
    *ret /= sqrt(fabs(p->pvm->R_val[p->n][p->n]));

    return;
}
//}}}
//{{{ routines

double
master(unsigned n, const double *x, double *grad, pvmm *p)
{//calculate the energy of a selected excitation state, also the wrapper for variational calculations.
    (void) n;
    (void) grad;
    (void) x;//we are smarter than the library. #YOLO

    double ret;
    size_t i;

    p->pvm->progress++;//keep count of how many iterations we make; TODO: make an array to keep check of each state!
    gs_learn(p);//orthonormalize the new variables.
    braket(p->pvm->pbs, p->pvm->H[p->n][p->n], &ret);

    for(i=0;i<var_nvars(p);i++)
        safe_doubles[i] = var_get_blind(p, i);

    if(ZERO_OVERLAP == p->error)
    { //Handle single point constraints.
        ret = HUGE_VAL;//this just werks. It's a minimization algorithm.
        p->error = NO_ERROR;
    }
    print_variational_params(p);


    if(p->pvm->progress==98)
        print_summary(p);

    return ret;
}




static void
params_variational_workspace_init(struct params_variational_workspace *pvw, pvmm *p)
{
    pvw->opt = nlopt_create(p->pvm->pbs->algorithm, var_nvars(p));
    nlopt_set_min_objective(pvw->opt, (nlopt_func)&master, p);
    nlopt_set_ftol_rel(pvw->opt, p->pvm->pbs->epsrel);
    pvw->lbounds = malloc(sizeof(*pvw->lbounds) * var_nvars(p));

    size_t i;
    for(i=0;i<var_nvars(p);i++)
        pvw->lbounds[i] = MIN_BOUND;
    pvw->lbounds[0] = -HUGE_VAL;//scale_i //TODO: this is dirty, make it good.
    nlopt_set_lower_bounds(pvw->opt, pvw->lbounds);

    pvw->p = p;

    return;
}

static void
params_variational_workspace_free(struct params_variational_workspace *pvw)
{
    free(pvw->lbounds);
    nlopt_destroy(pvw->opt);

    return;
}

int
variate_pvmm(struct params_variational_workspace *pvw)
{ /*Minimize a single energy level*/
    pvw->p->pvm->progress++;

    nlopt_optimize(pvw->opt, pvw->p->pvm->masks[pvw->p->n], &pvw->p->pvm->H_val[pvw->p->n]);
    printf("Convergence met for p->n==%zu.\n", pvw->p->n);

    size_t i;
    for(i = 0; i<var_nvars(pvw->p);i++)
        var_push_blind(pvw->p, i, safe_doubles[i]);


    return 0;
}


//}}}
//{{{ ui
#define USAGE(prog) "Usage: " prog " [Option] [File]\n\
Calculate a vibrational basis set and run calculations thereon.\n\
Options:\n\
    c\t\t\tcalculate a basis set.\n\
    r\t\t\trun a calculation.\n\
    z\t\t\trun an experimental method (hardcoded and unstable).\n\
\n\
File:\n\
    target file must be in the proper format (whatever that means).\n\
\n\
Exit Status:\n\
    0\t\t\tsuccess.\n\
    1\t\t\tprogram was aborted early because of a value handed to it.\n\
    2\t\t\tan internal error occured.\n"


//TODO: move this to header
struct algorithm_dict
{
    char *name;
    nlopt_algorithm algorithm;
};

struct algorithm_dict algorithm_dict[] =
{
    { "NLOPT_GN_DIRECT_L" , NLOPT_GN_DIRECT_L },
    { "NLOPT_GN_ORIG_DIRECT_L", NLOPT_GN_DIRECT_L },
    { "NLOPT_LN_NELDERMEAD", NLOPT_LN_NELDERMEAD },
    { "NLOPT_LN_SBPLX", NLOPT_LN_SBPLX },
    { "NLOPT_LN_PRAXIS", NLOPT_LN_PRAXIS },
    { "NLOPT_GN_CRS2_LM", NLOPT_GN_CRS2_LM },
    { "NLOPT_GN_ISRES", NLOPT_GN_ISRES },
    { "NLOPT_LN_COBYLA", NLOPT_LN_COBYLA },
    { "NLOPT_LN_NEWUOA", NLOPT_LN_NEWUOA },
    { "NLOPT_LN_BOBYQA", NLOPT_LN_BOBYQA },
};

int
translate_algorithm(char *name, struct params_basis_set *pbs)
{//Take the string version of the algorithm and
 //set the struct params_basis_set algorithm field

    size_t i;
    for(i=0; algorithm_dict[i].algorithm ; i++)
    {
        if(!strcmp(algorithm_dict[i].name, name))
        {
            pbs->algorithm = algorithm_dict[i].algorithm;
            return 0;
        }
    }

    fprintf(stderr, "Minimization algorithm %s not found.\n", name);
    return 1;
}

char *
report_algorithm(struct params_basis_set *pbs)
{ /*Return a pointer to the heap algoithm_dict string corresponding to the algorithm.*/
    size_t i;
    for(i=0; algorithm_dict[i].algorithm ; i++)
    {
        if(algorithm_dict[i].algorithm == pbs->algorithm)
        {
            return algorithm_dict[i].name;
        }
    }

    return NULL;
}

int
params_basis_set_init(struct params_basis_set *pbs, char *file)
{
    FILE *fp = fopen(file, "r");
    char buf[256]; //For translating nlopt_algorithm
    if(!fp)
        return 1;

    /*The order of parameters:
        name of basis set:
            name\n
        excitation states:
            maxn\n
        minimization algotihm:
            algorithm\n
        morse parameters:
            D b\n
        basis set parameters:
            reserved_vars, n_pol, n_exp, scaling_code\n
        output files:
            data, params\n
    */


    if(1 != fscanf(fp,"%s\n", pbs->name))
        goto error;
    if(1 != fscanf(fp,"%zu\n", &pbs->maxn))
        goto error;
    if(1 != fscanf(fp, "%s\n", buf))
        goto error;
    if(2 != fscanf(fp,"%lf %lf\n", &pbs->D, &pbs->b))
        goto error;
    if(2 != fscanf(fp,"%s %s\n", pbs->data_file, pbs->params_file))
        goto error;
    if(2 != fscanf(fp, "%zu %lf\n", &pbs->deriv_accuracy, &pbs->deriv_step))
        goto error;
    if(3 != fscanf(fp, "%lf %lf %lf\n", &pbs->xmin, &pbs->xmax, &pbs->epsrel))
        goto error;
    if(4 != fscanf(fp, "%zu %zu %zu %c", &pbs->var_growth, &pbs->var_npol_base,
            &pbs->var_nexp_base, &pbs->scaling_code))
        goto error;

    fclose(fp);


    translate_algorithm(buf, pbs);

    /*setup weights*/
    /*the abscissa*/
    size_t i;
    double *x = malloc(sizeof(*x) * pbs->deriv_accuracy);
    for(i = 0; i < pbs->deriv_accuracy; i++)
        x[i] = -(pbs->deriv_accuracy*pbs->deriv_step)/2 + i * pbs->deriv_step + pbs->deriv_step/2.0;

    /*the actual weights*/
    pbs->deriv_weights = malloc(sizeof(*pbs->deriv_weights)*pbs->deriv_accuracy*(DERIV_MAX+1));
    fornberg_populate_weights(0.0, x, pbs->deriv_accuracy - 1, DERIV_MAX, pbs->deriv_weights);

    free(x);

    return 0;

error:
    fprintf(stderr,"Invalid option file!\n");
    fclose(fp);

    return 1;
}

void
params_basis_set_free(struct params_basis_set *pbs)
{
    free(pbs->deriv_weights);

    return;
}

static inline int
input_parser(int argc, char **argv, struct params_basis_set *pbs)
{
    if(2 >= argc)
    {
        fprintf(stderr, USAGE("%s"), argv[0]);
        return 1;
    }
    /*set run mode*/
    switch(argv[1][0])
    {
        case 'c':
            pbs->run_mode = RM_CALCULATION;
            break;
        case 'r':
            pbs->run_mode = RM_RUN;
            break;
        case 'z':
            pbs->run_mode = RM_EXPERIMENTAL;
            break;
        default:
            fprintf(stderr, USAGE("%s"), argv[0]);
            return 1;
    }

    if(pbs->init(pbs, argv[2]))
        return 1;

    return 0;
}

void
method_calculation(struct params_variational_master *pvm)
{
    struct params_variational_workspace pvw;

    size_t i, j;
    for(i=0; i<=pvm->pbs->maxn; i++)
    {
        params_variational_workspace_init(&pvw, &pvm->pvmm[i]);
        variate_pvmm(&pvw);

        /*TODO: check this*/
        print_summary(&pvm->pvmm[i]);
        print_data(pvm);
        print_params(pvm);
        params_variational_workspace_free(&pvw);
    }

    return;
}

void
method_rugid(void)
{
    /*Defaults*/
//    double D, step_D;
//    double b, step_b;
//    double a;
    /**/

    return;
}

//}}}

double safe_doubles[100];


#include<stdio.h>
#include<signal.h>
#include<unistd.h>

void sig_handler(int signo)
{
    if (signo == SIGUSR1)
        printf("received 'SIGUSR1,' shutting down.\n");
    else if (signo == SIGINT)
        printf("Received 'SIGINT,' shutting down.\n");


    fprintf(stderr, "Caught signal and aborting run.\nA second signal will abort unconditionally.\n");
    /*save data and exit*/
    print_data(&global_pvm);
    print_params(&global_pvm);

    exit(0);
    return;
}

void
sig_init(void)
{
    if (signal(SIGUSR1, sig_handler) == SIG_ERR)
        fprintf(stderr, "Cannot catch 'SIGUSR1.'\n");
    if (signal(SIGINT, sig_handler) == SIG_ERR)
        fprintf(stderr, "Cannot catch 'SIGINT.'\n");
    // A long long wait so that we can easily issue a signal to this process

    return;
}


void
test_1(const double *x, void *kek, double *ret)
{
    *ret = 5.0* *x;
}

void
test_2(const double *x, void *kek, double *ret)
{
    *ret = 2.5* *x;
}

void
method_kek(struct params_basis_set *pbs)
{
    af_leaf *t = af_compile(2, &test_1, NULL, &test_2, NULL);
    size_t i;
    double x[1];
    double ret;

    for(i=0;i<10;i++)
    {
        x[0] = i/10.0;
        af_eval(x, t, &ret);
        printf("%g %g\n", *x, ret);
    }

    braket(pbs, t, &ret);
    printf("integral = %g\n", ret);

    af_free(t);

    return;
}

int
main(int argc, char **argv)
{
    extern struct params_basis_set global_pbs;
    extern struct params_variational_master global_pvm;

    if(input_parser(argc, argv, &global_pbs))
        return 1;
    global_pvm.init(&global_pbs, &global_pvm);

    printf("Startup successful. Entering run mode '%c'.\n", argv[1][0]);
    printf("File locations:\n\tData file:\t%s\n\tParameter file:\t%s\n",
        global_pbs.data_file,
        global_pbs.params_file
    );
    /*do stuff here*/



    sig_init(); //gracefully hadle signals
    print_summary(&global_pvm.pvmm[0]);
    sig_init(); //gracefully hadle signals
    switch(global_pbs.run_mode)
    {
        case RM_CALCULATION:
            method_calculation(&global_pvm);
            break;
        case RM_RUN:
            break;
        case RM_EXPERIMENTAL:
            method_kek(&global_pbs);
//            variational_basis_set(global_pvm.pvmm);
//            method_experimental(&global_pvm);
//            method_rugid();
            break;
    }


    global_pvm.destroy(&global_pvm);
    global_pbs.destroy(&global_pbs);

    return 0;
}
