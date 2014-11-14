/*Generic function type for use with anonymous functions*/
typedef void af_t(const double *, void *, double *);
typedef void af_free_t(void *);
typedef int af_init_t(void *, void *);
typedef void af_write_t(void *);
/*Basic leaf for creating products of functions.*/
typedef struct af_leaf
{
    af_t *f;
    void *p;
    struct af_leaf *next;
} af_leaf;


enum run_mode
{
    RM_CALCULATION = 1,
    RM_EXPERIMENTAL,
    RM_RUN
};

#define DERIV_MAX 2 //second order derivation is the most that we need.

/*Basis set implementation.*/
struct params_basis_set
{
    char name[256];//name of basis set
    char data_file[256];//plot data
    char params_file[256];//save paramaters after variation
    size_t maxn;
    size_t reserved_vars;//The number of vars of the basis set. Set at runtime?
    size_t n_pol;//terms of the polynomial expansion
    size_t n_exp;//terms of the exponential expansion
    char scaling_code;//for use in determing variable scaling schemes
    enum run_mode run_mode;//what are we doing
    double D;
    double b;//parameters of the basis set
    size_t deriv_accuracy;//points of integration, for derivative function
    double deriv_step;//distance between sampling points
    double *deriv_weights;//where we store the weights.

    double xmin;//box boundries
    double xmax;
    double epsrel;

    af_init_t *init;
    af_free_t *destroy;
};

//typedef struct options
//{
//    char c:1;
//    char r:1;
//    char z:1;
//} options;

/*Main data block.*/
struct params_variational_master
{
    struct params_basis_set *pbs;
    /*essential data values*/
//    size_t maxn; //Maximum excitation state.
    double **R_val; //Function overlaps used in the Grahm-Schmidt process. <f_i|e_j>
    double **masks; //Parameters to the basis set.
    double *S_val; //Overlap coefficients. As the diagonal is one, the diagonal is set here as the norm const.
    double *H_val; //Energy eigenvalues.
//    struct params_morse pm;


    /*auxiliaries*/
    size_t progress; //Count iterations within a calculation.
    struct params_variational_master_mask *pvmm;

    /*function objects*/
    af_leaf ***H; //Function to calculate energy expectant. <f_i|H|f_j>
    af_leaf ***S; //Function to calculate overlap. <f_i|f_j>

    /*simplify the syntax*/
    af_init_t *init;
    af_free_t *destroy;
};

typedef struct params_variational_master_mask
{
    size_t n;
    char error; //Catching special conditions.
    struct params_variational_master *pvm;
} pvmm;

struct params_variational_workspace
{
    nlopt_opt opt;
    double *lbounds;
    pvmm *p;
};


//Possible error codes:
#define NO_ERROR 0
#define ZERO_OVERLAP 1

int variate_pvmm(pvmm *p, struct params_variational_workspace *pvw);
int params_basis_set_init(struct params_basis_set *pbs, char *file);
double var_get(pvmm *p, size_t j);
double real_deriv_second(af_t *f, const double *x, pvmm *p);
double master(unsigned n, const double *x, double *grad, pvmm *p);
void print_data(struct params_variational_master *pvm);
void print_variational_params(pvmm *p);
void print_summary(pvmm *p);
void print_params(struct params_variational_master *pvm);
void read_params(struct params_variational_master *pvm);
void af_eval_integration_wrapper(unsigned ndim, const double *x, af_leaf *fdata, unsigned fdim, double *fval);
void af_free(af_leaf *tree);
void var_push(pvmm *p, size_t j, double x);
void params_variational_master_init(struct params_basis_set *pbs, struct params_variational_master *pvm);
void params_variational_master_free(struct params_variational_master *pvm);
void wave(const double *x, pvmm *p, double *ret);
void braket(struct params_basis_set *pbs, af_leaf *t, double *ret);
void ket_hamiltonian_wave(const double *x, pvmm *p, double *ret);
void morse_potential(const double *x, pvmm *p, double *ret);
void gs_learn(pvmm *p);
void gs_known(const double *x, pvmm *p, double *ret);
void method_experimental(struct params_variational_master *pvm);
void params_basis_set_free(struct params_basis_set *pbs);
void method_calculation(struct params_variational_master *pvm);
void method_rugid(void);
size_t var_getsize(pvmm *p);