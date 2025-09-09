#include "intracellularFBA.h"
// #include "buildModel.h"
#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <string>
#include <glpk.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <cstdarg>
#include <csetjmp>
#include <cstdio>
#include <cstdarg>

extern "C"
{
#include <matio.h>
    // Removed #include <matio_private.h> as it is not neededknm
}

using namespace std;

// static void my_glpk_error_hook(void *info, const char *msg);

int real_col;
bool model_initiated = false;

// Terminal output hook
static int my_term_hook(void *info, const char *s) {
    std::fputs("GLPK: ", stderr);
    std::fputs(s, stderr);
    std::fflush(stderr);
    return 0; // must return int
}

static void my_error_hook(void *info) {
    std::fprintf(stderr, "GLPK FATAL ERROR!\n");
    std::fflush(stderr);
    std::exit(EXIT_FAILURE);
}




class IntracellularFBA
{
private:
public:
    ~IntracellularFBA()
    {
        if (lp)
            glp_delete_prob(lp);
    }

    glp_prob *lp;
    glp_smcp parm;
    // Initialize GLPK parameters in the constructor
    IntracellularFBA() : lp(nullptr)
    {

        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_ERR;
        parm.meth = GLP_PRIMAL;

    // Install hooks once before solving
    glp_term_hook(my_term_hook, nullptr);
    glp_error_hook(my_error_hook, nullptr);

    }
    // IntracellularFBA(glp_prob *lp)
    // {
    //     this->lp = lp;
    //



    void buildModelFromMat(const char *file_path)
    {
        std::cout << "Reading MAT file...\n";
        mat_t *matfp = Mat_Open(file_path, MAT_ACC_RDONLY);
        if (!matfp)
        {
            std::cerr << "Error: Could not open MAT file: " << file_path << std::endl;
        }

        // Print MAT file version
        if (Mat_GetVersion(matfp) == MAT_FT_MAT4)
            std::cout << "MAT file version: v4\n";
        else if (Mat_GetVersion(matfp) == MAT_FT_MAT5)
            std::cout << "MAT file version: v5\n";
        else if (Mat_GetVersion(matfp) == MAT_FT_MAT73)
            std::cout << "MAT file version: v7.3\n";
        else
            std::cout << "MAT file version: unknown\n";
        // Mat_Close(matfp);
        // return nullptr;

        std::cout << "Reading variables from: " << file_path << "\n";
        matvar_t *var;
        var = Mat_VarReadNext(matfp); // Read the first variable (assumed struct)

        std::cout << "Variable name: " << var->name << "\n";

        // List all field names in the struct
        char *const *field_names = Mat_VarGetStructFieldnames(var);
        int num_fields = Mat_VarGetNumberOfFields(var);
        std::cout << "Number of fields: " << num_fields << "\n\n";
        for (int i = 0; i < num_fields; i++)
            std::cout << "Field " << i << ": " << field_names[i] << "\n";

        // std::cout << "Freeing variable: " << var->name << "...\n";

        matvar_t *S = Mat_VarGetStructFieldByName(var, "S", 0);
        matvar_t *lb = Mat_VarGetStructFieldByName(var, "lb", 0);
        matvar_t *ub = Mat_VarGetStructFieldByName(var, "ub", 0);
        matvar_t *rxn = Mat_VarGetStructFieldByName(var, "rxnNames", 0);
        // matvar_t *ver = Mat_VarGetStructFieldByName(var, "version", 0);
        if (!var)
        {
            Mat_Close(matfp);
            throw runtime_error("MAT variable read failed");
        }
        if (!S || !lb || !ub)
        {
            Mat_Close(matfp);
            throw runtime_error("Missing required fields S/lb/ub");
        }

        // char *version_str = (char *)ver->data;
        //  printf("Human-GEM version: %s\n", version_str);
        //  // Print model version to console
        //  if (ver && ver->data_type == MAT_T_DOUBLE && ver->data != nullptr)
        //  {
        //      double *ver_data = static_cast<double *>(ver->data);
        //      std::cout << "Model version: " << ver_data[0] << "\n";
        //  }
        //  else
        //  {
        //      std::cerr << "Version field not found or not of type double.\n";
        //  }

        
        // store reaction names
        std::vector<std::string> reaction_names;
        if (rxn->class_type == MAT_C_CELL)
        {
            size_t num_rxns = rxn->dims[0];
            matvar_t **rxn_cells = static_cast<matvar_t **>(rxn->data);
            for (size_t i = 0; i < num_rxns; ++i)
            {
                matvar_t *rxn_str = rxn_cells[i];
                if (rxn_str && (rxn_str->data_type == MAT_T_UTF8 || rxn_str->data_type == MAT_T_UINT8))
                {
                    std::string reaction_name(static_cast<char *>(rxn_str->data), rxn_str->nbytes);
                    reaction_names.push_back(reaction_name);
                    // std::cout << "Reaction " << i << ": " << reaction_name << "\n";
                }
                else if (rxn_str && rxn_str->data_type == MAT_T_DOUBLE)
                    reaction_names.push_back("[numeric data]");
                // std::cout << "Reaction " << i << ": [numeric data]" << "\n";
                else
                    reaction_names.push_back("[unknown type]");
                // std::cout << "Reaction " << i << ": [unknown type]" << "\n";
            }
        }
        else
            std::cout << "rxnNames is not a cell array.\n";

        // if (rxn->class_type == MAT_C_CELL)
        // {
        //     size_t num_rxns = rxn->dims[0];
        //     matvar_t **rxn_cells = static_cast<matvar_t **>(rxn->data);
        //     for (size_t i = 0; i < num_rxns; ++i)
        //     {
        //         matvar_t *rxn_str = rxn_cells[i];
        //         if (rxn_str && rxn_str->data_type == MAT_T_UTF8 || rxn_str->data_type == MAT_T_UINT8)
        //         {
        //             std::string reaction_name(static_cast<char *>(rxn_str->data), rxn_str->nbytes);
        //             std::cout << "Reaction " << i << ": " << reaction_name << "\n";
        //         }
        //         else if (rxn_str && rxn_str->data_type == MAT_T_DOUBLE)
        //         {
        //             std::cout << "Reaction " << i << ": [numeric data]" << "\n";
        //         }
        //         else
        //         {
        //             std::cout << "Reaction " << i << ": [unknown type]" << "\n";
        //         }
        //     }
        // }
        // else
        // {
        //     std::cout << "rxnNames is not a cell array.\n";
        // }

        // Build the stoichiometric matrix S_matrix from the sparse matrix S
        std::vector<std::vector<double>> S_matrix;

        // read matrix from mat file with matio and handle sparse and dense conditions
        if (S)
        {
            if (S->class_type == MAT_C_SPARSE)
            {
                mat_sparse_t *sp = static_cast<mat_sparse_t *>(S->data);
                size_t n_rows = S->dims[0];
                size_t n_cols = S->dims[1];
                size_t n_nz = sp->ndata;
                std::cout << "Sparse matrix detected.\n";
                std::cout << "\nDense dimensions of matrix: " << n_rows << " x " << n_cols << "\n";
                std::cout << "Number of non-zero elements: " << n_nz << "\n";
                std::cout << "1" << std::endl;
                double *values = static_cast<double *>(sp->data);
                std::cout << "2" << std::endl;
                mat_uint32_t *ir = sp->ir;
                mat_uint32_t *jc = sp->jc;
                std::cout << "3" << std::endl;
                S_matrix.clear(); // Clear any existing data in S_matrix
                std::cout << "444" << std::endl;
                // Convert sparse matrix to triplet format (row, col, value)
                for (size_t col = 0; col < n_cols; ++col)
                {
                    for (size_t idx = jc[col]; idx < jc[col + 1]; ++idx)
                        S_matrix.push_back({static_cast<double>(ir[idx]), static_cast<double>(col), values[idx]});
                }
                std::cout << "5" << std::endl;
            }
            else if ((S->class_type == MAT_C_DOUBLE || S->class_type == MAT_C_DOUBLE) && S->data != nullptr)
            {
                std::cout << "Dense matrix detected.\n";
                // handle dense matrix
                double *dense_data = static_cast<double *>(S->data);
                size_t n_rows = S->dims[0];
                size_t n_cols = S->dims[1];
                S_matrix.clear();
                for (size_t i = 0; i < n_rows; ++i)
                    for (size_t j = 0; j < n_cols; ++j)
                    {
                        double value = dense_data[i + j * n_rows]; // MATLAB column-major order
                        if (value != 0.0)
                            S_matrix.push_back({static_cast<double>(i), static_cast<double>(j), value});
                    }
            }
            else
                std::cout << "\nField is not a double or sparse matrix.\n";
        }
        else
        {
            std::cout << "Field 'S' not found in the variable.\n";
            // return 1;
        }
        std::cout << "6" << std::endl;
        // Determine the number of columns and rows in S_matrix
        int real_col = 0;
        for (const auto &col : S_matrix)
            if (!col.empty() && col[1] > real_col)
                real_col = static_cast<int>(col[1]);
        real_col += 1; // Adjust for 1-based indexing in GLPK
        std::cout << "7" << std::endl;
        int real_row = 0;
        for (const auto &col : S_matrix)
            if (!col.empty() && col[0] > real_row)
                real_row = static_cast<int>(col[0]);
        real_row += 1; // Adjust for 1-based indexing in GLPK
        std::cout << "8" << std::endl;
        // Extract lower and upper bounds for each reaction
        if (lb->dims[0] != real_col || ub->dims[0] != real_col)
        {
            std::cerr << "Lower/upper bounds array size does not match the number of columns in S_matrix.\n";
            std::cerr << "Expected: " << lb->dims[0] << ", but lb_size: " << ub->dims[0] << ", ub_size: " << real_col << "\n";
            // return 1;
        }
        double *lb_data = static_cast<double *>(lb->data);
        std::cout << "9" << std::endl;
        size_t lb_size = 1;
        if (lb && lb->data_type == MAT_T_DOUBLE && lb->data != nullptr)
            for (size_t i = 0; i < lb->rank; ++i)
                lb_size *= lb->dims[i];
        else
        {
            std::cerr << "Lower bounds not found or not of type double.\n";
            // return 1;
        }
        std::cout << "10" << std::endl;
        size_t ub_size = 1;
        if (ub && ub->data_type == MAT_T_DOUBLE && ub->data != nullptr)
            for (size_t i = 0; i < ub->rank; ++i)
                ub_size *= ub->dims[i];
        else
        {
            std::cerr << "Upper bounds not found or not of type double.\n";
            // return 1;
        }
        std::cout << "11" << std::endl;
        double *ub_data = static_cast<double *>(ub->data);
        std::cout << "12" << std::endl;
        Mat_Close(matfp);
        // Print some bounds for debugging
        // std::cout << "lb_size: " << lb_size << ", ub_size: " << ub_size << ", real_col: " << real_col << std::endl;
        // for (size_t i = 0; i < std::min(lb_size, size_t(10)); ++i)
        //     std::cout << "lb_data[" << i << "] = " << lb_data[i] << std::endl;
        // for (size_t i = 0; i < std::min(ub_size, size_t(10)); ++i)
        //     std::cout << "ub_data[" << i << "] = " << ub_data[i] << std::endl;

        // glp_term_hook(my_glpk_log, nullptr);
        // glp_error_hook(my_glpk_error, nullptr);

        // if (setjmp(GLPK_JMP) != 0)
        // {
        //     // We got here via longjmp from my_glpk_error
        //     // Do cleanup or just stop
        //     std::cerr << "Stopping after GLPK fatal error.\n"
        //               << std::flush;
        //     std::exit(EXIT_FAILURE);
        // }

        // glp_smcp parm;

        // ---- GLPK calls that may fail hard ----
        // Create GLPK problem for FBA

        // lp = glp_create_prob();
        //  ... set up rows/cols/matrix ...

        // ---------------------------------------



        std::cout << "13" << std::endl;

        this->lp = glp_create_prob();
        if (!this->lp)
            throw std::runtime_error("glp_create_prob failed");

        //glp_set_prob_name(lp, "FBA");

        // dimensions of the matrix
        std::cout << "S_matrix size: " << S_matrix.size() << " x " << S_matrix[0].size() << std::endl;
        std::cout << "14" << std::endl;
        std::cout << "Dense size: " << real_row << " x " << real_col << std::endl;
        glp_add_cols(lp, real_col); // +1 for the objective variable

        // use rxn names as glp col name
        std::cout << "15" << std::endl;
        for (int i = 1; i <= real_col; i++)
        {
            // glp_set_col_name(lp, i, ("v" + std::to_string(i)).c_str()); // Name each reaction variable
            glp_set_col_name(lp, i, reaction_names[i - 1].c_str()); // Name each reaction variable
            glp_set_col_kind(lp, i, GLP_CV);                        // Continuous variable
            glp_set_obj_coef(lp, i, 0);                             // Default objective coefficient
        }
        std::cout << "5" << std::endl;
        glp_set_obj_dir(lp, GLP_MAX);
        glp_set_obj_coef(lp, 25, 1); //  Cell Biomass

        // glp_set_obj_coef(lp, 5320, 1); // ATP
        //  glp_set_obj_coef(lp, 9281, 1);//Human Biomass
        //  glp_set_obj_coef(lp, 5320, 1);//ATP
        //   glp_set_obj_coef(lp, 18, -1);

        // read .mat matrix with matio
        // matvar_t *var2 = GetMatrixFromField(static_cast<const char *>("c:\\Users\\ozyur\\My Drive\\echan-metabolik-modeller\\Brain_models.mat"), "Result_ALL");
        // if (!var2)
        // {
        //     std::cerr << "Failed to read variable from MAT file.\n";
        //     //return 1;
        // }
        // // matvar_t *fluxes_matrix = Mat_VarGetStructFieldByName(var2, "fluxes", 0);

        // if (var2->data_type == MAT_T_DOUBLE && var2->data != nullptr && var2->rank == 2)
        // {
        //     double *data = static_cast<double *>(var2->data);
        //     size_t rows = var2->dims[0];
        //     size_t cols = var2->dims[1];

        //     // Print first column of the matrix (FBA fluxes for first condition)
        //     std::cout << "First column of FBA result matrix:" << std::endl;
        //     for (size_t i = 0; i < rows; ++i)
        //     {
        //         // std::cout << data[0 * rows + i] << " ";
        //         if (data[0 * rows + i] == 0) // Example condition
        //         {
        //             cout << i << " ";
        //             lb_data[i] = 0;
        //             ub_data[i] = 0;
        //         }
        //     }
        // }

        for (int i = 1; i <= real_col; i++)
        {
            // When lower bound equals upper bound, use GLP_FX instead of GLP_DB
            if (lb_data[i - 1] == ub_data[i - 1])
            {
                glp_set_col_bnds(lp, i, GLP_FX, lb_data[i - 1], lb_data[i - 1]);
            }
            else if (lb_data[i - 1] > ub_data[i - 1])
            {
                // If lower bound > upper bound, that's an error
                std::cerr << "Error: Column " << i << " has lb (" << lb_data[i - 1]
                          << ") > ub (" << ub_data[i - 1] << ")" << std::endl;
                // return 1;
            }
            else
            {
                // Normal case: lb < ub
                glp_set_col_bnds(lp, i, GLP_DB, lb_data[i - 1], ub_data[i - 1]);
            }
        }
        std::cout << "16" << std::endl;
        // glp_set_col_bnds(lp, 47, GLP_DB, -10, 1000);
        // glp_set_col_bnds(lp, 52, GLP_DB, -1, 0);
        // glp_set_col_bnds(lp, 55, GLP_DB, -1000, 1000);
        // glp_set_col_bnds(lp, 56, GLP_DB, -1000, 1000);
        // // glp_set_col_bnds(lp, 7735, GLP_FX, glutamate, glutamate);
        // glp_set_col_bnds(lp, 94, GLP_DB, -20, 0); // O2 difusion

        glp_add_rows(lp, real_row);
        std::cout << "5" << std::endl;
        for (int i = 1; i <= real_row; i++)
            glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
        int k = 0;

        // Declare arrays with appropriate sizes
        std::vector<int> ia(S_matrix.size() + 1, 0);
        std::vector<int> ja(S_matrix.size() + 1, 0);
        std::vector<double> ar(S_matrix.size() + 1, 0.0);

        for (size_t i = 0; i < S_matrix.size(); i++)
        {
            k++;
            ia[k] = static_cast<int>(S_matrix[i][0] + 1); // Convert to 1-based index
            ja[k] = static_cast<int>(S_matrix[i][1] + 1); // Convert to 1-based index
            ar[k] = S_matrix[i][2];
        }
        std::cout << "17" << std::endl;
        glp_load_matrix(lp, k, ia.data(), ja.data(), ar.data());
        std::cout << "18" << std::endl;

        auto start_opt = std::chrono::high_resolution_clock::now();
        glp_simplex(lp, &parm);
        auto end_opt = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_opt - start_opt;
        std::cout << "Optimization took " << elapsed.count() << " seconds.\n";
        // model_initiated = true;
        //  return lp; // Return the pointer as is
        //   Mat_VarFree(var);
        //   std::cout << "Variable freed.\n";
    }

    void changeBound(int col, double lb, double ub)
    {
        if (lb == ub)
            glp_set_col_bnds(lp, col, GLP_FX, lb, lb);
        else if (lb < ub)
            glp_set_col_bnds(lp, col, GLP_DB, lb, ub);
        else
            throw runtime_error("lb > ub in changeBound");
    }

    void optimizeModel()
    {
        auto start_opt = std::chrono::high_resolution_clock::now();

        // glp_init_smcp(&parm);
        // parm.msg_lev = GLP_MSG_ERR;
        // parm.meth = GLP_PRIMAL;

        // glp_iptcp iparm;
        // glp_init_iptcp(&iparm);
        // iparm.msg_lev = GLP_MSG_ERR;

        glp_simplex(lp, &parm);
        // glp_interior(lp, &iparm);
        auto end_opt = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed_opt = end_opt - start_opt;
        std::cout << "Model optimized in " << elapsed_opt.count() << " ms." << endl;
    }

    const std::vector<double> getFluxes() const
    {
        int real_col = glp_get_num_cols(lp);
        // print all flux values with their colnames
        for (int i = 1; i <= real_col; i++)
        {
            const char *reaction_name = glp_get_col_name(lp, i);
            std::cout << "flx_" << i << ": " << round(glp_get_col_prim(lp, i) * 100.0) / 100 << "\t" << reaction_name << endl;
        }

        std::vector<double> fluxes;
        for (int i = 1; i <= real_col; ++i)
            fluxes.push_back(glp_get_col_prim(lp, i));

        // glp_delete_prob(lp);
        std::cout << "GLPK problem deleted.\n";
        return fluxes;
    }
};