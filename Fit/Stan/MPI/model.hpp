// Code generated by Stan version 2.18.1

#include <stan/model/model_header.hpp>

namespace KATRIN_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model.stan");
    reader.add_event(72, 70, "end", "model.stan");
    return reader;
}

template <typename T0__>
std::vector<typename boost::math::tools::promote_args<T0__>::type>
signal(const std::vector<T0__>& pars, std::ostream* pstream__);

double
bkg(std::ostream* pstream__);

std::vector<int>
GetData(std::ostream* pstream__);

int
GetSubrunNum(std::ostream* pstream__);

class KATRIN : public prob_grad {
private:
    int nsubrun;
    vector<int> count;
public:
    KATRIN(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    KATRIN(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;

        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "KATRIN_namespace::KATRIN";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {

            // validate, data variables
            // initialize data variables
            current_statement_begin__ = 15;
            nsubrun = int(0);
            stan::math::fill(nsubrun, std::numeric_limits<int>::min());
            stan::math::assign(nsubrun,GetSubrunNum(pstream__));
            current_statement_begin__ = 16;
            validate_non_negative_index("count", "nsubrun", nsubrun);
            count = std::vector<int>(nsubrun,int(0));
            stan::math::fill(count, std::numeric_limits<int>::min());
            stan::math::assign(count,GetData(pstream__));


            // validate transformed data
            current_statement_begin__ = 15;
            current_statement_begin__ = 16;

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 22;
            ++num_params_r__;
            current_statement_begin__ = 23;
            ++num_params_r__;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~KATRIN() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("A_sig")))
            throw std::runtime_error("variable A_sig missing");
        vals_r__ = context__.vals_r("A_sig");
        pos__ = 0U;
        context__.validate_dims("initialization", "A_sig", "double", context__.to_vec());
        double A_sig(0);
        A_sig = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0,1000000.0,A_sig);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable A_sig: ") + e.what());
        }

        if (!(context__.contains_r("A_bkg")))
            throw std::runtime_error("variable A_bkg missing");
        vals_r__ = context__.vals_r("A_bkg");
        pos__ = 0U;
        context__.validate_dims("initialization", "A_bkg", "double", context__.to_vec());
        double A_bkg(0);
        A_bkg = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0,10,A_bkg);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable A_bkg: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            local_scalar_t__ A_sig;
            (void) A_sig;  // dummy to suppress unused var warning
            if (jacobian__)
                A_sig = in__.scalar_lub_constrain(0,1000000.0,lp__);
            else
                A_sig = in__.scalar_lub_constrain(0,1000000.0);

            local_scalar_t__ A_bkg;
            (void) A_bkg;  // dummy to suppress unused var warning
            if (jacobian__)
                A_bkg = in__.scalar_lub_constrain(0,10,lp__);
            else
                A_bkg = in__.scalar_lub_constrain(0,10);


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body
            {
            current_statement_begin__ = 36;
            local_scalar_t__ mass;
            (void) mass;  // dummy to suppress unused var warning

            stan::math::initialize(mass, DUMMY_VAR__);
            stan::math::fill(mass,DUMMY_VAR__);
            stan::math::assign(mass,0);
            current_statement_begin__ = 37;
            local_scalar_t__ endpoint;
            (void) endpoint;  // dummy to suppress unused var warning

            stan::math::initialize(endpoint, DUMMY_VAR__);
            stan::math::fill(endpoint,DUMMY_VAR__);
            stan::math::assign(endpoint,18574);
            current_statement_begin__ = 38;
            local_scalar_t__ e1;
            (void) e1;  // dummy to suppress unused var warning

            stan::math::initialize(e1, DUMMY_VAR__);
            stan::math::fill(e1,DUMMY_VAR__);
            stan::math::assign(e1,12.6);
            current_statement_begin__ = 39;
            local_scalar_t__ B_A;
            (void) B_A;  // dummy to suppress unused var warning

            stan::math::initialize(B_A, DUMMY_VAR__);
            stan::math::fill(B_A,DUMMY_VAR__);
            stan::math::assign(B_A,0.00026800000000000001);
            current_statement_begin__ = 40;
            local_scalar_t__ B_S;
            (void) B_S;  // dummy to suppress unused var warning

            stan::math::initialize(B_S, DUMMY_VAR__);
            stan::math::fill(B_S,DUMMY_VAR__);
            stan::math::assign(B_S,2.52);
            current_statement_begin__ = 41;
            local_scalar_t__ B_max;
            (void) B_max;  // dummy to suppress unused var warning

            stan::math::initialize(B_max, DUMMY_VAR__);
            stan::math::fill(B_max,DUMMY_VAR__);
            stan::math::assign(B_max,4.2000000000000002);
            current_statement_begin__ = 42;
            local_scalar_t__ A1;
            (void) A1;  // dummy to suppress unused var warning

            stan::math::initialize(A1, DUMMY_VAR__);
            stan::math::fill(A1,DUMMY_VAR__);
            stan::math::assign(A1,0.20399999999999999);
            current_statement_begin__ = 43;
            local_scalar_t__ A2;
            (void) A2;  // dummy to suppress unused var warning

            stan::math::initialize(A2, DUMMY_VAR__);
            stan::math::fill(A2,DUMMY_VAR__);
            stan::math::assign(A2,0.055599999999999997);
            current_statement_begin__ = 44;
            local_scalar_t__ w1;
            (void) w1;  // dummy to suppress unused var warning

            stan::math::initialize(w1, DUMMY_VAR__);
            stan::math::fill(w1,DUMMY_VAR__);
            stan::math::assign(w1,1.8500000000000001);
            current_statement_begin__ = 45;
            local_scalar_t__ w2;
            (void) w2;  // dummy to suppress unused var warning

            stan::math::initialize(w2, DUMMY_VAR__);
            stan::math::fill(w2,DUMMY_VAR__);
            stan::math::assign(w2,12.5);
            current_statement_begin__ = 46;
            local_scalar_t__ e2;
            (void) e2;  // dummy to suppress unused var warning

            stan::math::initialize(e2, DUMMY_VAR__);
            stan::math::fill(e2,DUMMY_VAR__);
            stan::math::assign(e2,14.300000000000001);
            current_statement_begin__ = 47;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning

            stan::math::initialize(sigma, DUMMY_VAR__);
            stan::math::fill(sigma,DUMMY_VAR__);
            stan::math::assign(sigma,3.3999999999999999);
            current_statement_begin__ = 48;
            validate_non_negative_index("pars", "12", 12);
            vector<local_scalar_t__> pars(12);
            stan::math::initialize(pars, DUMMY_VAR__);
            stan::math::fill(pars,DUMMY_VAR__);
            stan::math::assign(pars,static_cast<std::vector<local_scalar_t__> >(stan::math::array_builder<local_scalar_t__ >().add(mass).add(endpoint).add(B_A).add(B_S).add(B_max).add(A1).add(A2).add(w1).add(w2).add(e1).add(e2).add((sigma * 1.0000000000000001e-18)).array()));
            current_statement_begin__ = 49;
            validate_non_negative_index("sig", "nsubrun", nsubrun);
            vector<local_scalar_t__> sig(nsubrun);
            stan::math::initialize(sig, DUMMY_VAR__);
            stan::math::fill(sig,DUMMY_VAR__);
            stan::math::assign(sig,signal(pars, pstream__));
            current_statement_begin__ = 50;
            validate_non_negative_index("pred", "nsubrun", nsubrun);
            vector<local_scalar_t__> pred(nsubrun);
            stan::math::initialize(pred, DUMMY_VAR__);
            stan::math::fill(pred,DUMMY_VAR__);


            current_statement_begin__ = 51;
            if (pstream__) {
                stan_print(pstream__,"nsubrun: ");
                stan_print(pstream__,nsubrun);
                *pstream__ << std::endl;
            }
            current_statement_begin__ = 53;
            for (int n = 1; n <= nsubrun; ++n) {

                current_statement_begin__ = 54;
                stan::model::assign(pred, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            ((get_base1(sig,n,"sig",1) * A_sig) + (A_bkg * bkg(pstream__))), 
                            "assigning variable pred");
                current_statement_begin__ = 57;
                if (as_bool(logical_eq(n,110))) {
                    current_statement_begin__ = 57;
                    if (pstream__) {
                        stan_print(pstream__,"Nbin: ");
                        stan_print(pstream__,n);
                        stan_print(pstream__," Pred: ");
                        stan_print(pstream__,get_base1(pred,n,"pred",1));
                        stan_print(pstream__," count: ");
                        stan_print(pstream__,get_base1(count,n,"count",1));
                        stan_print(pstream__," A_sig: ");
                        stan_print(pstream__,A_sig);
                        *pstream__ << std::endl;
                    }
                }
            }
            current_statement_begin__ = 59;
            lp_accum__.add(poisson_log<propto__>(count, pred));
            }

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("A_sig");
        names__.push_back("A_bkg");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "KATRIN_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double A_sig = in__.scalar_lub_constrain(0,1000000.0);
        double A_bkg = in__.scalar_lub_constrain(0,10);
        vars__.push_back(A_sig);
        vars__.push_back(A_bkg);

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters
            if (include_tparams__) {
            }
            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "KATRIN";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "A_sig";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "A_bkg";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "A_sig";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "A_bkg";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
    }

}; // model

}

typedef KATRIN_namespace::KATRIN stan_model;

