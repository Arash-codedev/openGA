// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <string>
#include <vector>
#include "openga.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TdLibrary/loss_function.h"
#include <omp.h>

struct MySolution
{
	std::vector<double> x;

	std::string to_string() const
	{
		std::ostringstream out;
		out<<"{";
		for(unsigned long i=0;i<x.size();i++)
			out<<(i?",":"")<<std::setprecision(10)<<x[i];
		out<<"}";
		return out.str();
	}
};
struct MyMiddleCost
{
	// This is where the results of simulation
	// is stored but not yet finalized.
	double cost;
};
typedef EA::Genetic<MySolution,MyMiddleCost> GA_Type;
typedef EA::GenerationType<MySolution,MyMiddleCost> Generation_Type;


//bool eval_solution(
//	const MySolution& p,
//	MyMiddleCost &c)
//{
//	constexpr double pi=3.141592653589793238;
//	c.cost=10*double(p.x.size());
//	for(unsigned long i=0;i<p.x.size();i++)
//		c.cost+=p.x[i]*p.x[i]-10.0*cos(2.0*pi*p.x[i]);
//	return true;
//}



MySolution crossover(
        const MySolution& X1,
        const MySolution& X2,
        const std::function<double(void)> &rnd01)
{
    MySolution X_new;
    for(unsigned long i=0;i<X1.x.size();i++)
    {
        double r=rnd01();
        X_new.x.push_back(r*X1.x[i]+(1.0-r)*X2.x[i]);
    }
    return X_new;
}
class GeneticAlgorithm{
public:
    bool compute(const MySolution& p,MyMiddleCost &c){
        constexpr double pi=3.141592653589793238;
        c.cost=10*double(p.x.size());
        for(unsigned long i=0;i<p.x.size();i++)
            c.cost+=p.x[i]*p.x[i]-10.0*cos(2.0*pi*p.x[i]);
        return true;
    }
    void init_genes(MySolution& p,const std::function<double(void)> &rnd01)
    {
        for(int i=0;i<5;i++)
            p.x.push_back(5.12*2.0*(rnd01()-0.5));
    }
    MySolution mutate(
            const Generation_Type& last_generation,
            const MySolution& X_base,
            const std::function<double(void)> &rnd01,
            double shrink_scale)
    {
//        MySolution X_best = last_generation.chromosomes[last_generation.best_chromosome_index].genes;
//        double fitness_best = last_generation.chromosomes[last_generation.best_chromosome_index].middle_costs.cost;
//        std::cout<<"fitness_best = "<<fitness_best<<std::endl;
        MySolution X_new;
        bool out_of_range;
        do{
            out_of_range=false;
            X_new=X_base;

            for(unsigned long i=0;i<X_new.x.size();i++)
            {
                double mu=1.7*rnd01()*shrink_scale;
                X_new.x[i]+=mu*(rnd01()-rnd01());
                if(std::abs(X_new.x[i])>5.12)
                    out_of_range=true;
            }
        } while(out_of_range);
        return X_new;
    }

};


double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
{
	// finalize the cost
//    td::slam::backend::CauchyLoss robust_function(5);
//    Eigen::Vector3d rho = Eigen::Vector3d::Zero();
//    double cost_inverse = 1 / X.middle_costs.cost;
//	robust_function.Compute(cost_inverse, rho) ;
//	cost_inverse = rho[0];
//	double cost_robust = 1 / cost_inverse;
//    std::cout<<"cost = "<<X.middle_costs.cost<<"  cost_robust = "<<cost_robust<<std::endl;
//    return cost_robust;
	return X.middle_costs.cost;
}

std::ofstream output_file;

void SO_report_generation(
	int generation_number,
	const EA::GenerationType<MySolution,MyMiddleCost> &last_generation,
	const MySolution& best_genes)
{
	std::cout
		<<"Generation ["<<generation_number<<"], "
		<<"Best="<<last_generation.best_total_cost<<", "
		<<"Average="<<last_generation.average_cost<<", "
		<<"Best genes=("<<best_genes.to_string()<<")"<<", "
		<<"Exe_time="<<last_generation.exe_time
		<<std::endl;

	output_file
		<<generation_number<<"\t"
		<<last_generation.average_cost<<"\t"
		<<last_generation.best_total_cost<<"\t"
		<<best_genes.x[0]<<"\t"
		<<best_genes.x[1]<<"\t"
		<<best_genes.x[2]<<"\t"
		<<best_genes.x[3]<<"\t"
		<<best_genes.x[4]<<"\t"
		<<"\n";
}

int main()
{
    EA::Chronometer timer_pipe;
    timer_pipe.tic();
    omp_set_num_threads(2);
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < 2; ++i) {
            std::string output_file_path = "../bin/result_so-rastrigin" + std::to_string(i) + ".txt";
            output_file.open(output_file_path,std::fstream::trunc | std::fstream::out);
            output_file
                    <<"step"<<"\t"
                    <<"cost_avg"<<"\t"
                    <<"cost_best"<<"\t"
                    <<"x_best0"<<"\t"
                    <<"x_best1"<<"\t"
                    <<"x_best2"<<"\t"
                    <<"x_best3"<<"\t"
                    <<"x_best4"
                    <<"\n";

            EA::Chronometer timer;
            timer.tic();

            MySolution x1,x2;
            x1.x.push_back(0.001);
            x1.x.push_back(0.001);
            x1.x.push_back(0.001);
            x1.x.push_back(0.001);
            x1.x.push_back(0.001);
            x2.x.push_back(0.2);
            x2.x.push_back(0.2);
            x2.x.push_back(0.1);
            x2.x.push_back(0.2);
            x2.x.push_back(0.1);
            std::vector<MySolution> init_genes_manually;
            init_genes_manually.push_back(x1);
            init_genes_manually.push_back(x2);
            GeneticAlgorithm computer;
            using std::placeholders::_1;
            using std::placeholders::_2;
            using std::placeholders::_3;
            using std::placeholders::_4;
            std::function<bool(const MySolution&,MyMiddleCost&)> eval_solution = std::bind( &GeneticAlgorithm::compute, computer, _1, _2);
            std::function<void(MySolution& p,const std::function<double(void)> &rnd01)> init_genes = std::bind( &GeneticAlgorithm::init_genes, computer, _1, _2);
            std::function<MySolution(const Generation_Type& last_generation,
                                     const MySolution& X_base,
                                     const std::function<double(void)> &rnd01,
                                     double shrink_scale)> mutate = std::bind( &GeneticAlgorithm::mutate, computer, _1, _2,_3,_4);
            GA_Type ga_obj;
            ga_obj.SetInitPopulationManually(init_genes_manually);
            ga_obj.problem_mode=EA::GA_MODE::SOGA;
            ga_obj.multi_threading=false;
            ga_obj.dynamic_threading= false;
            ga_obj.idle_delay_us=0; // switch between threads quickly
            ga_obj.verbose=false;
            ga_obj.population=10000;
            ga_obj.generation_max=1000;
            ga_obj.calculate_SO_total_fitness=calculate_SO_total_fitness;
            ga_obj.init_genes=init_genes;
            ga_obj.eval_solution=eval_solution;
            ga_obj.mutate=mutate;
            ga_obj.crossover=crossover;
            ga_obj.SO_report_generation=SO_report_generation;
            ga_obj.best_stall_max=20;
            ga_obj.average_stall_max=20;
            ga_obj.tol_stall_best=1e-6;
            ga_obj.tol_stall_average=1e-6;
            ga_obj.elite_count=10;
            ga_obj.crossover_fraction=0.7;
            ga_obj.mutation_rate=0.1;
            ga_obj.solve();

            std::cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<std::endl;

            output_file.close();
        }

    }
    std::cout<<"all problem is optimized in "<<timer_pipe.toc()<<" seconds."<<std::endl;
	return 0;
}
