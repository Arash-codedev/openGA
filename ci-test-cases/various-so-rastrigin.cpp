// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <string>
#include <vector>
#include "openGA.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::vector;

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

void init_genes(MySolution& p,const std::function<double(void)> &rnd01)
{
	for(int i=0;i<5;i++)
		p.x.push_back(5.12*2.0*(rnd01()-0.5));
}

bool eval_solution(
	const MySolution& p,
	MyMiddleCost &c)
{
	constexpr double pi=3.141592653589793238;
	c.cost=10*double(p.x.size());
	for(unsigned long i=0;i<p.x.size();i++)
		c.cost+=p.x[i]*p.x[i]-10.0*cos(2.0*pi*p.x[i]);
	return true;
}

MySolution mutate(
	const MySolution& X_base,
	const std::function<double(void)> &rnd01,
	double shrink_scale)
{
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

double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
{
	// finalize the cost
	return X.middle_costs.cost;
}

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
}

struct TestResult
{
	double duration;
	string title;

	TestResult(double duration, string title):
		duration(duration),
		title(title)
	{
	}
};

vector<TestResult> test_results;

void run_test(
	bool multi_threading,
	bool dynamic_threading,
	int idle_delay_us,
	string title
	)
{
	cout<<"-------------------------------"<<endl;
	cout<<"Running the test: "<<title<<endl;

	EA::Chronometer timer;
	timer.tic();

	GA_Type ga_obj;
	ga_obj.problem_mode=EA::GA_MODE::SOGA;
	ga_obj.multi_threading=multi_threading;
	ga_obj.dynamic_threading=dynamic_threading;
	ga_obj.idle_delay_us=idle_delay_us;
	ga_obj.verbose=false;
	ga_obj.population=5000;
	ga_obj.generation_max=50;
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

	double duration = timer.toc();
	std::cout<<"The problem is optimized in "<<duration<<" seconds."<<endl;
	test_results.emplace_back(duration,title);
}

int main()
{
	run_test(false,false,0,"single-threaded");
	run_test(true,true,0,"multi-threaded with dynamic threading");
	run_test(true,false,0,"multi-threaded with static threading");
	run_test(true,true,100,"multi-threaded with dynamic threading with thread wait delay");

	cout<<endl;
	cout<<"Summary:"<<endl;
	for(auto &tr:test_results)
		cout<<tr.duration<<" seconds: "<<tr.title<<endl;
	cout<<endl;
	return 0;
}
