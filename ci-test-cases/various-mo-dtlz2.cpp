// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "openGA.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;

struct MySolution
{
	double x1;
	double x2;

	std::string to_string() const
	{
		return 
			"{x1:"+std::to_string(x1)+
			", x2:"+std::to_string(x2)+
			"}";
	}
};

struct MyMiddleCost
{
	// This is where the results of simulation
	// is stored but not yet finalized.
	double cost_f1;
	double cost_f2;
};

typedef EA::Genetic<MySolution,MyMiddleCost> GA_Type;
typedef EA::GenerationType<MySolution,MyMiddleCost> Generation_Type;

void init_genes(MySolution& p,const std::function<double(void)> &rnd01)
{
	p.x1=1.0*rnd01();
	p.x2=1.0*rnd01();
}

bool eval_solution(
	const MySolution& p,
	MyMiddleCost &c)
{
	double x1=p.x1;
	double x2=p.x2;
	// the middle comupations of cost:
	double g=(x1-0.5)*(x1-0.5)+(x2-0.5)*(x2-0.5);
	constexpr double pi=3.14159265359;
	c.cost_f1=(1+g)*cos(x1*pi/2);
	c.cost_f2=(1+g)*sin(x1*pi/2);
	return true; // genes are accepted
}

MySolution mutate(
	const MySolution& X_base,
	const std::function<double(void)> &rnd01,
	double shrink_scale)
{
	MySolution X_new;
	bool in_range_x1,in_range_x2;
	do{
		X_new=X_base;
		X_new.x1+=0.1*(rnd01()-rnd01())*shrink_scale;
		X_new.x2+=0.1*(rnd01()-rnd01())*shrink_scale;
		in_range_x1= (X_new.x1>=0.0 && X_new.x1<1.0);
		in_range_x2= (X_new.x2>=0.0 && X_new.x2<1.0);
	} while(!in_range_x1 || !in_range_x2);
	return X_new;
}

MySolution crossover(
	const MySolution& X1,
	const MySolution& X2,
	const std::function<double(void)> &rnd01)
{
	MySolution X_new;
	double r;
	r=rnd01();
	X_new.x1=r*X1.x1+(1.0-r)*X2.x1;
	r=rnd01();
	X_new.x2=r*X1.x2+(1.0-r)*X2.x2;
	return X_new;
}

std::vector<double> calculate_MO_objectives(const GA_Type::thisChromosomeType &X)
{
	return {
		X.middle_costs.cost_f1,
		X.middle_costs.cost_f2
	};
}

void MO_report_generation(
	int generation_number,
	const EA::GenerationType<MySolution,MyMiddleCost> &last_generation,
	const std::vector<unsigned int>& pareto_front)
{
	(void) last_generation;

	cout<<"Generation ["<<generation_number<<"], ";
	cout<<"Pareto-Front {";
	for(unsigned int i=0;i<pareto_front.size();i++)
	{
		cout<<(i>0?",":"");
		cout<<pareto_front[i];
	}
	cout<<"}"<<endl;
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
	ga_obj.problem_mode= EA::GA_MODE::NSGA_III;
	ga_obj.dynamic_threading=dynamic_threading;
	ga_obj.multi_threading=multi_threading;
	ga_obj.idle_delay_us=idle_delay_us;
	ga_obj.verbose=false;
	ga_obj.population=100;
	ga_obj.generation_max=100;
	ga_obj.calculate_MO_objectives= calculate_MO_objectives;
	ga_obj.init_genes=init_genes;
	ga_obj.eval_solution=eval_solution;
	ga_obj.mutate=mutate;
	ga_obj.crossover=crossover;
	ga_obj.MO_report_generation=MO_report_generation;
	ga_obj.crossover_fraction=0.7;
	ga_obj.mutation_rate=0.4;
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
