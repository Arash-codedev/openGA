// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <string>
#include <iostream>
#include <fstream>
#include "openga.hpp"

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

	std::cout<<"Generation ["<<generation_number<<"], ";
	std::cout<<"Pareto-Front {";
	for(unsigned int i=0;i<pareto_front.size();i++)
	{
		std::cout<<(i>0?",":"");
		std::cout<<pareto_front[i];
	}
	std::cout<<"}"<<std::endl;
}

void save_results(const GA_Type &ga_obj)
{
	std::ofstream output_file;
	output_file.open("./bin/result_mo2.txt");
	output_file<<"N"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"cost1"<<"\t"<<"cost2"<<"\n";
	std::vector<unsigned int> paretofront_indices=ga_obj.last_generation.fronts[0];
	for(unsigned int i:paretofront_indices)
	{
		const auto &X=ga_obj.last_generation.chromosomes[i];
		output_file
			<<i<<"\t"
			<<X.genes.x1<<"\t"
			<<X.genes.x2<<"\t"
			<<X.middle_costs.cost_f1<<"\t"
			<<X.middle_costs.cost_f2<<"\n";

	}
	output_file.close();
}

int main()
{
	EA::Chronometer timer;
	timer.tic();

	GA_Type ga_obj;
	ga_obj.problem_mode= EA::GA_MODE::NSGA_III;
	ga_obj.dynamic_threading=false;
	ga_obj.multi_threading=true;
	ga_obj.idle_delay_us=1; // switch between threads quickly
	ga_obj.verbose=false;
	ga_obj.population=40;
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

	std::cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<std::endl;

	save_results(ga_obj);
	return 0;
}
