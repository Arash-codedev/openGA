// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <string>
#include <iostream>
#include <fstream>
#include "genetic.hpp"

struct MyGenes
{
	double x;
	double y;

	std::string to_string() const
	{
		return 
			"{x:"+std::to_string(x)+
			", y:"+std::to_string(y)+
			"}";
	}
};

struct MyMiddleCost
{
	// This is where the results of simulation
	// is stored but not yet finalized.
	double cost_A;
	double cost_B;
};

typedef EA::Genetic<MyGenes,MyMiddleCost> GA_Type;
typedef EA::GenerationType<MyGenes,MyMiddleCost> Generation_Type;

void init_genes(MyGenes& p,const std::function<double(void)> &rand)
{
	p.x=10.0*rand();
	p.y=10.0*rand();
}

bool eval_genes(
	const MyGenes& p,
	MyMiddleCost &c)
{
	double x=p.x;
	double y=p.y;
	// the middle comupations of cost:
	c.cost_A=log(1.0+x*sqrt(x*y));
	c.cost_B=98.0-100.0*(1.0-1.0/(1.0+y*sqrt(x*y)));
	return true; // genes are accepted
}

MyGenes mutate(
	const MyGenes& X_base,
	const std::function<double(void)> &rand,
	double shrink_scale)
{
	MyGenes X_new;
	double loca_scale=shrink_scale;
	if(rand()<0.4)
		loca_scale*=loca_scale;
	else if(rand()<0.1)
		loca_scale=1.0;
	bool in_range_x,in_range_y;
	double local_scale=shrink_scale;
	if(rand()<0.4)
		local_scale*=local_scale;
	else if(rand()<0.1)
		local_scale=1.0;
	do{
		X_new=X_base;
		X_new.x+=0.2*(rand()-rand())*local_scale;
		X_new.y+=0.2*(rand()-rand())*local_scale;
		in_range_x= (X_new.x>=0.0 && X_new.x<10.0);
		in_range_y= (X_new.y>=0.0 && X_new.y<10.0);
	} while(!in_range_x || !in_range_y);
	return X_new;
}

MyGenes crossover(
	const MyGenes& X1,
	const MyGenes& X2,
	const std::function<double(void)> &rand)
{
	MyGenes X_new;
	double r;
	r=rand();
	X_new.x=r*X1.x+(1.0-r)*X2.x;
	r=rand();
	X_new.y=r*X1.y+(1.0-r)*X2.y;
	return X_new;
}

std::vector<double> calculate_MO_objectives(const GA_Type::thisChromosomeType &X)
{
	return {
		X.middle_costs.cost_A,
		X.middle_costs.cost_B
	};
}

std::vector<double> distribution_objective_reductions(const std::vector<double> &objs)
{
	return objs;
}

void MO_report_generation(
	int generation_number,
	const EA::GenerationType<MyGenes,MyMiddleCost> &last_generation,
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
	output_file.open("./bin/result_mo1.txt");
	output_file<<"N"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"cost1"<<"\t"<<"cost2"<<"\n";
	std::vector<unsigned int> paretofront_indices=ga_obj.last_generation.fronts[0];
	for(unsigned int i:paretofront_indices)
	{
		const auto &X=ga_obj.last_generation.chromosomes[i];
		output_file
			<<i<<"\t"
			<<X.genes.x<<"\t"
			<<X.genes.y<<"\t"
			<<X.middle_costs.cost_A<<"\t"
			<<X.middle_costs.cost_B<<"\n";

	}
	output_file.close();
}

int main()
{
	EA::Chronometer timer;
	timer.tic();

	GA_Type ga_obj;
	ga_obj.problem_mode= EA::GA_MODE::NSGA_III;
	ga_obj.multi_threading=true;
	ga_obj.idle_delay_us=1; // switch between threads quickly
	ga_obj.verbose=false;
	ga_obj.population=40;
	ga_obj.generation_max=100;
	ga_obj.calculate_MO_objectives= calculate_MO_objectives;
	ga_obj.init_genes=init_genes;
	ga_obj.eval_genes=eval_genes;
	ga_obj.distribution_objective_reductions=distribution_objective_reductions;
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
