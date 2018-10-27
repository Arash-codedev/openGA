// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <string>
#include "openga.hpp"
#include "gui.hpp"
#include <fstream>

struct MySolution
{
	double R,G,B;

	std::string to_string() const
	{
		const unsigned red = (unsigned)R;
		const unsigned green = (unsigned)G;
		const unsigned blue = (unsigned)B;
		char hexstr[16];
		snprintf(hexstr,sizeof(hexstr),"%02x%02x%02x",red,green,blue);
		std::string retstr=hexstr;
		return retstr;
	}
};

struct MyMiddleCost
{
	double R,G,B;
	double cost_user_score;
};

typedef EA::Genetic<MySolution,MyMiddleCost> GA_Type;
typedef EA::GenerationType<MySolution,MyMiddleCost> Generation_Type;

void init_genes(MySolution& p,const std::function<double(void)> &rnd01)
{
	p.R=255.0*rnd01();
	p.G=255.0*rnd01();
	p.B=255.0*rnd01();
}

bool eval_solution_IGA(
	const MySolution& p,
	MyMiddleCost &c,
	const EA::GenerationType<MySolution,MyMiddleCost>&)
{
	c.R=p.R;
	c.G=p.G;
	c.B=p.B;
	return true; // genes are accepted
}

MySolution mutate(
	const MySolution& X_base,
	const std::function<double(void)> &rnd01,
	double shrink_scale)
{
	MySolution X_new;
	(void) shrink_scale;
	bool in_range_R,in_range_G,in_range_B;
	do{
		X_new=X_base;
		X_new.R+=100*(rnd01()-rnd01());
		X_new.G+=100*(rnd01()-rnd01());
		X_new.B+=100*(rnd01()-rnd01());
		in_range_R= (X_new.R>=0.0 && X_new.R<255.0);
		in_range_G= (X_new.G>=0.0 && X_new.G<255.0);
		in_range_B= (X_new.B>=0.0 && X_new.B<255.0);
	} while(!in_range_R || !in_range_G || !in_range_B);
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
	X_new.R=r*X1.R+(1.0-r)*X2.R;
	r=rnd01();
	X_new.G=r*X1.G+(1.0-r)*X2.G;
	r=rnd01();
	X_new.B=r*X1.B+(1.0-r)*X2.B;
	return X_new;
}

void calculate_IGA_total_fitness(GA_Type::thisGenerationType &g)
{


	for(unsigned int i=0;i<g.chromosomes.size();i++)
	{
		GA_Type::thisChromosomeType &X=g.chromosomes[i];
		// X.total_cost=100.0-X.middle_costs.cost_user_score;
		gui_subject_R=X.middle_costs.R;
		gui_subject_G=X.middle_costs.G;
		gui_subject_B=X.middle_costs.B;
		refresh_gui();
		refresh_gui();
		std::cout<<"How much do you like this ("<<X.genes.to_string()<<") blue color (0-100%)? ";
		std::cin>>X.middle_costs.cost_user_score;
		X.total_cost=100.0-X.middle_costs.cost_user_score;
	}
}

std::ofstream output_file;

void SO_report_generation(
	int generation_number,
	const EA::GenerationType<MySolution,MyMiddleCost> &last_generation,
	const MySolution& best_genes)
{
	std::cout
		<<"Generation ["<<generation_number<<"], "
		<<"Best="<<100.0-last_generation.best_total_cost<<", "
		<<"Average="<<100.0-last_generation.average_cost<<", "
		<<"Best genes=("<<best_genes.to_string()<<")"<<", "
		<<"Exe_time="<<last_generation.exe_time
		<<std::endl;

	output_file
		<<generation_number<<"\t"
		<<best_genes.to_string()<<"\t"
		<<100.0-last_generation.average_cost<<"\t"
		<<100.0-last_generation.best_total_cost<<"\n";
}

int main()
{
	output_file.open("./bin/result_iga1.txt");
	output_file<<"step"<<"\t"<<"color_best"<<"\t"<<"cost_avg"<<"\t"<<"cost_best"<<"\n";
	init_gui();

	GA_Type ga_obj;
	ga_obj.problem_mode= EA::GA_MODE::IGA;
	ga_obj.verbose=false;
	ga_obj.population=15;
	ga_obj.generation_max=20;
	ga_obj.calculate_IGA_total_fitness=calculate_IGA_total_fitness;
	ga_obj.init_genes= init_genes;
	ga_obj.eval_solution_IGA= eval_solution_IGA;
	ga_obj.mutate= mutate;
	ga_obj.crossover= crossover;
	ga_obj.SO_report_generation= SO_report_generation;
	ga_obj.elite_count=3;
	double non_elit_fraction=1-double(ga_obj.elite_count)/double(ga_obj.population);
	ga_obj.crossover_fraction=non_elit_fraction;
	ga_obj.mutation_rate=0.1;
	ga_obj.solve();

	output_file.close();
	return 0;
}
