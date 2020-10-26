// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <string>
#include "openGA.hpp"
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
