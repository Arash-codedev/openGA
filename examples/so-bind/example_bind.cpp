// This library is free and distributed under
// Mozilla Public License Version 2.0.

// This example is demonstrating how to use 
// genetic operations which belong to a class object.


#include <string>
#include <iostream>
#include <fstream>
#include <functional>
#include "openGA.hpp"

struct MySolution
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
	double cost_distance2;
	double cost_sqsin;
};

typedef EA::Genetic<MySolution,MyMiddleCost> GA_Type;
typedef EA::GenerationType<MySolution,MyMiddleCost> Generation_Type;

class MyClass
{
public:
	std::ofstream output_file;

	MyClass(std::string output_location)
	{
		output_file.open(output_location);
		output_file<<"step"<<"\t"<<"x_best"<<"\t"<<"y_best"<<"\t"<<"cost_avg"<<"\t"<<"cost_best"<<"\n";
	}

	~MyClass()
	{
		output_file.close();
	}

	void init_genes(MySolution& p,const std::function<double(void)> &rnd01)
	{
		p.x=20.0*rnd01()-10.0;
		p.y=20.0*rnd01()-10.0;
	}

	bool eval_solution(
		const MySolution& p,
		MyMiddleCost &c)
	{
		double x=p.x;
		double y=p.y;
		// see the surface plot at:
		// https://academo.org/demos/3d-surface-plotter/?expression=x*x%2By*y%2B30.0*sin(x*100.0*sin(y)%2By*100.0*cos(x))%2B125%2B45.0*sqrt(x%2By)*sin((15.0*(x%2By))%2F(x*x%2By*y))&xRange=-10%2C%2B10&yRange=-10%2C%2B10&resolution=100
		// 
		// the middle comupations of cost:
		if(x+y>0)
		{
			double predictable_noise=30.0*sin(x*100.0*sin(y)+y*100.0*cos(x));
			c.cost_distance2=x*x+y*y+predictable_noise;
			c.cost_sqsin=125+45.0*sqrt(x+y)*sin((15.0*(x+y))/(x*x+y*y));
			return true; // genes are accepted
		}
		else
			return false; // genes are rejected
	}

	MySolution mutate(
		const MySolution& X_base,
		const std::function<double(void)> &rnd01,
		double shrink_scale)
	{
		MySolution X_new;
		bool in_range_x,in_range_y;
		const double mu=0.2*shrink_scale; // mutation radius
		do{
			X_new=X_base;
			X_new.x+=mu*(rnd01()-rnd01());
			X_new.y+=mu*(rnd01()-rnd01());
			in_range_x= (X_new.x>=-10.0 && X_new.x<10.0);
			in_range_y= (X_new.y>=-10.0 && X_new.y<10.0);
		} while(!in_range_x || !in_range_y);
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
		X_new.x=r*X1.x+(1.0-r)*X2.x;
		r=rnd01();
		X_new.y=r*X1.y+(1.0-r)*X2.y;
		return X_new;
	}

	double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
	{
		// finalize the cost
		double cost1,cost2;
		cost1=X.middle_costs.cost_distance2;
		cost2=X.middle_costs.cost_sqsin;
		return cost1+cost2;
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

		output_file
			<<generation_number<<"\t"
			<<best_genes.x<<"\t"
			<<best_genes.y<<"\t"
			<<last_generation.average_cost<<"\t"
			<<last_generation.best_total_cost<<"\n";
	}
};

int main()
{
	MyClass myobject("./bin/result_so_bind.txt");

	EA::Chronometer timer;
	timer.tic();

	using std::bind;
	using std::placeholders::_1;
	using std::placeholders::_2;
	using std::placeholders::_3;

	GA_Type ga_obj;
	ga_obj.problem_mode= EA::GA_MODE::SOGA;
	ga_obj.multi_threading=true;
	ga_obj.idle_delay_us=1; // switch between threads quickly
	ga_obj.verbose=false;
	ga_obj.population=20;
	ga_obj.generation_max=1000;
	ga_obj.calculate_SO_total_fitness=bind(&MyClass::calculate_SO_total_fitness, &myobject, _1);
	ga_obj.init_genes= bind(&MyClass::init_genes, &myobject, _1, _2);
	ga_obj.eval_solution= bind(&MyClass::eval_solution, &myobject, _1, _2);
	ga_obj.mutate= bind(&MyClass::mutate, &myobject, _1, _2, _3);
	ga_obj.crossover= bind(&MyClass::crossover, &myobject, _1, _2, _3);
	ga_obj.SO_report_generation= bind(&MyClass::SO_report_generation, &myobject, _1, _2, _3);
	ga_obj.best_stall_max=10;
	ga_obj.elite_count=10;
	ga_obj.crossover_fraction=0.7;
	ga_obj.mutation_rate=0.4;
	ga_obj.solve();

	std::cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<std::endl;
	std::cout<<"This is a trivial example of openGA usage with std::bind."<<std::endl;

	return 0;
}
