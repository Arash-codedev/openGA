// This library is free and distributed under
// Mozilla Public License Version 2.0.

#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <armadillo>
#include <thread>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include <ctime>
#include <string>
#include <iostream>
#include <stdexcept>

#ifndef NS_EA_BEGIN
#define NS_EA_BEGIN namespace EA {
#define NS_EA_END	 }
#endif

NS_EA_BEGIN;

enum class GA_MODE
{
	SOGA,
	IGA,
	NSGA_III
};

template<typename GeneType,typename MiddleCostType>
struct ChromosomeType
{
	GeneType genes;
	MiddleCostType middle_costs; 	// individual costs
	double total_cost;				// for single objective
	arma::vec objectives;			// for multi-objective
};

template<typename GeneType,typename MiddleCostType>
struct GenerationType
{
	std::vector<ChromosomeType<GeneType,MiddleCostType>> chromosomes;
	double best_total_cost= (std::numeric_limits<double>::infinity()); // for single objective
	double average_cost= 0.0; // for single objective

	int best_chromosome_index=-1; // for single objective
	std::vector<int> sorted_indices; // for single objective
	std::vector<std::vector<uint>> fronts; // for multi-objective
	std::vector<double> selection_chance_cumulative;
	double exe_time;
};

template<typename GeneType,typename MiddleCostType>
struct GenerationType_SO_abstract
{
	double best_total_cost= (std::numeric_limits<double>::infinity()); // for single objective
	double average_cost= 0.0;// for single objective

	GenerationType_SO_abstract(const GenerationType<GeneType,MiddleCostType> &generation):
		best_total_cost(generation.best_total_cost),
		average_cost(generation.average_cost)
	{
	}
};

enum class StopReason
{
	Undefined,
	MaxGenerations,
	StallAverage,
	StallBest,
	UserRequest
};

class Chronometer
{
protected:
	timespec time_start, time_stop;
	bool initialized;
public:

	Chronometer() : 
			initialized(false)
	{
	}

	void tic()
	{
		initialized=true;
		clock_gettime(CLOCK_MONOTONIC, &time_start);
	}

	double toc()
	{
		if(!initialized)
			throw std::runtime_error("Chronometer is not initialized!");
		clock_gettime(CLOCK_MONOTONIC, &time_stop);
		long diff_sec=time_stop.tv_sec-time_start.tv_sec;
		long diff_nsec=time_stop.tv_nsec-time_start.tv_nsec;
		double time_diff=static_cast<double>(diff_sec)+(1e-9)*static_cast<double>(diff_nsec);
		return time_diff;
	}

};

template<typename GeneType,typename MiddleCostType>
class Genetic
{
private:
	std::mt19937_64 rng; // random generator
	std::uniform_real_distribution<double> unif_dist;
	int average_stall_count;
	int best_stall_count;
	arma::vec ideal_objectives;		// for multi-objective
	arma::mat extreme_objectives;	// for multi-objective
	arma::vec scalarized_objectives_min;	// for multi-objective
	std::vector<arma::vec> reference_vectors;
	double shrink_scale;

public:

	typedef ChromosomeType<GeneType,MiddleCostType> thisChromosomeType;
	typedef GenerationType<GeneType,MiddleCostType> thisGenerationType;
	typedef GenerationType_SO_abstract<GeneType,MiddleCostType> thisGenSOAbs;

	////////////////////////////////////////////////////

	GA_MODE problem_mode;
	uint population;
	double crossover_fraction;
	double mutation_rate;
	bool verbose;
	int generation_step;
	int elite_count;
	int generation_max;
	double tol_stall_average;
	int average_stall_max;
	double tol_stall_best;
	int best_stall_max;
	uint reference_vector_divisions;
	bool enable_reference_vectors;
	bool multi_threading;
	bool dynamic_threading;
	int N_threads;
	bool user_request_stop;
	long idle_delay_us;

	std::function<void(thisGenerationType&)> calculate_IGA_total_fitness;
	std::function<double(const thisChromosomeType&)> calculate_SO_total_fitness;
	std::function<arma::vec(thisChromosomeType&)> calculate_MO_objectives;
	std::function<arma::vec(const arma::vec&)> distribution_objective_reductions;
	std::function<void(GeneType&,const std::function<double(void)> &rand)> init_genes;
	std::function<bool(const GeneType&,MiddleCostType&)> eval_genes;
	std::function<bool(const GeneType&,MiddleCostType&,const thisGenerationType&)> eval_genes_IGA;
	std::function<GeneType(const GeneType&,const std::function<double(void)> &rand,double shrink_scale)> mutate;
	std::function<GeneType(const GeneType&,const GeneType&,const std::function<double(void)> &rand)> crossover;
	std::function<void(int,const thisGenerationType&,const GeneType&)> SO_report_generation;
	std::function<void(int,const thisGenerationType&,const std::vector<uint>&)> MO_report_generation;
	std::function<void(void)> custom_refresh;
	std::function<double(int)> set_shrink_scale=[](int n){return (n<=5?1.0:1.0/sqrt(n-5+1));};

	std::vector<thisGenSOAbs> generations_so_abs;
	thisGenerationType last_generation;

	////////////////////////////////////////////////////

	Genetic() :
		unif_dist(0.0,1.0),
		problem_mode(GA_MODE::SOGA),
		population(50),
		crossover_fraction(0.7),
		mutation_rate(0.1),
		verbose(false),
		generation_step(-1),
		elite_count(5),
		generation_max(100),
		tol_stall_average(1e-4),
		average_stall_max(10),
		tol_stall_best(1e-6),
		best_stall_max(10),
		reference_vector_divisions(10),
		enable_reference_vectors(true),
		multi_threading(true),
		dynamic_threading(true),
		N_threads(std::thread::hardware_concurrency()),
		user_request_stop(false),
		idle_delay_us(1000),
		calculate_IGA_total_fitness(nullptr),
		calculate_SO_total_fitness(nullptr),
		calculate_MO_objectives(nullptr),
		distribution_objective_reductions(nullptr),
		init_genes(nullptr),
		eval_genes(nullptr),
		eval_genes_IGA(nullptr),
		mutate(nullptr),
		crossover(nullptr),
		SO_report_generation(nullptr),
		MO_report_generation(nullptr),
		custom_refresh(nullptr)
	{
		// initialize the random number generator with time-dependent seed
		uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
		rng.seed(ss);
		std::uniform_real_distribution<double> unif(0, 1);
		if(N_threads==0) // number of CPU cores not detected.
			N_threads=8;
	}

	void solve_init()
	{
		check_settings();
		shrink_scale=1.0;
		average_stall_count=0;
		best_stall_count=0;
		generation_step=-1;
		if(verbose)
		{
			std::cout<<"**************************************"<<std::endl;
			std::cout<<"*             GA started             *"<<std::endl;
			std::cout<<"**************************************"<<std::endl;
			std::cout<<"population: "<<population<<std::endl;
			std::cout<<"elite_count: "<<elite_count<<std::endl;
			std::cout<<"crossover_fraction: "<<crossover_fraction<<std::endl;
			std::cout<<"mutation_rate: "<<mutation_rate<<std::endl;
			std::cout<<"**************************************"<<std::endl;
		}
		Chronometer timer;
		timer.tic();

		thisGenerationType generation0;
		init_population(generation0);
		generation_step=0;
		finalize_objectives(generation0);
		rank_population(generation0); // used for ellite tranfre, crossover and mutation
		finalize_generation(generation0);
		if(!is_single_objective())
		{ // muti-objective
			update_ideal_objectives(generation0,true);
			extreme_objectives.clear();
			scalarized_objectives_min.clear();
		}
		generation0.exe_time=timer.toc();
		if(!user_request_stop)
		{
			generations_so_abs.push_back(thisGenSOAbs(generation0));
			report_generation(generation0);
		}
		last_generation=generation0;
	}

	StopReason solve_next_generation()
	{
		Chronometer timer;
		timer.tic();
		generation_step++;
		shrink_scale=set_shrink_scale(generation_step);
		thisGenerationType new_generation;
		transfer(new_generation);
		crossover_and_mutation(new_generation);

		finalize_objectives(new_generation);
		rank_population(new_generation);  // used for selection
		thisGenerationType selected_generation;
		select_population(new_generation,selected_generation);
		new_generation=selected_generation;
		rank_population(new_generation); // used for elite tranfre, crossover and mutation
		finalize_generation(new_generation);
		new_generation.exe_time=timer.toc();

		if(!user_request_stop)
		{
			generations_so_abs.push_back(thisGenSOAbs(new_generation));
			report_generation(new_generation);
		}
		last_generation=new_generation;

		return stop_critera();
	}

	StopReason solve()
	{
		StopReason stop=StopReason::Undefined;
		solve_init();
		while(stop==StopReason::Undefined)
			stop=solve_next_generation();
		show_stop_reason(stop);
		return stop;
	}

	std::string stop_reason_to_string(StopReason stop)
	{
		switch(stop)
		{
			case StopReason::Undefined:
				return "No-stop";
				break;
			case StopReason::MaxGenerations:
				return "Maximum generation reached";
				break;
			case StopReason::StallAverage:
				return "Average stalled";
				break;
			case StopReason::StallBest:
				return "Best stalled";
				break;
			case StopReason::UserRequest:
				return "User request";
				break;
			default:
				return "Unknown reason";
		}
	}

protected:

	void report_generation(const thisGenerationType &new_generation)
	{
		if(is_single_objective())
		{ // SO (including IGA)
			SO_report_generation(
				generation_step,
				new_generation,
				new_generation.chromosomes[new_generation.best_chromosome_index].genes
				);
		}
		else
		{
			MO_report_generation(
				generation_step,
				new_generation,
				new_generation.fronts[0]
				);
		}
	}

	void show_stop_reason(StopReason stop)
	{
		if(verbose)
		{
			std::cout<<"Stop criteria: ";
			if(stop==StopReason::Undefined)
				std::cout<<"There is a bug in this function";
			else
				std::cout<<stop_reason_to_string(stop);
			std::cout<<std::endl;
			std::cout<<"**************************************"<<std::endl;
		}
	}

	void transfer(thisGenerationType &new_generation)
	{
		if(user_request_stop)
			return ;

		if(!is_interactive())
		{ // add all members
			for(thisChromosomeType c:last_generation.chromosomes)
				new_generation.chromosomes.push_back(c);
		}
		else
		{
			// in IGA, the final evaluation is expensive
			// therefore, only elites would be transfered.
			for(int i=0;i<elite_count;i++)
				new_generation.chromosomes.push_back(last_generation.chromosomes[last_generation.sorted_indices[i]]);
		}
	}

	void finalize_generation(thisGenerationType &new_generation)
	{
		if(user_request_stop)
			return ;

		if(is_single_objective())
		{
			double best=new_generation.chromosomes[0].total_cost;
			double sum=0;
			new_generation.best_chromosome_index=0;

			for(uint i=0;i<new_generation.chromosomes.size();i++)
			{
				double current_cost=new_generation.chromosomes[i].total_cost;
				sum+=current_cost;
				if(current_cost<=best)
				{
					new_generation.best_chromosome_index=i;
					best=current_cost;
				}
				best=std::min(best,current_cost);
			}

			new_generation.best_total_cost=best;
			new_generation.average_cost=sum/double(new_generation.chromosomes.size());
		}
	}

	void check_settings()
	{
		if(is_interactive())
		{
			if(calculate_IGA_total_fitness==nullptr)
				throw std::runtime_error("calculate_IGA_total_fitness is null in interactive mode!");
			if(calculate_SO_total_fitness!=nullptr)
				throw std::runtime_error("calculate_SO_total_fitness is not null in interactive mode!");
			if(calculate_MO_objectives!=nullptr)
				throw std::runtime_error("calculate_MO_objectives is not null in interactive mode!");
			if(distribution_objective_reductions!=nullptr)
				throw std::runtime_error("distribution_objective_reductions is not null in interactive mode!");
			if(MO_report_generation!=nullptr)
				throw std::runtime_error("MO_report_generation is not null in interactive mode!");
			if(eval_genes_IGA==nullptr)
				throw std::runtime_error("eval_genes_IGA is null in interactive mode!");
			if(eval_genes!=nullptr)
				throw std::runtime_error("eval_genes is not null in interactive mode (use eval_genes_IGA instead)!");
		}
		else
		{
			if(calculate_IGA_total_fitness!=nullptr)
				throw std::runtime_error("calculate_IGA_total_fitness is not null in non-interactive mode!");
			if(eval_genes_IGA!=nullptr)
				throw std::runtime_error("eval_genes_IGA is not null in non-interactive mode!");
			if(eval_genes==nullptr)
				throw std::runtime_error("eval_genes is null!");
			if(is_single_objective())
			{
				if(calculate_SO_total_fitness==nullptr)
					throw std::runtime_error("calculate_SO_total_fitness is null in single objective mode!");
				if(calculate_MO_objectives!=nullptr)
					throw std::runtime_error("calculate_MO_objectives is not null in single objective mode!");
				if(distribution_objective_reductions!=nullptr)
					throw std::runtime_error("distribution_objective_reductions is not null in single objective mode!");
				if(MO_report_generation!=nullptr)
					throw std::runtime_error("MO_report_generation is not null in single objective mode!");
			}
			else
			{
				if(calculate_SO_total_fitness!=nullptr)
					throw std::runtime_error("calculate_SO_total_fitness is no null in multi-objective mode!");
				if(calculate_MO_objectives==nullptr)
					throw std::runtime_error("calculate_MO_objectives is null in multi-objective mode!");
				if(distribution_objective_reductions==nullptr)
					throw std::runtime_error("distribution_objective_reductions is null in multi-objective mode!");
				if(MO_report_generation==nullptr)
					throw std::runtime_error("MO_report_generation is null in multi-objective mode!");
			}
		}

		if(init_genes==nullptr)
			throw std::runtime_error("init_genes is not adjusted.");
		if(mutate==nullptr)
			throw std::runtime_error("mutate is not adjusted.");
		if(crossover==nullptr)
			throw std::runtime_error("crossover is not adjusted.");
		if(N_threads<1)
			throw std::runtime_error("Number of threads is below 1.");
		if(population<1)
			throw std::runtime_error("population is below 1.");
		if(is_single_objective())
		{ // SO (including IGA)
			if(SO_report_generation==nullptr)
				throw std::runtime_error("SO_report_generation is not adjusted while problem mode is single-objective");
			if(MO_report_generation!=nullptr)
				throw std::runtime_error("MO_report_generation is adjusted while problem mode is single-objective");
		}
		else
		{
			if(SO_report_generation!=nullptr)
				throw std::runtime_error("SO_report_generation is adjusted while problem mode is multi-objective");
			if(MO_report_generation==nullptr)
				throw std::runtime_error("MO_report_generation is not adjusted while problem mode is multi-objective");
		}
	}

	void select_population(const thisGenerationType &g,thisGenerationType &g2)
	{
		if(user_request_stop)
			return ;

		if(is_single_objective())
			select_population_SO(g,g2);
		else
			select_population_MO(g,g2);
	}

	void update_ideal_objectives(const thisGenerationType &g,bool reset)
	{
		if(user_request_stop)
			return ;

		if(is_single_objective())
			throw std::runtime_error("Wrong code A0812473247.");
		if(reset)
			ideal_objectives=distribution_objective_reductions(g.chromosomes[0].objectives);
		uint N_r_objectives=uint(ideal_objectives.n_rows);
		for(thisChromosomeType x:g.chromosomes)
		{
			arma::vec obj_reduced=distribution_objective_reductions(x.objectives);
			for(uint i=0;i<N_r_objectives;i++)
				if(obj_reduced(i)<ideal_objectives(i))
					ideal_objectives(i)=obj_reduced(i);
		}
	}

	void select_population_MO(const thisGenerationType &g,thisGenerationType &g2)
	{
		update_ideal_objectives(g,false);
		if(generation_step<=0)
		{
			g2=g;
			return ;
		}
		g2.chromosomes.clear();
		std::vector<arma::vec> zb_objectives;
		for(thisChromosomeType x:g.chromosomes)
			zb_objectives.push_back(distribution_objective_reductions(x.objectives)-ideal_objectives);
		scalarize_objectives(zb_objectives);
		arma::vec intercepts=build_hyperplane_intercepts();
		std::vector<arma::vec> norm_objectives;
		norm_objectives.reserve(g.chromosomes.size());
		for(uint i=0;i<g.chromosomes.size();i++)
		{
			norm_objectives.push_back(zb_objectives[i]/intercepts);
		}
		if(g.chromosomes.size()==population)
		{
			g2=g;
			return ;
		}
		if(reference_vectors.empty())
		{
			uint obj_dept=uint(distribution_objective_reductions(g.chromosomes[0].objectives).n_rows);
			reference_vectors=generate_referenceVectors(obj_dept,reference_vector_divisions);
		}
		std::vector<uint> associated_ref_vector;
		std::vector<double> distance_ref_vector;

		std::vector<uint> niche_count;
		arma::mat distances; // row: pop, col: ref_vec
		associate_to_references(
			g,
			norm_objectives,
			associated_ref_vector,
			distance_ref_vector,
			niche_count,
			distances);

		uint last_front_index=0;
		// select from best fronts as long as they are accommodated in the population
		while(g2.chromosomes.size()+g.fronts[last_front_index].size()<=population)
		{
			for(uint i:g.fronts[last_front_index])
				g2.chromosomes.push_back(g.chromosomes[i]);
			last_front_index++;
		}
		std::vector<uint> last_front=g.fronts[last_front_index];
		// select randomly from the next front
		std::vector<uint> to_add;
		while(g2.chromosomes.size()+to_add.size()<population)
		{
			if(!enable_reference_vectors)
			{ // disabling reference points
				uint msz=uint(last_front.size());
				uint to_add_index=uint(std::floor(msz*rand()));
				if(to_add_index>=msz)
					to_add_index=0;
				to_add.push_back(last_front[to_add_index]);
				last_front.erase(last_front.begin()+to_add_index);
				continue ;
			}

			uint min_niche_index=index_of_min(niche_count);
			std::vector<uint> min_vec_neighbors;
			for(uint i:last_front)
			{
				if(associated_ref_vector[i]==min_niche_index)
					min_vec_neighbors.push_back(i);
			}
			if(min_vec_neighbors.size()==0)
			{
				niche_count[min_niche_index]=uint(10*g.chromosomes.size()); // inf
				continue;
			}
			uint next_member_index=0;
			if(niche_count[min_niche_index]==0)
			{
				arma::vec nc=distances.col(min_niche_index);
				for(uint i:min_vec_neighbors)
					if(nc(i)<distances(next_member_index))
						next_member_index=i;
			}
			else
			{
				uint msz=uint(min_vec_neighbors.size());
				next_member_index=uint(std::floor(msz*rand()));
				if(next_member_index>=msz)
					next_member_index=0;
			}
			uint to_add_index=min_vec_neighbors[next_member_index];
			to_add.push_back(to_add_index);
			int to_del_front=-1;
			for(uint i=0;i<last_front.size();i++)
				if(last_front[i]==to_add_index)
					to_del_front=i;

			if(to_del_front>=0)
				last_front.erase(last_front.begin()+to_del_front);

			niche_count[min_niche_index]++;
		}
		for(uint i:to_add)
			g2.chromosomes.push_back(g.chromosomes[i]);
	}

	void associate_to_references(
		const thisGenerationType &gen,
		const std::vector<arma::vec> &norm_objectives,
		std::vector<uint> &associated_ref_vector,
		std::vector<double> &distance_ref_vector,
		std::vector<uint> &niche_count,
		arma::mat distances)
	{
		uint N_ref=uint(reference_vectors.size());
		uint N_x=uint(gen.chromosomes.size());
		niche_count.assign(N_ref, 0);
		distances.zeros(N_x,N_ref); // row: pop, col: ref_vec
		associated_ref_vector.assign(gen.chromosomes.size(),0);
		distance_ref_vector.assign(gen.chromosomes.size(),0.0);
		for(uint i=0;i<N_x;i++)
		{
			double dist_min=0.0;   // to avoid uninitialization warning
			uint dist_min_index=0; // to avoid uninitialization warning
			for(uint j=0;j<N_ref;j++)
			{
				arma::vec w=reference_vectors[j]/arma::norm(reference_vectors[j]);
				arma::vec norm_obj=norm_objectives[i];
				arma::mat wtnorm=w.t()*norm_obj;
				if(wtnorm.n_rows!=1 || wtnorm.n_cols!=1)
					throw std::runtime_error("unexpected matrix size A087624293!");
				double scalar_wtnorm=wtnorm(0,0);
				double dist=arma::norm(norm_obj-scalar_wtnorm*w);
				distances(i,j)=dist;
				if(j==0 || dist<dist_min)
				{
					dist_min=dist;
					dist_min_index=j;
				}
			}
			associated_ref_vector[i]=dist_min_index;
			distance_ref_vector[i]=dist_min;
			niche_count[dist_min_index]++;
		}
	}

	arma::vec build_hyperplane_intercepts()
	{
		uint N_objectives=uint(extreme_objectives.n_rows);
		arma::vec intercepts;
		arma::mat ones_vec=arma::ones<arma::vec>(N_objectives);
		intercepts=(1.0/arma::solve(extreme_objectives.t(),ones_vec));
		return intercepts;
	}

	template<typename T>
	uint index_of_min(const std::vector<T> &v)
	{
		return uint(std::distance(v.begin(), std::min_element(v.begin(), v.end())));
	}

	void scalarize_objectives(const std::vector<arma::vec> &zb_objectives)
	{
		uint N_objectives=uint(zb_objectives[0].n_rows);
		if(scalarized_objectives_min.is_empty())
		{
			extreme_objectives.zeros(N_objectives,N_objectives);
			scalarized_objectives_min.ones(N_objectives);
			scalarized_objectives_min*=arma::datum::inf;
		}
		for(uint i=0;i<N_objectives;i++)
		{
			arma::vec w;
			w.ones(N_objectives);
			w*=(1e-10);
			w(i)=1.0;
			std::vector<double> s;
			int N=int(zb_objectives.size());
			for(int j=0;j<N;j++)
				s.push_back((zb_objectives[j]/w).max());
			int min_sc_idx=index_of_min(s);
			double min_sc=s[min_sc_idx];

			if(min_sc<scalarized_objectives_min(i))
			{
				scalarized_objectives_min(i)=min_sc;
				extreme_objectives.col(i)=zb_objectives[min_sc_idx];
			}
		}

	}

	void select_population_SO(const thisGenerationType &g,thisGenerationType &g2)
	{
		if(generation_step<=0)
		{
			g2=g;
			return ;
		}

		if(verbose)
			std::cout<<"Transfered elites: ";
		std::vector<int> blocked;
		for(int i=0;i<elite_count;i++)
		{
			g2.chromosomes.push_back(g.chromosomes[g.sorted_indices[i]]);
			blocked.push_back(g.sorted_indices[i]);
			if(verbose)
			{
				std::cout<<(i==0?"":", ");
				std::cout<<(g.sorted_indices[i]+1);
			}
		}
		if(verbose)
			std::cout<<std::endl;
		for(int i=0;i<int(population)-elite_count;i++)
		{
			int j;
			bool allowed;
			do
			{
				allowed=true;
				j=select_parent(g);
				for(int k=0;k<int(blocked.size()) && allowed;k++)
					if(blocked[k]==j)
						allowed=false;
			} while(!allowed);
			g2.chromosomes.push_back(g.chromosomes[j]);
			blocked.push_back(g.sorted_indices[j]);
		}
		if(verbose)
			std::cout<<"Selection done."<<std::endl;
	}

	void rank_population(thisGenerationType &gen)
	{
		if(user_request_stop)
			return ;

		if(is_single_objective())
			rank_population_SO(gen);
		else
			rank_population_MO(gen);
	}

	void quicksort_indices_SO(std::vector<int> &array_indices,const thisGenerationType &gen,int left ,int right)
	{
		if(left<right)
		{
			int middle;
			double x=gen.chromosomes[array_indices[left]].total_cost;
			int l=left;
			int r=right;
			while(l<r)
			{
				while((gen.chromosomes[array_indices[l]].total_cost<=x)&&(l<right)) l++ ;
				while((gen.chromosomes[array_indices[r]].total_cost>x)&&(r>=left)) r-- ;
				if(l<r)
				{
					int temp = array_indices[l];
					array_indices[l]=array_indices[r];
					array_indices[r]=temp ;
				}
			}
			middle=r;
				int temp=array_indices[left];
				array_indices[left]=array_indices[middle];
				array_indices[middle]=temp;

			quicksort_indices_SO(array_indices,gen,left,middle-1);
			quicksort_indices_SO(array_indices,gen,middle+1,right);
		}
	}

	void rank_population_SO(thisGenerationType &gen)
	{
		int N=int(gen.chromosomes.size());
		gen.sorted_indices.clear();
		gen.sorted_indices.reserve(N);
		for(int i=0;i<N;i++)
			gen.sorted_indices.push_back(i);

		quicksort_indices_SO(gen.sorted_indices,gen,0,int(gen.sorted_indices.size())-1);

		std::vector<int> ranks;
		ranks.assign(gen.chromosomes.size(),0);
		for(uint i=0;i<gen.chromosomes.size();i++)
				ranks[gen.sorted_indices[i]]=i;

		generate_selection_chance(gen,ranks);
	}

	void generate_selection_chance(thisGenerationType &gen,const std::vector<int> &rank)
	{
		double chance_cumulative=0.0;
		uint N=uint(gen.chromosomes.size());
		gen.selection_chance_cumulative.clear();
		gen.selection_chance_cumulative.reserve(N);
		for(uint i=0;i<N;i++)
		{
			chance_cumulative+=1.0/sqrt(double(rank[i]+1));
			gen.selection_chance_cumulative.push_back(chance_cumulative);
		}
		for(uint i=0;i<N;i++)
		{	// normalizing
			gen.selection_chance_cumulative[i]=gen.selection_chance_cumulative[i]/gen.selection_chance_cumulative[population-1];
		}
	}

	void rank_population_MO(thisGenerationType &gen)
	{
		std::vector<std::vector<uint>> domination_set;
		std::vector<int> dominated_count;
		domination_set.reserve(gen.chromosomes.size());
		dominated_count.reserve(gen.chromosomes.size());
		for(uint i=0;i<gen.chromosomes.size();i++)
		{
			domination_set.push_back({});
			dominated_count.push_back(0);
		}
		std::vector<uint> pareto_front;

		for(uint i=0;i<gen.chromosomes.size();i++)
		{
			for(uint j=i+1;j<gen.chromosomes.size();j++)
			{
				if(Dominates(gen.chromosomes[i],gen.chromosomes[j]))
				{
					domination_set[i].push_back(j);
					dominated_count[j]++;
				}
				if(Dominates(gen.chromosomes[j],gen.chromosomes[i]))
				{
					domination_set[j].push_back(i);
					dominated_count[i]++;
				}
			}
			if(dominated_count[i]==0)
				pareto_front.push_back(i);
		}
		gen.fronts.clear();
		gen.fronts.push_back(pareto_front);
		std::vector<uint> next_front;
		do
		{
			next_front.clear();
			std::vector<uint> &last_front=gen.fronts[gen.fronts.size()-1];
			for(uint i:last_front)
				for(uint j:domination_set[i])
					if(--dominated_count[j]==0)
						next_front.push_back(j);
			if(!next_front.empty())
				gen.fronts.push_back(next_front);
		} while (!next_front.empty());
		std::vector<int> ranks;
		ranks.assign(gen.chromosomes.size(),0);
		for(uint i=0;i<gen.fronts.size();i++)
			for(uint j=0;j<gen.fronts[i].size();j++)
				ranks[gen.fronts[i][j]]=i;
		generate_selection_chance(gen,ranks);
	}

	bool Dominates(thisChromosomeType a,thisChromosomeType b)
	{
		if(a.objectives.n_rows!=b.objectives.n_rows)
			throw std::runtime_error("vector size mismatch A73592753!");
		for(uint i=0;i<a.objectives.n_rows;i++)
			if(a.objectives(i)>b.objectives(i))
				return false;
		for(uint i=0;i<a.objectives.n_rows;i++)
			if(a.objectives(i)<b.objectives(i))
				return true;
		return false;
	}

	std::vector<arma::vec> generate_integerReferenceVectors(int dept,int N_division)
	{
		std::vector<arma::vec> result;
		if(dept<1)
			throw std::runtime_error("wrong vector dept!");
		if(dept==1)
		{
			arma::vec v(1);
			v(0)=N_division;
			result.push_back(v);

			return result;
		}
		for(int i=0;i<=N_division;i++)
		{
			std::vector<arma::vec> tail;
			tail=generate_integerReferenceVectors(dept-1,N_division-i);

			for(int j=0;j<int(tail.size());j++)
			{
				arma::vec v1=tail[j];
				arma::vec v2(v1.n_rows+1);
				v2(0)=i;
				for(int k=0;k<int(v1.size());k++)
				{
					v2(k+1)=v1(k);
				}
				result.push_back(v2);
			}
		}
		return result;
	}

	std::vector<arma::vec> generate_referenceVectors(int dept,int N_division)
	{
		std::vector<arma::vec> A=generate_integerReferenceVectors(dept,N_division);
		for(int i=0;i<int(A.size());i++)
			A[i]=A[i]/arma::sum(A[i]);
		return A;
	}

	bool is_single_objective()
	{
		switch(problem_mode)
		{
			case GA_MODE::SOGA:		return true;
			case GA_MODE::IGA:		return true;
			case GA_MODE::NSGA_III:	return false;
			default:
				throw std::runtime_error("Code should not reach here!");
		}
	}

	bool is_interactive()
	{
		switch(problem_mode)
		{
			case GA_MODE::SOGA:		return false;
			case GA_MODE::IGA:		return true;
			case GA_MODE::NSGA_III:	return false;
			default:
				throw std::runtime_error("Code should not reach here!");
		}
	}

	void init_population_range(
		thisGenerationType *p_generation0,
		int index_begin,
		int index_end,
		uint *attemps,
		int *active_thread)
	{
		int dummy;
		for(int i=index_begin;i<=index_end;i++)
			init_population_single(p_generation0,i,attemps,&dummy);
		*active_thread=0; // false
	}

	void init_population_single(
		thisGenerationType *p_generation0,
		int index,
		uint *attemps,
		int *active_thread)
	{
		bool accepted=false;
		while(!accepted)
		{
			thisChromosomeType X;
			init_genes(X.genes,[this](){return rand();});
			if(is_interactive())
			{
				if(eval_genes_IGA(X.genes,X.middle_costs,*p_generation0))
				{
					// in IGA mode, code cannot run in parallel.
					p_generation0->chromosomes.push_back(X);
					accepted=true;
				}
			}
			else
			{
				if(eval_genes(X.genes,X.middle_costs))
				{
					if(index>=0)
						p_generation0->chromosomes[index]=X;
					else
						p_generation0->chromosomes.push_back(X);
					accepted=true;
				}
			}
			(*attemps)++;
		}
		*active_thread=0; //false
	}

	void idle()
	{
		if(custom_refresh!=nullptr)
			custom_refresh();
		if(idle_delay_us>0)
			boost::this_thread::sleep(boost::posix_time::microseconds(idle_delay_us));
	}

	void init_population(thisGenerationType &generation0)
	{
		generation0.chromosomes.clear();

		uint total_attempts=0;
		if(!multi_threading || N_threads==1 || is_interactive())
		{
			int dummy;
			for(uint i=0;i<population && !user_request_stop;i++)
				init_population_single(&generation0,-1,&total_attempts,&dummy);
		}
		else
		{
			for(uint i=0;i<population;i++)
				generation0.chromosomes.push_back(thisChromosomeType());
			std::vector<int> active_threads; // std::vector<bool> is broken
			active_threads.assign(N_threads,0);
			std::vector<uint> attempts;
			attempts.assign(N_threads,0);

			std::vector<std::thread> thread_pool;
			for(int i=0;i<N_threads;i++)
				thread_pool.push_back(std::thread());
			for(std::thread& th : thread_pool)
				if(th.joinable())
					th.join();

			if(dynamic_threading)
			{
				uint x_index=0;
				while(x_index<population && !user_request_stop)
				{
					int free_thread=-1;
					for(int i=0;i<N_threads && free_thread<0;i++)
					{
						if(!active_threads[i])
						{
							free_thread=i;
							if(thread_pool[free_thread].joinable())
								thread_pool[free_thread].join();
						}
					}
					if(free_thread>-1)
					{
						active_threads[free_thread]=1;
						thread_pool[free_thread]=
							std::thread(
								&std::remove_reference<decltype(*this)>::type::init_population_single,
								this,
								&generation0,
								int(x_index),
								&attempts[free_thread],
								&(active_threads[free_thread])
									);
						x_index++;
					}
					else
						idle();
				} // while
			} // endif: dynamic threading
			else
			{// on static threading
				int x_index_start=0;
				int x_index_end=0;
				int pop_chunk=population/N_threads;
				pop_chunk=std::max(pop_chunk,1);
				for(int i=0;i<N_threads;i++)
				{
					x_index_end=x_index_start+pop_chunk;
					if(i+1==N_threads) // last chunk
						x_index_end=population-1;
					else
						x_index_end=std::min(x_index_end,int(population)-1);

					if(x_index_end>=x_index_start)
					{
						active_threads[i]=1;
						thread_pool[i]=
						std::thread(
							&std::remove_reference<decltype(*this)>::type::init_population_range,
							this,
							&generation0,
							x_index_start,
							x_index_end,
							&attempts[i],
							&(active_threads[i])
								);
					}
					x_index_start=x_index_end+1;
				}
			} // endif: static threading

			bool all_tasks_finished;
			do
			{
				all_tasks_finished=true;
				for(int i=0;i<N_threads;i++)
					if(active_threads[i])
						all_tasks_finished=false;
				if(!all_tasks_finished)
					idle();
			}while(!all_tasks_finished);
			// wait for tasks to finish
			for(std::thread& th : thread_pool)
				if(th.joinable())
					th.join();

			for(uint ac:attempts)
				total_attempts+=ac;
		}

		/////////////////////

		if(verbose)
		{
			std::cout<<"Initial population of "<<population<<" was created with "<<total_attempts<<" attemps."<<std::endl;
		}
	}

	double rand()
	{
		return unif_dist(rng);
	}

	int select_parent(const thisGenerationType &g)
	{
		int N_max=int(g.chromosomes.size());
		double r=rand();
		int position=0;
		while(position<N_max && g.selection_chance_cumulative[position]<r)
			position++;
		return position;
	}

	void crossover_and_mutation_range(
		thisGenerationType *p_new_generation,
		uint pop_previous_size,
		int x_index_begin,
		int x_index_end,
		int *active_thread)
	{
		int dummy;
		for(int i=x_index_begin;i<=x_index_end;i++)
			crossover_and_mutation_single(p_new_generation,pop_previous_size,i,&dummy);
		*active_thread=0; // false
	}

	void crossover_and_mutation_single(
		thisGenerationType *p_new_generation,
		uint pop_previous_size,
		int index,
		int *active_thread)
	{

		if(verbose)
			std::cout<<"Action: crossover"<<std::endl;

		bool successful=false;
		while(!successful)
		{
			thisChromosomeType X;

			int pidx_c1=select_parent(last_generation);
			int pidx_c2=select_parent(last_generation);
			if(pidx_c1==pidx_c2)
				continue ;
			if(verbose)
				std::cout<<"Crossover of chromosomes "<<pidx_c1<<","<<pidx_c2<<std::endl;
			GeneType Xp1=last_generation.chromosomes[pidx_c1].genes;
			GeneType Xp2=last_generation.chromosomes[pidx_c2].genes;
			X.genes=crossover(Xp1,Xp2,[this](){return rand();});
			if(rand()<=mutation_rate)
			{
				if(verbose)
					std::cout<<"Mutation of chromosome "<<std::endl;
				X.genes=mutate(X.genes,[this](){return rand();},shrink_scale);
			}
			if(is_interactive())
			{
				if(eval_genes_IGA(X.genes,X.middle_costs,*p_new_generation))
				{
					p_new_generation->chromosomes.push_back(X);
					successful=true;
				}
			}
			else
			{
				if(eval_genes(X.genes,X.middle_costs))
				{
					if(index>=0)
						p_new_generation->chromosomes[pop_previous_size+index]=X;
					else
						p_new_generation->chromosomes.push_back(X);
					successful=true;
				}
			}
		}
		*active_thread=0; // false
	}

	void crossover_and_mutation(thisGenerationType &new_generation)
	{
		if(user_request_stop)
			return ;

		if(crossover_fraction<=0.0 || crossover_fraction>1.0)
			throw std::runtime_error("Wrong crossover fractoin");
		if(mutation_rate<0.0 || mutation_rate>1.0)
			throw std::runtime_error("Wrong mutation rate");
		if(generation_step<=0)
			return ;
		uint N_add=uint(std::round(double(population)*(crossover_fraction)));
		uint pop_previous_size=uint(new_generation.chromosomes.size());
		if(is_interactive())
		{
			if(N_add+elite_count!=population)
				throw std::runtime_error("In IGA mode, elite fraction + crossover fraction + mutation fraction must be equal to 1.0 !");
		}

		if(!multi_threading || N_threads==1 || is_interactive())
		{
			int dummy;
			for(uint i=0;i<N_add && !user_request_stop;i++)
				crossover_and_mutation_single(&new_generation,pop_previous_size,-1,&dummy);
		}
		else
		{
			for(uint i=0;i<N_add;i++)
				new_generation.chromosomes.push_back(thisChromosomeType());
			std::vector<int> active_threads; // std::vector<bool> is broken
			active_threads.assign(N_threads,0);

			std::vector<std::thread> thread_pool;
			for(int i=0;i<N_threads;i++)
				thread_pool.push_back(std::thread());
			for(std::thread& th : thread_pool)
				if(th.joinable())
					th.join();

			if(dynamic_threading)
			{
				uint x_index=0;
				while(x_index<N_add && !user_request_stop)
				{
					int free_thread=-1;
					for(int i=0;i<N_threads && free_thread<0;i++)
					{
						if(!active_threads[i])
						{
							free_thread=i;
							if(thread_pool[free_thread].joinable())
								thread_pool[free_thread].join();
						}
					}
					if(free_thread>-1)
					{
						active_threads[free_thread]=1;
						thread_pool[free_thread]=
							std::thread(
								&std::remove_reference<decltype(*this)>::type::crossover_and_mutation_single,
								this,
								&new_generation,
								pop_previous_size,
								int(x_index),
								&(active_threads[free_thread])
									);
						x_index++;
					}
					else
						idle();
				}
			}// endif: dynamic threading
			else
			{// on static threading
				int x_index_start=0;
				int x_index_end=0;
				int pop_chunk=N_add/N_threads;
				pop_chunk=std::max(pop_chunk,1);
				for(int i=0;i<N_threads;i++)
				{
					x_index_end=x_index_start+pop_chunk;
					if(i+1==N_threads) // last chunk
						x_index_end=N_add-1;
					else
						x_index_end=std::min(x_index_end,int(N_add)-1);

					if(x_index_end>=x_index_start)
					{
						active_threads[i]=1;

						thread_pool[i]=
							std::thread(
								&std::remove_reference<decltype(*this)>::type::crossover_and_mutation_range,
								this,
								&new_generation,
								pop_previous_size,
								x_index_start,
								x_index_end,
								&(active_threads[i])
									);
					}
					x_index_start=x_index_end+1;
				}
			}// endif: static threading

			bool all_tasks_finished;
			do
			{
				all_tasks_finished=true;
				for(int i=0;i<N_threads;i++)
					if(active_threads[i])
						all_tasks_finished=false;
				if(!all_tasks_finished)
					idle();
			}while(!all_tasks_finished);

			// wait for tasks to finish
			for(std::thread& th : thread_pool)
				if(th.joinable())
					th.join();
		}
	}

	StopReason stop_critera()
	{
		if(generation_step<2 && !user_request_stop)
			return StopReason::Undefined;



		if(is_single_objective())
		{
			const thisGenSOAbs &g1=generations_so_abs[int(generations_so_abs.size())-2];
			const thisGenSOAbs &g2=generations_so_abs[int(generations_so_abs.size())-1];

			if(std::abs(g1.best_total_cost-g2.best_total_cost)<tol_stall_best)
				best_stall_count++;
			else
				best_stall_count=0;
			if(std::abs(g1.average_cost-g2.average_cost)<tol_stall_average)
				average_stall_count++;
			else
				average_stall_count=0;
		}

		if(generation_step>=generation_max)
			return StopReason::MaxGenerations;

		if(average_stall_count>=average_stall_max)
			return StopReason::StallAverage;

		if(best_stall_count>=best_stall_max)
			return StopReason::StallBest;

		if(user_request_stop)
			return StopReason::UserRequest;

		return StopReason::Undefined;
	}

	void finalize_objectives(thisGenerationType &g)
	{
		if(user_request_stop)
			return ;

		switch(problem_mode)
		{
			case GA_MODE::SOGA:
				for(int i=0;i<int(g.chromosomes.size());i++)
					g.chromosomes[i].total_cost=calculate_SO_total_fitness(g.chromosomes[i]);
				break;
			case GA_MODE::IGA:
				calculate_IGA_total_fitness(g);
				break;
			case GA_MODE::NSGA_III:
				for(uint i=0;i<g.chromosomes.size();i++)
					g.chromosomes[i].objectives=calculate_MO_objectives(g.chromosomes[i]);
				break;
			default:
				throw std::runtime_error("Code should not reach here!");
		}

	}

};

NS_EA_END
