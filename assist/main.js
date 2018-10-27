
$(function() {
    add_more_var();
    add_more_var();
    add_more_obj();
    update_codes();
});

var var_count=0;
var obj_count=0;

function add_more_var()
{
	var content='';
	content+=var_options();
	var_count++;
	content+='<input type="text" onchange="update_codes()" class="var_name" placeholder="e.g. x" value="var'+var_count.toString()+'">';
	content+='<input type="text" onchange="update_codes()" class="var_min" placeholder="min" value="0.0">';
	content+='<input type="text" onchange="update_codes()" class="var_max" placeholder="max" value="10.0">';
	content+='<button onclick="delete_var(this)">x</button>';
	content='<div class="var">'+content+'</div>';
	$('#opt_var').append(content);
	update_codes();
}

function add_more_obj()
{
	var content='';
	obj_count++;
	content+='<input class="obj_cell_name" onchange="update_codes()" type="text" placeholder="e.g. my_objective" value="objective'+obj_count.toString()+'">';
	content+='= <input class="obj_cell_val" onchange="update_codes()" type="text" placeholder="e.g. x*y+tan(x)+sin(y)" value="'+random_eval()+'"><button onclick="delete_obj(this)">x</button>';
	content='<div class="obj">'+content+'</div>';
	$('#opt_obj').append(content);
	update_codes();
}

function delete_var(button)
{
	$(button).parent().remove();
	update_codes();
}

function delete_obj(button)
{
	$(button).parent().remove();
	update_codes();
}

function var_options()
{
	var content='';
	content+='<select onchange="update_codes()">';
	content+='	<option value="double" selected>double</option>';
	content+='	<option value="int">integer</option>';
	content+='	<option value="bool">boolean</option>';
	content+='	<option value="string">string</option>';
	content+='</select>';
	return content;
}

function rand_range(min, max)
{
    return Math.floor(Math.random() * (max - min + 1)) + min;
}

function random_eval()
{
	var variables=[];
	$('#opt_var').children().find('input.var_name').each(function(index,elem){
		var name=$(elem).val();
		if(name)
			variables.push($(elem).val());
	});
	if(variables.length<2)
	{
		variables.push('x');
		variables.push('y');
	}
	var functions=['sin','cos','sqrt','log','exp'];
	var phrases='';
	var phrase_count=rand_range(2,3);
	for(var i=0;i<phrase_count;i++)
	{
		if(i>0)
			phrases+='+';
		var idx1=rand_range(0,variables.length-1);
		var idx2=idx1;
		while(idx2==idx1) /* make sure operands are different*/
			idx2=rand_range(0,variables.length-1);

		var opr1=variables[idx1];
		var opr2=variables[idx2];
		if(Math.random()>0.1)
			opr1=functions[rand_range(0,functions.length-1)]+'('+opr1+')';
		if(Math.random()>0.1)
			opr2=functions[rand_range(0,functions.length-1)]+'('+opr2+')';
		phrases+=opr1+'*'+opr2;
	}
	return phrases;
}

function escapeHTML(unsafe_text)
{
    let div = document.createElement('div');
    div.innerText = unsafe_text;
    return div.innerHTML;
}

function update_codes()
{
	var code=[];
	code.push('// main.cpp');
	code.push('');
	code.push('#include <string>');
	code.push('#include <iostream>');
	code.push('#include <fstream>');
	code.push('#include "openga.hpp"');
	code.push('');
	code.push('using std::string;');
	code.push('using std::cout;');
	code.push('using std::endl;');
	code.push('');
	var multiobj=false;
	switch($('select[name=objnum]').val())
	{
		case "single":
			multiobj=false;
			break;
		case "multiple":
			multiobj=true;
			break;
		default:
			multiobj=false;
	}
	var solution_name=$('input[name=solution_name]').val();
	solution_name=solution_name||'MySolution';
	var eval_name=$('input[name=eval_name]').val();
	eval_name=eval_name||'MyMiddleCost';
	// solution
	code.push('struct '+solution_name);
	code.push('{');
	var var_list=[];
	$('#opt_var').children().each(function(index,elem){
		var var_name=$(elem).find('input.var_name').val();
		var var_type=$(elem).find('select').val();
		var_name=var_name||'(unknown variable name)';
		var_type=var_type||'(unknown variable type)';
		var_list.push(var_name);
		code.push('\t'+var_type+' '+var_name+';');
	});
	code.push('');
	code.push('	string to_string() const');
	code.push('	{');
	code.push('		return ');
	code.push('			string("{")');
	var first_var=true;
	var_list.forEach(function(name){
		var prefix=first_var?'+  "':'+", ';
		code.push('\t\t\t'+prefix+name+':"+std::to_string('+name+')');
		first_var=false;
	});
	code.push('			+"}";');
	code.push('	}');
	code.push('};');
	code.push('');
	// middle costs
	code.push('struct '+eval_name);
	code.push('{');
	code.push('	// This is where the results of simulation');
	code.push('	// is stored but not yet finalized.');
	var obj_list=[];
	$('#opt_obj').children().each(function(index,elem){
		var obj_name=$(elem).find('input.obj_cell_name').val();
		obj_name=obj_name||'(unknown variable name)';
		obj_list.push(obj_name);
		code.push('\t'+'double '+obj_name+';');
	});
	code.push('};');
	code.push('');
	code.push('typedef EA::Genetic<'+solution_name+','+eval_name+'> GA_Type;');
	code.push('typedef EA::GenerationType<'+solution_name+','+eval_name+'> Generation_Type;');
	code.push('');
	// init_genes
	code.push('void init_genes('+solution_name+'& p,const std::function<double(void)> &rnd01)');
	code.push('{');
	code.push('\t// rnd01() gives a random number in 0~1');
	$('#opt_var').children().each(function(index,elem){
		var var_name=$(elem).find('input.var_name').val();
		var var_min=$(elem).find('input.var_min').val();
		var var_max=$(elem).find('input.var_max').val();
		var_name=var_name||'(unknown variable name)';
		var_min=var_min||'0.0';
		var_max=var_max||'10.0';
		var diff=Number(var_max)-Number(var_min);
		var diff_string=diff.toString();
		// if(Math.floor(diff)==diff)
		// 	diff_string+='.0';
		if(diff<0)
			code.push('\t'+'warning: for variable '+var_name+' max<min');
		if(diff==0)
			code.push('\t'+'warning: for variable '+var_name+' max==min');
		code.push('\t'+'p.'+var_name+'='+var_min+'+'+diff_string+'*rnd01();');
	});
	code.push('}');
	code.push('');
	// eval_solution
	code.push('bool eval_solution(');
	code.push('	const '+solution_name+'& p,');
	code.push('	'+eval_name+' &c)');
	code.push('{');
	$('#opt_var').children().each(function(index,elem){
		var var_name=$(elem).find('input.var_name').val();
		var var_type=$(elem).find('select').val();
		var_name=var_name||'(unknown variable name)';
		var_type=var_type||'(unknown variable type)';
		var_list.push(var_name);
		code.push('\t'+'const '+var_type+'& '+var_name+'=p.'+var_name+';');
	});
	code.push('');
	$('#opt_obj').children().each(function(index,elem){
		var obj_name=$(elem).find('input.obj_cell_name').val();
		var obj_val=$(elem).find('input.obj_cell_val').val();
		obj_name=obj_name||'(unknown variable name)';
		obj_val=obj_val||'(unknown variable value)';
		code.push('\t'+'c.'+obj_name+'='+obj_val+';');
	});
	code.push('	return true; // solution is accepted');
	code.push('}');	
	code.push('');
	// mutate
	code.push(solution_name+' mutate(');
	code.push('	const '+solution_name+'& X_base,');
	code.push('	const std::function<double(void)> &rnd01,');
	code.push('	double shrink_scale)');
	code.push('{');
	code.push('	'+solution_name+' X_new;');
	code.push('	bool in_range;');
	code.push('	do{');
	code.push('		in_range=true;');
	code.push('		X_new=X_base;');
	$('#opt_var').children().each(function(index,elem){
		var var_name=$(elem).find('input.var_name').val();
		var var_min=$(elem).find('input.var_min').val();
		var var_max=$(elem).find('input.var_max').val();
		var_name=var_name||'(unknown variable name)';
		var_min=var_min||'0.0';
		var_max=var_max||'10.0';
		// if(Math.floor(Number(var_min))==Number(var_min))
		// 	var_min+='.0';
		// if(Math.floor(Number(var_max))==Number(var_max))
		// 	var_max+='.0';
		code.push('		X_new.'+var_name+'+=0.2*(rnd01()-rnd01())*shrink_scale;');
		code.push('		in_range=in_range&&(X_new.'+var_name+'>='+var_min+' && X_new.'+var_name+'<'+var_max+');');
	});
	code.push('	} while(!in_range);');
	code.push('	return X_new;');
	code.push('}');
	code.push('');
	// crossover
	code.push(solution_name+' crossover(');
	code.push('	const '+solution_name+'& X1,');
	code.push('	const '+solution_name+'& X2,');
	code.push('	const std::function<double(void)> &rnd01)');
	code.push('{');
	code.push('	'+solution_name+' X_new;');
	code.push('	double r;');
	$('#opt_var').children().each(function(index,elem){
		var var_name=$(elem).find('input.var_name').val();
		var_name=var_name||'(unknown variable name)';
		code.push('	r=rnd01();');
		code.push('\t'+'X_new.'+var_name+'=r*X1.'+var_name+'+(1.0-r)*X2.'+var_name+';');
	});
	code.push('	return X_new;');
	code.push('}');
	code.push('');
	// calculate_SO_total_fitness
	if(multiobj)
	{
		code.push('std::vector<double> calculate_MO_objectives(const GA_Type::thisChromosomeType &X)');
		code.push('{');
		code.push('	return {');
		var array_objs=[];
		$('#opt_obj').children().each(function(index,elem){
			var obj_name=$(elem).find('input.obj_cell_name').val();
			obj_name=obj_name||'(unknown variable name)';
			array_objs.push('\t\t'+'X.middle_costs.'+obj_name);
		});
		code.push(array_objs.join(',\n'));
		code.push('	};');
		code.push('}');
	}
	else
	{
		code.push('double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)');
		code.push('{');
		code.push('	// finalize the cost');
		code.push('	double final_cost=0.0;');
		$('#opt_obj').children().each(function(index,elem){
			var obj_name=$(elem).find('input.obj_cell_name').val();
			var obj_val=$(elem).find('input.obj_cell_val').val();
			obj_name=obj_name||'(unknown variable name)';
			obj_val=obj_val||'(unknown variable value)';
			code.push('\t'+'final_cost+=X.middle_costs.'+obj_name+';');
		});
		code.push('	return final_cost;');
		code.push('}');
	}

	code.push('');
	code.push('std::ofstream output_file;');
	code.push('');
	// Report Generation
	if(multiobj)
	{
		code.push('void MO_report_generation(');
		code.push('	int generation_number,');
		code.push('	const EA::GenerationType<'+solution_name+','+eval_name+'> &last_generation,');
		code.push('	const std::vector<unsigned int>& pareto_front)');
		code.push('{');
		code.push('	(void) last_generation;');
		code.push('');
		code.push('	cout<<"Generation ["<<generation_number<<"], ";');
		code.push('	cout<<"Pareto-Front {";');
		code.push('	for(unsigned int i=0;i<pareto_front.size();i++)');
		code.push('	{');
		code.push('		cout<<(i>0?",":"");');
		code.push('		cout<<pareto_front[i];');
		code.push('	}');
		code.push('	cout<<"}"<<endl;');
		code.push('}');
		code.push('');
		code.push('void save_results(const GA_Type &ga_obj)');
		code.push('{');
		code.push('	std::ofstream output_file;');
		code.push('	output_file.open("paretofront.txt");');
		code.push('	output_file');
		code.push('			<<"N"');
		$('#opt_obj').children().each(function(index,elem){
			var obj_name=$(elem).find('input.obj_cell_name').val();
			obj_name=obj_name||'(unknown variable name)';
			code.push('\t\t\t'+'<<"\\t"<<"'+obj_name+'"');
		});
		code.push('			<<"\\t"<<"solution"<<"\\n";');
		code.push('	std::vector<unsigned int> paretofront_indices=ga_obj.last_generation.fronts[0];');
		code.push('	for(unsigned int i:paretofront_indices)');
		code.push('	{');
		code.push('		const auto &X=ga_obj.last_generation.chromosomes[i];');
		code.push('		output_file');
		code.push('			<<i<<"\\t"');
		$('#opt_obj').children().each(function(index,elem){
			var obj_name=$(elem).find('input.obj_cell_name').val();
			obj_name=obj_name||'(unknown variable name)';
			code.push('\t\t\t'+'<<X.middle_costs.'+obj_name+'<<"\\t"');
		});
		code.push('			<<X.genes.to_string()<<"\\n";');
		code.push('');
		code.push('	}');
		code.push('	output_file.close();');
		code.push('}		');
	}
	else
	{
		code.push('void SO_report_generation(');
		code.push('	int generation_number,');
		code.push('	const EA::GenerationType<'+solution_name+','+eval_name+'> &last_generation,');
		code.push('	const '+solution_name+'& best_genes)');
		code.push('{');
		code.push('	cout');
		code.push('		<<"Generation ["<<generation_number<<"], "');
		code.push('		<<"Best="<<last_generation.best_total_cost<<", "');
		code.push('		<<"Average="<<last_generation.average_cost<<", "');
		code.push('		<<"Best genes=("<<best_genes.to_string()<<")"<<", "');
		code.push('		<<"Exe_time="<<last_generation.exe_time');
		code.push('		<<endl;');
		code.push('');
		code.push('	output_file');
		code.push('		<<generation_number<<"\\t"');
		code.push('		<<last_generation.average_cost<<"\\t"');
		code.push('		<<last_generation.best_total_cost<<"\\t"');
		code.push('		<<best_genes.to_string()<<"\\n";');
		code.push('}');
	}
	code.push('');
	// main.cpp
	code.push('int main()');
	code.push('{');
	if(!multiobj)
	{
		code.push('	output_file.open("results.txt");');
		code.push('	output_file<<"step"<<"\\t"<<"cost_avg"<<"\\t"<<"cost_best"<<"\\t"<<"solution_best"<<"\\n";');
		code.push('');
	}
	code.push('	EA::Chronometer timer;');
	code.push('	timer.tic();');
	code.push('');
	code.push('	GA_Type ga_obj;');
	if(multiobj)
		code.push('	ga_obj.problem_mode=EA::GA_MODE::NSGA_III;');
	else
		code.push('	ga_obj.problem_mode=EA::GA_MODE::SOGA;');
	var heavyness=-1;
	switch($('select[name=heavyness]').val())
	{
		case "light":
			heavyness=0;
			break;
		case "medium":
			heavyness=1;
			break;
		case "heavy":
			heavyness=2;
			break;
		default:
			heavyness=-1;
	}
	var multi_threading=(heavyness>=1?'true':'false');
	code.push('	ga_obj.multi_threading='+multi_threading+';');
	if(heavyness==1)
		code.push('	ga_obj.idle_delay_us=1; // switch between threads quickly');
	if(heavyness==2)
		code.push('	ga_obj.idle_delay_us=10; // switch between threads quickly');
	if(heavyness>=1)
	{
		var dynamic_threading=(heavyness>1?'true':'false');
		code.push('	ga_obj.dynamic_threading='+dynamic_threading+';');
	}
	var verbose=($('select[name=verbose]').val()==="true" ? 'true' : 'false');
	code.push('	ga_obj.verbose='+verbose+';');
	var population=-1;
	switch($('select[name=population]').val())
	{

		case "small":
			population=50;
			break;
		case "medium":
			population=200;
			break;
		case "large":
			population=1000;
			break;
		default:
			population=-1;
	}
	if(multiobj && population<100)
		population=100;
	code.push('	ga_obj.population='+population.toString()+';');
	code.push('	ga_obj.generation_max=1000;');
	if(multiobj)
		code.push('	ga_obj.calculate_MO_objectives=calculate_MO_objectives;');
	else
		code.push('	ga_obj.calculate_SO_total_fitness=calculate_SO_total_fitness;');
	code.push('	ga_obj.init_genes=init_genes;');
	code.push('	ga_obj.eval_solution=eval_solution;');
	code.push('	ga_obj.mutate=mutate;');
	code.push('	ga_obj.crossover=crossover;');
	if(multiobj)
		code.push('	ga_obj.MO_report_generation=MO_report_generation;');
	else
		code.push('	ga_obj.SO_report_generation=SO_report_generation;');
	code.push('	ga_obj.best_stall_max=10;');
	code.push('	ga_obj.elite_count=10;');
	code.push('	ga_obj.crossover_fraction=0.7;');
	code.push('	ga_obj.mutation_rate=0.2;');
	if(!multiobj)
	{
		code.push('	ga_obj.best_stall_max=10;');
		code.push('	ga_obj.elite_count=10;');
	}
	code.push('	ga_obj.solve();');
	code.push('');
	code.push('	cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<endl;');
	code.push('');
	if(!multiobj)
		code.push('	output_file.close();');
	else
		code.push('	save_results(ga_obj);');
	code.push('	return 0;');
	code.push('}');
	code.push('');

	$('#code').html('<pre>'+escapeHTML(code.join('\n'))+'</pre>');
}