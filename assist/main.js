document.addEventListener('DOMContentLoaded', function() { // Wait for DOM to be ready
    add_more_var();
    add_more_var();
    add_more_obj();
    update_codes();

    // Listen to all inputs
    const inputs = document.querySelectorAll('select, input');
    inputs.forEach(input => {
        input.addEventListener('change', update_codes);
        if(input.tagName === 'INPUT') {
            input.addEventListener('keyup', update_codes);
        }
    });
});

let var_count = 0;
let obj_count = 0;

function add_more_var()
{
    let content = '';
    content += var_options();
    var_count++;
    content += '<input type="text" onchange="update_codes()" onkeyup="update_codes()" class="var_name" placeholder="e.g. x" value="var' + var_count.toString() + '">';
    content += '<input type="text" onchange="update_codes()" onkeyup="update_codes()" class="var_min" placeholder="min" value="0.0">';
    content += '<input type="text" onchange="update_codes()" onkeyup="update_codes()" class="var_max" placeholder="max" value="10.0">';
    content += '<button onclick="delete_var(this)" class="btn-delete">✖</button>';
    content = '<div class="var">' + content + '</div>';
    document.getElementById('opt_var').appendChild(htmlToElement(content));
    update_codes();
}

function add_more_obj()
{
    let content = '';
    obj_count++;
    content += '<input class="obj_cell_name" onchange="update_codes()" onkeyup="update_codes()" type="text" placeholder="e.g. my_objective" value="objective' + obj_count.toString() + '">';
    content += ' = <input class="obj_cell_val" onchange="update_codes()" onkeyup="update_codes()" type="text" placeholder="e.g. x*y+tan(x)+sin(y)" value="' + random_eval() + '">';
    content += '<button onclick="delete_obj(this)" class="btn-delete">✖</button>';
    content = '<div class="obj">' + content + '</div>';
    document.getElementById('opt_obj').appendChild(htmlToElement(content));
    update_codes();
}

function delete_var(button)
{
    button.parentElement.remove();
    update_codes();
}

function delete_obj(button)
{
    button.parentElement.remove();
    update_codes();
}

function var_options()
{
    let content = '';
    content += '<select onchange="update_codes()" class="var-type-select">';
    content += '  <option value="double" selected>double</option>';
    content += '  <option value="int">integer</option>';
    content += '  <option value="bool">boolean</option>';
    content += '  <option value="string">string</option>';
    content += '</select>';
    return content;
}

function rand_range(min, max)
{
    return Math.floor(Math.random() * (max - min + 1)) + min;
}

function random_eval()
{
    const variables = [];
    const varElements = document.querySelectorAll('#opt_var .var');
    varElements.forEach(function(elem) {
        const var_name = elem.querySelector('.var_name').value;
        if(var_name)
            variables.push(var_name);
    });
    if(variables.length < 2)
    {
        variables.push('x');
        variables.push('y');
    }
    const functions = ['sin', 'cos', 'sqrt', 'log', 'exp'];
    let phrases = '';
    const phrase_count = rand_range(2, 3);
    for(let i = 0; i < phrase_count; i++)
    {
        if(i > 0)
            phrases += '+';
        let idx1 = rand_range(0, variables.length - 1);
        let idx2 = idx1;
        while(idx2 == idx1)
            idx2 = rand_range(0, variables.length - 1);

        let opr1 = variables[idx1];
        let opr2 = variables[idx2];
        if(Math.random() > 0.1)
            opr1 = functions[rand_range(0, functions.length - 1)] + '(' + opr1 + ')';
        if(Math.random() > 0.1)
            opr2 = functions[rand_range(0, functions.length - 1)] + '(' + opr2 + ')';
        phrases += opr1 + '*' + opr2;
    }
    return phrases;
}

function htmlToElement(html)
{
    const template = document.createElement('template');
    html = html.trim();
    template.innerHTML = html;
    return template.content.firstChild;
}

function escapeHTML(unsafe_text)
{
    const div = document.createElement('div');
    div.innerText = unsafe_text;
    return div.innerHTML;
}

function copyCode()
{
    const codeElement = document.getElementById('code');
    const codeText = codeElement.innerText || codeElement.textContent;

    navigator.clipboard.writeText(codeText).then(function() {
        const status = document.getElementById('copy-status');
        status.textContent = '✓ Copied!';
        status.classList.add('show');
        setTimeout(function() {
            status.classList.remove('show');
        }, 2000);
    }).catch(function() {
        const status = document.getElementById('copy-status');
        status.textContent = '✗ Failed to copy';
        status.classList.add('show');
        setTimeout(function() {
            status.classList.remove('show');
        }, 2000);
    });
}

function getSelectValue(id)
{
    const select = document.getElementById(id);
    return select ? select.value : null;
}

function getInputValue(id)
{
    const input = document.getElementById(id);
    return input ? input.value : null;
}

function update_codes()
{
    const code = [];
    code.push('// main.cpp');
    code.push('');
    code.push('#include <string>');
    code.push('#include <iostream>');
    code.push('#include <fstream>');
    code.push('#include "openGA.hpp"');
    code.push('');
    code.push('using std::string;');
    code.push('using std::cout;');
    code.push('using std::endl;');
    code.push('');
    const multiobj = (getSelectValue('objnum') === 'multiple');
    const solution_name = getInputValue('solution_name') || 'MySolution';
    const eval_name = getInputValue('eval_name') || 'MyMiddleCost';
    // solution
    code.push('struct ' + solution_name);
    code.push('{');
    const var_list = [];
    const varElements = document.querySelectorAll('#opt_var .var');
    varElements.forEach(function(elem) {
        const var_name = elem.querySelector('.var_name').value || '(unknown variable name)';
        const var_type = elem.querySelector('.var-type-select').value || '(unknown variable type)';
        var_list.push(var_name);
        code.push('    ' + var_type + ' ' + var_name + ';');
    });
    code.push('');
    code.push('    std::string to_string() const');
    code.push('    {');
    code.push('        return');
    let first_var = true;
    var_list.forEach(function(var_name) {
        const prefix = first_var ? '            std::string("{")' : '            + ", ';
        code.push(prefix + var_name + ':" + std::to_string(' + var_name + ')');
        first_var = false;
    });
    code.push('            + "}";');
    code.push('    }');
    code.push('};');
    code.push('');
    // middle costs
    code.push('struct ' + eval_name);
    code.push('{');
    code.push('    // This is where the results of simulation');
    code.push('    // is stored but not yet finalized.');
    const obj_list = [];
    const objElements = document.querySelectorAll('#opt_obj .obj');
    objElements.forEach(function(elem) {
        const obj_name = elem.querySelector('.obj_cell_name').value || '(unknown variable name)';
        obj_list.push(obj_name);
        code.push('    double ' + obj_name + ';');
    });
    code.push('};');
    code.push('');
    code.push('typedef EA::Genetic<' + solution_name + ', ' + eval_name + '> GA_Type;');
    code.push('typedef EA::GenerationType<' + solution_name + ', ' + eval_name + '> Generation_Type;');
    code.push('');
    // init_genes
    code.push('void init_genes(' + solution_name + '& p, const std::function<double(void)> &rnd01)');
    code.push('{');
    code.push('    // rnd01() gives a random number in 0~1');
    varElements.forEach(function(elem) {
        const var_name = elem.querySelector('.var_name').value || '(unknown variable name)';
        const var_min = elem.querySelector('.var_min').value || '0.0';
        const var_max = elem.querySelector('.var_max').value || '10.0';
        const diff = Number(var_max) - Number(var_min);
        const diff_string = diff.toString();
        if(diff < 0)
            code.push('    // warning: for variable ' + var_name + ' max<min');
        if(diff == 0)
            code.push('    // warning: for variable ' + var_name + ' max==min');
        code.push('    p.' + var_name + ' = ' + var_min + ' + ' + diff_string + ' * rnd01();');
    });
    code.push('}');
    code.push('');
    // eval_solution
    code.push('bool eval_solution(');
    code.push('    const ' + solution_name + '& p,');
    code.push('    ' + eval_name + ' &c)');
    code.push('{');
    varElements.forEach(function(elem) {
        const var_name = elem.querySelector('.var_name').value || '(unknown variable name)';
        const var_type = elem.querySelector('.var-type-select').value || '(unknown variable type)';
        var_list.push(var_name);
        code.push('    const ' + var_type + '& ' + var_name + ' = p.' + var_name + ';');
    });
    code.push('');
    objElements.forEach(function(elem) {
        const obj_name = elem.querySelector('.obj_cell_name').value || '(unknown variable name)';
        const obj_val = elem.querySelector('.obj_cell_val').value || '(unknown variable value)';
        code.push('    c.' + obj_name + ' = ' + obj_val + ';');
    });
    code.push('    return true; // solution is accepted');
    code.push('}');
    code.push('');
    // mutate
    code.push(solution_name + ' mutate(');
    code.push('    const ' + solution_name + '& X_base,');
    code.push('    const std::function<double(void)> &rnd01,');
    code.push('    double shrink_scale)');
    code.push('{');
    code.push('    ' + solution_name + ' X_new;');
    code.push('    const double mu = 0.2 * shrink_scale; // mutation radius (adjustable)');
    code.push('    bool in_range;');
    code.push('    do {');
    code.push('        in_range = true;');
    code.push('        X_new = X_base;');
    varElements.forEach(function(elem) {
        const var_name = elem.querySelector('.var_name').value || '(unknown variable name)';
        const var_min = elem.querySelector('.var_min').value || '0.0';
        const var_max = elem.querySelector('.var_max').value || '10.0';
        code.push('        X_new.' + var_name + ' += mu * (rnd01() - rnd01());');
        code.push('        in_range = in_range && (X_new.' + var_name + ' >= ' + var_min + ' && X_new.' + var_name + ' < ' + var_max + ');');
    });
    code.push('    } while(!in_range);');
    code.push('    return X_new;');
    code.push('}');
    code.push('');
    // crossover
    code.push(solution_name + ' crossover(');
    code.push('    const ' + solution_name + '& X1,');
    code.push('    const ' + solution_name + '& X2,');
    code.push('    const std::function<double(void)> &rnd01)');
    code.push('{');
    code.push('    ' + solution_name + ' X_new;');
    code.push('    double r;');
    varElements.forEach(function(elem) {
        const var_name = elem.querySelector('.var_name').value || '(unknown variable name)';
        code.push('    r = rnd01();');
        code.push('    X_new.' + var_name + ' = r * X1.' + var_name + ' + (1.0 - r) * X2.' + var_name + ';');
    });
    code.push('    return X_new;');
    code.push('}');
    code.push('');
    // calculate_SO_total_fitness
    if(multiobj)
    {
        code.push('std::vector<double> calculate_MO_objectives(const GA_Type::thisChromosomeType &X)');
        code.push('{');
        code.push('    return {');
        const array_objs = [];
        objElements.forEach(function(elem) {
            const obj_name = elem.querySelector('.obj_cell_name').value || '(unknown variable name)';
            array_objs.push('        X.middle_costs.' + obj_name);
        });
        code.push(array_objs.join(',\n'));
        code.push('    };');
        code.push('}');
    }
    else
    {
        code.push('double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)');
        code.push('{');
        code.push('    // finalize the cost');
        code.push('    double final_cost = 0.0;');
        objElements.forEach(function(elem) {
            const obj_name = elem.querySelector('.obj_cell_name').value || '(unknown variable name)';
            code.push('    final_cost += X.middle_costs.' + obj_name + ';');
        });
        code.push('    return final_cost;');
        code.push('}');
    }

    code.push('');
    code.push('std::ofstream output_file;');
    code.push('');
    // Report Generation
    if(multiobj)
    {
        code.push('void MO_report_generation(');
        code.push('    int generation_number,');
        code.push('    const EA::GenerationType<' + solution_name + ', ' + eval_name + '> &last_generation,');
        code.push('    const std::vector<unsigned int>& pareto_front)');
        code.push('{');
        code.push('    (void) last_generation;');
        code.push('');
        code.push('    cout << "Generation [" << generation_number << "], ";');
        code.push('    cout << "Pareto-Front {";');
        code.push('    for(unsigned int i = 0; i < pareto_front.size(); i++)');
        code.push('    {');
        code.push('        cout << (i > 0 ? "," : "");');
        code.push('        cout << pareto_front[i];');
        code.push('    }');
        code.push('    cout << "}" << endl;');
        code.push('}');
        code.push('');
        code.push('void save_results(const GA_Type &ga_obj)');
        code.push('{');
        code.push('    std::ofstream output_file;');
        code.push('    output_file.open("paretofront.txt");');
        code.push('    output_file');
        code.push('        << "N"');
        objElements.forEach(function(elem) {
            const obj_name = elem.querySelector('.obj_cell_name').value || '(unknown variable name)';
            code.push('        << "\\t" << "' + obj_name + '"');
        });
        code.push('        << "\\t" << "solution" << "\\n";');
        code.push('    std::vector<unsigned int> paretofront_indices = ga_obj.last_generation.fronts[0];');
        code.push('    for(unsigned int i : paretofront_indices)');
        code.push('    {');
        code.push('        const auto &X = ga_obj.last_generation.chromosomes[i];');
        code.push('        output_file');
        code.push('            << i << "\\t"');
        objElements.forEach(function(elem) {
            const obj_name = elem.querySelector('.obj_cell_name').value || '(unknown variable name)';
            code.push('            << X.middle_costs.' + obj_name + ' << "\\t"');
        });
        code.push('            << X.genes.to_string() << "\\n";');
        code.push('    }');
        code.push('    output_file.close();');
        code.push('}');
    }
    else
    {
        code.push('void SO_report_generation(');
        code.push('    int generation_number,');
        code.push('    const EA::GenerationType<' + solution_name + ', ' + eval_name + '> &last_generation,');
        code.push('    const ' + solution_name + '& best_genes)');
        code.push('{');
        code.push('    cout');
        code.push('        << "Generation [" << generation_number << "], "');
        code.push('        << "Best=" << last_generation.best_total_cost << ", "');
        code.push('        << "Average=" << last_generation.average_cost << ", "');
        code.push('        << "Best genes=(" << best_genes.to_string() << ")" << ", "');
        code.push('        << "Exe_time=" << last_generation.exe_time');
        code.push('        << endl;');
        code.push('');
        code.push('    output_file');
        code.push('        << generation_number << "\\t"');
        code.push('        << last_generation.average_cost << "\\t"');
        code.push('        << last_generation.best_total_cost << "\\t"');
        code.push('        << best_genes.to_string() << "\\n";');
        code.push('}');
    }
    code.push('');
    // main.cpp
    code.push('int main()');
    code.push('{');
    if(!multiobj)
    {
        code.push('    output_file.open("results.txt");');
        code.push('    output_file << "step" << "\\t" << "cost_avg" << "\\t" << "cost_best" << "\\t" << "solution_best" << "\\n";');
        code.push('');
    }
    code.push('    EA::Chronometer timer;');
    code.push('    timer.tic();');
    code.push('');
    code.push('    GA_Type ga_obj;');
    if(multiobj)
        code.push('    ga_obj.problem_mode = EA::GA_MODE::NSGA_III;');
    else
        code.push('    ga_obj.problem_mode = EA::GA_MODE::SOGA;');
    const heavyness = getSelectValue('heavyness');
    const multi_threading = (heavyness !== 'light');
    code.push('    ga_obj.multi_threading = ' + multi_threading + ';');
    if(heavyness === 'medium')
        code.push('    ga_obj.idle_delay_us = 1; // switch between threads quickly');
    if(heavyness === 'heavy')
        code.push('    ga_obj.idle_delay_us = 10; // switch between threads quickly');
    if(heavyness !== 'light')
    {
        const dynamic_threading = (heavyness === 'heavy');
        code.push('    ga_obj.dynamic_threading = ' + dynamic_threading + ';');
    }
    const verbose = (getSelectValue('verbose') === "true");
    code.push('    ga_obj.verbose = ' + verbose + ';');
    var population;
    switch(getSelectValue('population'))
    {
        case "small":
            population = 50;
            break;
        case "medium":
            population = 200;
            break;
        case "large":
            population = 1000;
            break;
        default:
            population = 200;
    }
    if(multiobj && population < 100)
        population = 100;
    code.push('    ga_obj.population = ' + population.toString() + ';');
    code.push('    ga_obj.generation_max = 1000;');
    if(multiobj)
        code.push('    ga_obj.calculate_MO_objectives = calculate_MO_objectives;');
    else
        code.push('    ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;');
    code.push('    ga_obj.init_genes = init_genes;');
    code.push('    ga_obj.eval_solution = eval_solution;');
    code.push('    ga_obj.mutate = mutate;');
    code.push('    ga_obj.crossover = crossover;');
    if(multiobj)
        code.push('    ga_obj.MO_report_generation = MO_report_generation;');
    else
        code.push('    ga_obj.SO_report_generation = SO_report_generation;');
    code.push('    ga_obj.crossover_fraction = 0.7;');
    code.push('    ga_obj.mutation_rate = 0.2;');
    if(!multiobj)
    {
        code.push('    ga_obj.best_stall_max = 10;');
        code.push('    ga_obj.elite_count = 10;');
    }
    code.push('    ga_obj.solve();');
    code.push('');
    code.push('    cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;');
    code.push('');
    if(!multiobj)
        code.push('    output_file.close();');
    else
        code.push('    save_results(ga_obj);');
    code.push('    return 0;');
    code.push('}');
    code.push('');

    const codeElement = document.getElementById('code');
    codeElement.innerHTML = '<pre>' + escapeHTML(code.join('\n')) + '</pre>';
}
