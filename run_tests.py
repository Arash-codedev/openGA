from pathlib import Path
import platform
import subprocess
import sys




# Colors
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
BLUE = '\033[0;34m'
NC = '\033[0m'

# Global test definitions
TESTS = [
    {
        'name': 'Example: single-objective',
        'sources': ['examples/so-1/example_so1.cpp'],
        'output': 'example_so1'
    },
    {
        'name': 'Example: single-objective Rastrigin',
        'sources': ['examples/so-rastrigin/so-rastrigin.cpp'],
        'output': 'so-rastrigin'
    },
    {
        'name': 'Example: single-objective bind',
        'sources': ['examples/so-bind/example_bind.cpp'],
        'output': 'example_bind'
    },
    {
        'name': 'Example: single-objective with multiple linked sources',
        'sources': [
            'examples/so-multi-source/multi-source-part1.cpp',
            'examples/so-multi-source/multi-source-part2.cpp',
            'examples/so-multi-source/multi-source-part3.cpp'
        ],
        'output': 'example_mulisource'
    },
    {
        'name': 'Example: multi-objective',
        'sources': ['examples/mo-1/example_mo1.cpp'],
        'output': 'example_mo1'
    },
    {
        'name': 'Example: multi-objective DTLZ2',
        'sources': ['examples/mo-dtlz2/mo-dtlz2.cpp'],
        'output': 'mo-dtlz2'
    },
    {
        'name': 'Example: GA with solution initialization',
        'sources': ['examples/so-init-solutions/example_so-init-solutions.cpp'],
        'output': 'example_so-init-solutions'
    }
]

# Sanitizer configurations
SANITIZERS = [
    {'name': 'normal', 'flags': []},
    {'name': 'address', 'flags': ['-fsanitize=address']},
    {'name': 'thread', 'flags': ['-fsanitize=thread']},
    {'name': 'undefined', 'flags': ['-fsanitize=undefined']}
]

# Base compiler flags
BASE_FLAGS = [
    '-std=c++11',
    '-pthread',
    '-Isrc',
    '-Wall',
    '-Wconversion',
    '-Werror',
    '-Wfatal-errors',
    '-Wextra'
]




def get_command(test, os_type, sanitizer):
    sources = test['sources']
    output = test['output']
    Path('bin').mkdir(exist_ok=True)

    flags = BASE_FLAGS + sanitizer['flags']

    if os_type == 'windows':
        if sanitizer['name'] == 'normal':
            output_path = f'bin\\{output}.exe'
        else:
            output_path = f'bin\\{output}_{sanitizer["name"]}.exe'
        return ['cl', '/EHsc'] + flags + [f'/Fe:{output_path}'] + sources
    elif os_type == 'darwin':
        if sanitizer['name'] == 'normal':
            output_path = f'bin/{output}'
        else:
            output_path = f'bin/{output}_{sanitizer["name"]}'
        return ['clang++'] + flags + ['-o', output_path] + sources
    else:
        if sanitizer['name'] == 'normal':
            output_path = f'bin/{output}'
        else:
            output_path = f'bin/{output}_{sanitizer["name"]}'
        return ['g++'] + flags + ['-o', output_path] + sources



def run_test(test, os_type, test_num, sanitizer):
    command = get_command(test, os_type, sanitizer)
    print(f'{YELLOW}[{test_num}] Running: {test["name"]} [{sanitizer["name"]}]{NC}')
    print(f'Command: {" ".join(command)}')

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f'{GREEN}✓ PASSED: {test["name"]} [{sanitizer["name"]}]{NC}')
        return True
    except subprocess.CalledProcessError as e:
        print(f'{RED}✗ FAILED: {test["name"]} [{sanitizer["name"]}]{NC}')
        if e.stderr:
            print(f'{RED}Error: {e.stderr}{NC}')
        return False



def main():
    os_type = platform.system().lower()
    print(f'Detected OS: {os_type}')
    print('=' * 41)

    total_tests = 0
    passed_tests = 0
    failed_tests = 0

    for test in TESTS:
        print(f'\n{BLUE}Testing: {test["name"]}{NC}')
        print('-' * 41)

        for sanitizer in SANITIZERS:
            total_tests += 1
            if run_test(test, os_type, total_tests, sanitizer):
                passed_tests += 1
            else:
                failed_tests += 1

    print()
    print('=' * 41)
    print(f'Total: {total_tests} | {GREEN}Passed: {passed_tests}{NC} | {RED}Failed: {failed_tests}{NC}')
    print('=' * 41)

    sys.exit(0 if failed_tests == 0 else 1)




if __name__ == '__main__':
    main()
