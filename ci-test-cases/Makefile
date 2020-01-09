# This library is free and distributed under
# Mozilla Public License Version 2.0.

CXX:=g++
CFLAGS_R:=-O3 -s
CFLAGS_D:=-Og -DNDEBUG -D_GLIBCXX_SANITIZE_VECTOR
CLANG_CFLAGS:=-std=c++11 -pthread -I../src -Wall -Wfatal-errors -Wextra
CFLAGS:=-std=c++11 -pthread -I../src -Wall -Wconversion -Werror -Wfatal-errors -Wextra
SANITIZER_A:=-fsanitize=address
SANITIZER_T:=-fsanitize=thread
SANITIZER_UB:=-fsanitize=undefined
SRC_S:=various-so-rastrigin.cpp
SRC_M:=various-mo-dtlz2.cpp
OUT:=./a.out

all:
	$(error only the designated CI targets should be called)

so-skim:
	$(CXX) $(CFLAGS_R) $(CFLAGS) $(SRC_S) -o $(OUT)
	$(OUT)

mo-skim:
	$(CXX) $(CFLAGS_R) $(CFLAGS) $(SRC_M) -o $(OUT)
	$(OUT)

so-sanitizer-address:
	$(CXX) $(CFLAGS_D) $(CFLAGS) $(SANITIZER_A) $(SRC_S) -o $(OUT) $(SANITIZER_A)
	env ASAN_OPTIONS="detect_stack_use_after_return=1" $(OUT)

mo-sanitizer-address:
	$(CXX) $(CFLAGS_D) $(CFLAGS) $(SANITIZER_A) $(SRC_M) -o $(OUT) $(SANITIZER_A)
	env ASAN_OPTIONS="detect_stack_use_after_return=1" $(OUT)

so-sanitizer-thread:
	$(CXX) $(CFLAGS_D) $(CFLAGS) $(SANITIZER_T) $(SRC_S) -o $(OUT) $(SANITIZER_T)
	$(OUT)

mo-sanitizer-thread:
	$(CXX) $(CFLAGS_D) $(CFLAGS) $(SANITIZER_T) $(SRC_M) -o $(OUT) $(SANITIZER_T)
	$(OUT)

so-sanitizer-ub:
	$(CXX) $(CFLAGS_D) $(CFLAGS) $(SANITIZER_UB) $(SRC_S) -o $(OUT) $(SANITIZER_UB) -fuse-ld=gold
	$(OUT)

mo-sanitizer-ub:
	$(CXX) $(CFLAGS_D) $(CFLAGS) $(SANITIZER_UB) $(SRC_M) -o $(OUT) $(SANITIZER_UB) -fuse-ld=gold
	$(OUT)

so-clang:
	clang++ $(CFLAGS_R) $(CLANG_CFLAGS) $(SRC_S) -o $(OUT)
	$(OUT)

mo-clang:
	clang++ $(CFLAGS_R) $(CLANG_CFLAGS) $(SRC_M) -o $(OUT)
	$(OUT)
