CXX=g++
#OMPFlags=-fopenmp

# Her skrives de forskellige "path's" som makefilen kan bruge til at finde de forskellige pakker. 
# Har f.eks. udkommenteret de "rigtige" paths, men man kan også finde pathen ved at bruge gsl-config.  
CFlags=-c -Wall -O3 -pedantic -MMD -std=c++11 -Werror
CFlags+=-ffast-math
CFlags+=$(OMPFlags)
CFlags+=`gsl-config --cflags`
#CFlags+=-I/opt/homebrew/Cellar/gsl/2.6/include
#GSLFlags=-L/opt/homebrew/Cellar/gsl/2.6/lib -lgsl -lgslcblas
GSLFlags=`gsl-config --libs`

# Her defineres de filer vi gerne vil have executed:
Sources=$(wildcard src/*.cpp)

# Her defineres hvordan vi kommer den til Include pakken: 
IncludeDir=-I./include


# Vi executer main filen: 
Executables=main

# Når man executer en fil, skaber compileren en object fil (.o), som indeholder den binære representation af koden. Denne køres så direkte af 
# computeren, som så giver os vores resultater. 
AllObjects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))


# Her køres alt shitet, men her skal jeg lige have Mauricios hjælp:
all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data obj fig
	$(CXX) $(Objects) $(addprefix obj/,$@.o) $(GSLFlags) $(OMPFlags) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

test: $(Executables)
	$(foreach exe,$(Executables),./$(exe);)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)

