.PHONY: all clean

JULIA := julia

all:
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=. main.jl

clean:
	rm -rf results/*