
models = $(wildcard *.ode)
models_c = $(patsubst %.ode,%.c,$(models))

#$(info models: $(models))
#$(info models_c: $(models_c))

gotran2c_args = --solvers.explicit_euler.generate=1 --solvers.explicit_euler.function_name=step_FE --solvers.generalized_rush_larsen.generate=1 --solvers.generalized_rush_larsen.function_name=step_GRL1 --code.body.use_enum=1
gotran2c_args += --solvers.explicit_euler.function_name=step_FE_singlecell
gotran2c_args += --solvers.generalized_rush_larsen.function_name=step_GRL1_singlecell
gotran2c_args += --functions.monitored.generate=0 --functions.rhs.generate=0

all: $(models_c)

%.c : %.ode
	gotran2c $< $(gotran2c_args) --output=$@

