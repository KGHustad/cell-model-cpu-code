gotran2c TP06.ode --code.body.use_enum=1 --functions.monitored.generate=0 --functions.rhs.generate=0 --solvers.hybrid_generalized_rush_larsen.stiff_states=m,h,j,d,f,f2,fCass,Xr1,Xr2,Xs,s,r --solvers.hybrid_generalized_rush_larsen.generate=1 --solvers.hybrid_generalized_rush_larsen.function_name=step_RL_singlecell --output=TP06_RL.c

