function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.initial_condition = data(2);
parms.W_f = data(3);
parms.W_i = data(4);
parms.beta = data(5);
parms.x0 = data(6);
parms.y0 = data(7)

end