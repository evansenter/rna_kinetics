dyn.load("rnaeval_r.so")
print('.C("energy_of_struct_r", sequence = "GGGGGCCCCC", structure = "..........", energy = as.double(1))')
.C("energy_of_struct_r", sequence = "GGGGGCCCCC", structure = "..........", energy = as.double(1))$energy
print('.C("energy_of_struct_r", sequence = "GGGGGCCCCC", structure = "(((....)))", energy = as.double(1))')
.C("energy_of_struct_r", sequence = "GGGGGCCCCC", structure = "(((....)))", energy = as.double(1))$energy