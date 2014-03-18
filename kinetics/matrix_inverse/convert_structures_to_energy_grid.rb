require File.join(File.dirname(__FILE__), "structural_enumeration_to_probability_transition_matrix.rb")

raise ArgumentError.new("ruby convert_structures_to_energy_grid.rb [input_file] [hastings / no_hastings]") unless ARGV.size == 2

move_klass = case (algorithm = ARGV[1])
when "hastings"    then MoveWithHastings
when "no_hastings" then MoveWithoutHastings
else raise ArgumentError.new("Second argument must be one of hastings / no_hastings") end

Converter.to_matrix(ARGV[0], move_klass, ARGV[0])
