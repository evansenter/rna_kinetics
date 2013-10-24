require "set"

raise ArgumentError.new("ruby convert_structures_to_energy_grid.rb [input_file] [hastings / no_hastings]") unless ARGV.size == 2

class Structure
  RT = 1e-3 * 1.9872041 * (273.15 + 37) # kcal / K / mol @ 37C
  
  attr_reader :structure, :mfe, :index
  
  def initialize((structure, mfe), index)
    @structure, @mfe, @index = structure, mfe.to_f, index
  end
  
  def distance(other_str)
    self.class.bp_distance(structure, other_str.structure)
  end
  
  class << self
    def base_pairs(structure)
    	stack = []
  
      structure.each_char.each_with_index.inject(Set.new) do |set, (symbol, index)|
    	  set.tap do      
    	    case symbol
    	    when ?( then stack.push(index)
    	    when ?) then 
            set << [stack.pop, index]
    	    end
    	  end
    	end
    end
    
    def bp_distance(structure_1, structure_2)
      bp_set_1, bp_set_2 = base_pairs(structure_1), base_pairs(structure_2)
  
      ((bp_set_1 - bp_set_2) + (bp_set_2 - bp_set_1)).count
    end
  end
end

class Move
  attr_reader :from, :to, :move_set, :p
  
  def initialize(from, to, move_set, p: nil)
    @from, @to, @move_set = from, to, move_set
    @p                    = p ? p : probability
  end
  
  def to_csv
    "%d,%d,%.#{Float::DIG}f" % [from.index, to.index, p]
  end
  
  def to_debug
    "from: %s (%+.2f)\tto: %s (%+.2f)\tp: %.#{Float::DIG}f" % [from.structure, from.mfe, to.structure, to.mfe, p]
  end
  
  def <=>(other_move)
    from.index != other_move.from.index ? from.index <=> other_move.from.index : to.index <=> other_move.to.index
  end
end

class MoveWithoutHastings < Move
  def probability
    [1.0, Math.exp(-(to.mfe - from.mfe) / Structure::RT)].min / move_set[from].size
  end
end

class MoveWithHastings < Move
  def probability
    [1.0, (move_set[from].size.to_f / move_set[to].size.to_f) * (Math.exp(-(to.mfe - from.mfe) / Structure::RT))].min / move_set[from].size
  end
end

move_klass = case (algorithm = ARGV[1])
when "hastings"    then MoveWithHastings
when "no_hastings" then MoveWithoutHastings
else raise ArgumentError.new("Second argument must be one of hastings / no_hastings") end

structures  = File.read(ARGV[0]).split(?\n).map { |line| line.split(?\t) }.each_with_index.map(&Structure.method(:new))
empty_index = structures.find { |structure| structure.mfe.zero? }.index
mfe_index   = structures.min { |a, b| a.mfe <=> b.mfe }.index
move_set    = structures.inject({}) do |hash, structure_1|
  hash.tap do
    hash[structure_1] = structures.select do |structure_2|
      (structure_1.structure.scan(?.).size - structure_2.structure.scan(?.).size).abs == 2 && structure_1.distance(structure_2) == 1
    end
  end
end

move_list = move_set.inject([]) do |list, (from, to_array)|
  outgoing  = to_array.map { |to| move_klass.new(from, to, move_set) }
  all_moves = (outgoing + [move_klass.new(from, from, move_set, p: 1 - outgoing.map(&:p).inject(&:+))]).select { |move| move.p > 0 }
  list + all_moves 
end.sort

File.open("%s__%d_%d_%s.csv" % [File.basename(ARGV[0], ".txt"), empty_index, mfe_index, algorithm], ?w) do |file|
  file.write(move_list.map(&:to_csv).join(?\n) + ?\n)
end
