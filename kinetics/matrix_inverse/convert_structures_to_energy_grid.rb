require "set"

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
      get_pairings(structure).each_with_index.inject(Set.new) do |set, (j, i)|
        j >= 0 ? set << Set[i, j] : set
      end
    end

    def get_pairings(structure)
    	stack = []
  
      structure.each_char.each_with_index.inject(Array.new(structure.length, -1)) do |array, (symbol, index)|
    	  array.tap do      
    	    case symbol
    	    when "(" then stack.push(index)
    	    when ")" then 
    	      if stack.empty?
    	        raise "Too many ')' in '#{structure}'"
    	      else
    	        stack.pop.tap do |opening|
    	          array[opening] = index
    	          array[index]   = opening
    	        end
    	      end
    	    end
    	  end
    	end.tap do
    	  raise "Too many '(' in '#{structure}'" unless stack.empty?
    	end
    end
    
    def bp_distance(structure_1, structure_2)
      # Takes two structures and calculates the distance between them by |symmetric difference(bp_in_a, bp_in_b)|
      raise "The two structures are not the same length" unless structure_1.length == structure_2.length
  
      bp_set_1, bp_set_2 = base_pairs(structure_1), base_pairs(structure_2)
  
      ((bp_set_1 - bp_set_2) + (bp_set_2 - bp_set_1)).count
    end
  end
end

class MetropolisMove
  attr_reader :from, :to, :p
  
  def initialize(from, to, p: nil, num_moves: nil)
    @from, @to   = from, to
    @p           = p ? p : probability(num_moves)
  end
  
  def probability(num_moves)
    [1.0, Math.exp(-(to.mfe - from.mfe) / Structure::RT)].min / num_moves
  end
  
  def to_csv
    "%d,%d,%.#{Float::DIG}f" % [from.index, to.index, p]
  end
  
  def to_debug_string
    "from: %s (%+.2f)\tto: %s (%+.2f)\tp: %.#{Float::DIG}f" % [from.structure, from.mfe, to.structure, to.mfe, p]
  end
  
  def <=>(other_move)
    from.index != other_move.from.index ? from.index <=> other_move.from.index : to.index <=> other_move.to.index
  end
end

structures  = File.read(ARGV[0]).split(?\n).map { |line| line.split(?\t) }.each_with_index.map(&Structure.method(:new))
empty_index = structures.find { |structure| structure.mfe.zero? }.index
mfe_index   = structures.min { |a, b| a.mfe <=> b.mfe }.index
move_set    = structures.inject({}) do |hash, structure_1|
  hash.tap do
    hash[structure_1] = structures.select do |structure_2|
      structure_1.distance(structure_2) == 1
    end
  end
end

move_list = move_set.inject([]) do |list, (from, to_array)|
  outgoing  = to_array.map { |to| MetropolisMove.new(from, to, num_moves: to_array.size) }
  all_moves = (outgoing + [MetropolisMove.new(from, from, p: 1 - outgoing.map(&:p).inject(&:+))]).select { |move| move.p > 0 }
  list + all_moves 
end.sort

File.open("%s__%d_%d.csv" % [File.basename(ARGV[0], ".txt"), empty_index, mfe_index], ?w) do |file|
  file.write(move_list.map(&:to_csv).join(?\n) + ?\n)
end