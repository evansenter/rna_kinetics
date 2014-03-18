require "set"
require "active_support/inflector"

module Enumerable
  def inner_inject(default = :not_used, &block)
    map { |object| default == :not_used ? object.inject(&block) : object.inject(default, &block) }
  end
end

class Structure
  RT = 1e-3 * 1.9872041 * (273.15 + 37) # kcal / K / mol @ 37C

  attr_reader :structure, :base_pairs, :unpaired_count, :mfe, :index

  def initialize((structure, mfe), index)
    @structure, @base_pairs, @unpaired_count, @mfe, @index = structure, self.class.bp_set(structure), structure.scan(?.).size, mfe.to_f, index
  end

  def distance(other_str)
    self.class.bp_distance(base_pairs, other_str.base_pairs)
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

    def bp_set(structure)
      base_pairs(structure).map { |a| [2, 3].zip(a).inner_inject(&:**).inject(&:*) }.to_set
    end

    def bp_distance(bp_set_1, bp_set_2)
      ((bp_set_1 - bp_set_2) + (bp_set_2 - bp_set_1)).count
    end
  end
end

class Move
  attr_reader :from, :to, :from_size, :to_size, :p

  def initialize(from, to, from_size, to_size, p: nil)
    @from, @to, @from_size, @to_size = from, to, from_size.to_f, to_size.to_f
    @p                               = p ? p : probability
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
    # return MIN(1., exp(-(klp_matrix.p[j] - klp_matrix.p[i]) / RT)) / number_of_adjacent_moves[i];
    [1.0, Math.exp(-(to.mfe - from.mfe) / Structure::RT)].min / from_size
  end
end

class MoveWithHastings < Move
  def probability
    # return MIN(1., (number_of_adjacent_moves[i] / number_of_adjacent_moves[j]) * exp(-(klp_matrix.p[j] - klp_matrix.p[i]) / RT)) / number_of_adjacent_moves[i];
    [1.0, (from_size / to_size) * (Math.exp(-(to.mfe - from.mfe) / Structure::RT))].min / from_size
  end
end

module Converter
  class << self
    def to_matrix(structures, move_klass = MoveWithHastings, basename = false)
      structures = parse_structures(structures)
      moves      = move_set(structures)
      matrix     = move_list(moves, move_klass)

      if basename
        writer(structures, matrix, move_klass, basename)
      else
        matrix.tap do
          puts matrix.map(&:to_debug).join(?\n)
        end
      end
    end

    def parse_structures(structures)
      if structures.is_a?(String) && File.exists?(structures)
        File.read(structures).split(?\n).map { |line| line.split(?\t) }.each_with_index.map(&Structure.method(:new))
      else
        structures.each_with_index.map { |subopt, index| Structure.new([subopt.structure, subopt.mfe], index) }
      end
    end

    def move_set(structures)
      structures.inject({}) do |hash, structure_1|
        hash.tap do
          hash[structure_1] = structures.select do |structure_2|
            (structure_1.unpaired_count - structure_2.unpaired_count).abs == 2 && structure_1.distance(structure_2) == 1
          end
        end
      end
    end

    def move_list(moves, move_klass)
      moves.inject([]) do |list, (from, to_array)|
        outgoing  = to_array.map { |to| move_klass.new(from, to, moves[from].size, moves[to].size) }
        all_moves = (outgoing + [move_klass.new(from, from, moves[from].size, moves[from].size, p: 1 - outgoing.map(&:p).inject(&:+))]).select { |move| move.p > 0 }
        list + all_moves
      end.sort
    end

    def writer(structures, matrix, move_klass, basename)
      File.open("%s__%d_%d_%s.csv" % [
        File.basename(basename, ".txt"),
        structures.find { |structure| structure.mfe.zero? }.index,
        structures.min { |a, b| a.mfe <=> b.mfe }.index,
        move_klass.name.underscore
      ], ?w) do |file|
        file.write(matrix.map(&:to_csv).join(?\n) + ?\n)
      end
    end
  end
end
