# Example
# GAGCCTACTAACGGGAT
# CATCGTAATGACGGCCT
# ^ ^ ^  ^ ^    ^^

function hamming_distance(a, b)
    if length(a) === length(b)
        d = 0
        for (index, x) in enumerate(a)
            (x !== b[index]) && (d += 1)
        end
        return d
    else
        throw("The sequenes are not equal in length!")
    end
end

hamming_distance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")