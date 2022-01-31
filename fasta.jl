using FASTX
using BioSequences

function fasta_read(fasta_input)
    reader = open(FASTA.Reader, fasta_input)
    num = 1
    for record in reader
        sequence = FASTA.sequence(record)
        id = FASTA.identifier(record)
        (description(record) !== nothing ) && (id = join([id, description(record)], " "))
        #println(id)
        names = split(id, ";")
        if length(names) < 3
            println(">$id", names[1], " Accession-$num;")
        else
            println(">$id", names[end-2], " ", names[end-1], " Accession-$num;")
        end
        println(sequence)

        num +=1
    end
    close(reader)
end

fasta_read(ARGS[1])
